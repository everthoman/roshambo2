"""Consensus pharmacophore from several dopaminergic agonists.

roshambo2 only aligns pairwise, so a multi-ligand model is built on top:

  1. pick the most rigid ligand as the reference (its lowest-MMFF-energy conformer)
  2. roshambo2-align every other ligand's conformers to that reference, keep the
     best-scoring pose (ProjectedPointPharmacophoreGenerator colour model)
  3. run the pharmacophore generator on every aligned ligand
  4. cluster same-family feature points across all ligands; a cluster seen in
     >= MIN_LIGANDS of them becomes a consensus point at its centroid

Ligands are read from dopamine_ligands.txt (one "name SMILES" per line, neutral);
the basic amine of each is protonated before use.

Outputs (cwd):
  consensus_model.sdf              the consensus points
  aligned_<ligand>.sdf             each ligand in the reference frame (polar H)
  consensus_pharmacophore.pml/.png PyMOL loader + render
"""
import os
import sys
from collections import defaultdict

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolDescriptors
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator
from ligtools import pose_with_H, polar_only
import overlaptools as ot

LIGAND_FILE = os.path.join(os.path.dirname(__file__), "dopamine_ligands.txt")
N_CONFS = 50
CLUSTER_RADIUS = 1.6                                   # Angstrom

# trivalent neutral N, not amide / amidine / N-oxide / N-N / aniline
_BASIC_N = Chem.MolFromSmarts(
    "[#7;X3;+0;!$([#7][#6]=[O,N,S]);!$([#7]=*);!$([#7]~[O,N,P,S]);!$([#7]a)]")


def protonate(smiles):
    m = Chem.MolFromSmiles(smiles)
    rw = Chem.RWMol(m)
    for (i,) in m.GetSubstructMatches(_BASIC_N):
        a = rw.GetAtomWithIdx(i)
        a.SetFormalCharge(1)
        a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
        a.SetNoImplicit(True)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return m


LIGANDS = {}
for line in open(LIGAND_FILE):
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    name, smi = line.split(None, 1)
    LIGANDS[name] = Chem.MolToSmiles(protonate(smi))

N = len(LIGANDS)
MIN_LIGANDS = int(sys.argv[1]) if len(sys.argv) > 1 else max(2, round(0.6 * N))

cg = ProjectedPointPharmacophoreGenerator()
idx2fam = cg.get_index_to_feature()


def prep(smiles, name):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(m, numConfs=N_CONFS, params=p)
    ff = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(m)
    m.SetProp("_Name", name)
    return m, [e for _, e in ff]


mols = {name: prep(smi, name) for name, smi in LIGANDS.items()}

ref_name = min(LIGANDS, key=lambda n: rdMolDescriptors.CalcNumRotatableBonds(
    Chem.MolFromSmiles(LIGANDS[n])))
ref_mol, ref_e = mols[ref_name]
ref_ci = int(np.argmin(ref_e))
ref_single = Chem.Mol(ref_mol, confId=ref_ci)
ref_single.SetProp("_Name", ref_name)
print(f"reference = {ref_name} (conf {ref_ci});  consensus needs >= {MIN_LIGANDS}/{N} ligands\n")

# ---- 1-2. align every ligand into the reference frame -----------------------
aligned = {ref_name: ref_single}
for name, (m, _) in mols.items():
    if name == ref_name:
        continue
    calc = Roshambo2(ref_single, [m], color=True,
                     remove_Hs_before_color_assignment=False, color_generator=cg)
    row = next(iter(calc.compute(backend="cuda", reduce_over_conformers=False,
                                 optim_mode="combination", combination_param=0.5,
                                 write_scores=False).values())).iloc[0]
    hit_noH = calc.get_best_fit_structures(top_n=1)[f"{ref_name}_0"][0]
    k = int(hit_noH.GetProp("name").rsplit("_", 1)[-1])
    aligned[name] = pose_with_H(hit_noH, Chem.Mol(m, confId=k))
    print(f"  aligned {name:12}  shape={row['tanimoto_shape']:.2f} "
          f"color={row['tanimoto_color']:.2f} combo={row['tanimoto_combination']:.2f}")

# ---- 3. feature points for every aligned ligand ---------------------------
by_family = defaultdict(list)                         # family -> [(ligand, xyz), ...]
for name, mol in aligned.items():
    coords, types = cg.generate_color_atoms(mol)
    for xyz, t in zip(coords, types):
        by_family[idx2fam[int(t)]].append((name, np.asarray(xyz, float)))

# ---- 4. cluster within each family, keep the consensus ones ---------------
consensus = []
for fam, pts in by_family.items():
    clusters = []
    for lig, p in pts:
        for c in clusters:
            if any(np.linalg.norm(p - q) <= CLUSTER_RADIUS for _, q in c):
                c.append((lig, p))
                break
        else:
            clusters.append([(lig, p)])
    for c in clusters:
        ligs = sorted({lig for lig, _ in c})
        if len(ligs) >= MIN_LIGANDS:
            consensus.append(dict(feature=fam, symbol=ot.FEATURE_TO_SYMBOL[fam], d=0.0,
                                  xyz=np.mean([p for _, p in c], axis=0),
                                  n=len(ligs), ligs=ligs))

consensus.sort(key=lambda x: (-x["n"], x["feature"]))
print(f"\n{'feature':<16}{'#ligands':>9}   ligands")
for p in consensus:
    print(f"{p['feature']:<16}{p['n']:>9}   {', '.join(p['ligs'])}")
print(f"\n{len(consensus)} consensus points")

kept = {p["feature"] for p in consensus}
missed = sorted(f for f in by_family if f not in kept)
if missed:
    print(f"below threshold (scattered after alignment): {', '.join(missed)}")

ot.write_model(consensus, "consensus_model.sdf", name="dopaminergic_consensus_pharmacophore")
_m = next(iter(Chem.SDMolSupplier("consensus_model.sdf", sanitize=False)))
_m.SetProp("ligand_counts", ", ".join(f"{p['feature']}={p['n']}" for p in consensus))
with Chem.SDWriter("consensus_model.sdf") as w:
    w.write(_m)

# ---- write aligned ligands + a PyMOL loader ------------------------------
lig_colors = ["cyan", "salmon", "palegreen", "lightblue", "wheat", "violet",
              "yellow", "lightpink", "aquamarine", "paleyellow"]
fnames = {}
for name, mol in aligned.items():
    fn = f"aligned_{name}.sdf"
    fnames[name] = fn
    with Chem.SDWriter(fn) as w:
        w.write(polar_only(mol))

pml = ["# dopaminergic consensus pharmacophore (roshambo2 pairwise alignment to "
       f"{ref_name}, then feature clustering)"]
for i, name in enumerate(aligned):
    obj = name.replace("-", "_")
    pml += [f"load {fnames[name]}, {obj}",
            f"hide everything, {obj}", f"show sticks, {obj}",
            f"set stick_radius, 0.07, {obj}",
            f"color {lig_colors[i % len(lig_colors)]}, {obj} and elem C"]
pml += [
    "set valence, 0",
    "show spheres, (elem H and neighbor (elem N+O))",
    "set sphere_scale, 0.10, (elem H)",
    "color grey40, elem H",
    "",
    "load consensus_model.sdf, pharm",
    "unbond pharm, pharm",
    "show spheres, pharm",
    "set sphere_scale, 0.45, pharm",
    "set sphere_transparency, 0.1, pharm",
]
present = []
for p in consensus:
    if p["feature"] not in present:
        present.append(p["feature"])
for fam in present:
    e, col, _ = ot.FAMILY[fam]
    pml.append(f"color {col}, pharm and elem {e}")
lab = {e.upper(): ot.FAMILY[f][2] for f, (e, _, _) in ot.FAMILY.items()}
pml += [
    "python",
    f"lab = {lab!r}",
    "seen = set()",
    "for at in cmd.get_model('pharm').atom:",
    "    nm = lab.get(at.symbol.upper())",
    "    if nm and nm not in seen:",
    "        seen.add(nm)",
    "        cmd.label(f'pharm and index {at.index}', repr(nm))",
    "        cmd.set('label_position', (0, 0, 3), f'pharm and index {at.index}')",
    "python end",
    "set label_size, 15", "set label_color, black",
    "bg_color white", "set orthoscopic, 1", "set ray_shadows, 0",
    "orient", "zoom all, 2.5",
]
open("consensus_pharmacophore.pml", "w").write("\n".join(pml) + "\n")
print("\nwrote consensus_model.sdf, aligned_*.sdf and consensus_pharmacophore.pml")
