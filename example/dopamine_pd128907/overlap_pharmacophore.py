"""Build a 3D pharmacophore model from ONLY the pharmacophore ("color") features
that overlap between dopamine (query) and the roshambo2-aligned PD128907 (hit),
using the directional model (ProjectedPointPharmacophoreGenerator: base points +
projected donor / acceptor lone-pair / ring-normal points).

roshambo2's colour term is a sum of Gaussian overlaps between same-family feature
points (cpp_helper_functions.cpp::volume_color). For the default interaction map
every same-family pair has width a=1, height p=1, point weights w=1, so a single
query/hit feature pair separated by d (Angstrom) contributes

    v(d) = (pi/2)**1.5 * exp(-d**2 / 2)          # (pi/2)**1.5 = 1.9687 = v(0)

A pair is "overlapping" when its overlap fraction f = exp(-d**2/2) exceeds
--min-overlap (argv[1], default 0.15, i.e. d <= ~1.95 A). Features are matched
query<->hit greedily by nearest distance within each family and written at the
pair midpoint as a standalone model.

Outputs (cwd):
  overlap_hits_dopamine_H+_0.sdf                query + aligned PD128907
  overlap_hits_dopamine_H+_0_color_features.sdf all feature points (per molecule)
  ligand_dopamine_ovl.sdf / ligand_PD128907_ovl.sdf     the two ligands, separate
  pharmacophore_overlap_model.sdf               ONLY the overlapping features
  overlap_pharmacophore.pml                     PyMOL loader
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator

MIN_OVERLAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15   # d <= ~1.95 A
V0 = (np.pi / 2) ** 1.5
N_QUERY_CONFS = 25
N_DATASET_CONFS = 100

FEATURE_TO_SYMBOL = {
    "Donor": "N", "Acceptor": "O", "PosIonizable": "Fe", "NegIonizable": "Cl",
    "Aromatic": "S", "Hydrophobe": "Br",
    "DonorProj": "F", "AcceptorProj": "I", "AromaticProj": "B",
    "PosIonizableProj": "P",
}
SYMBOL_TO_FEATURE = {v: k for k, v in FEATURE_TO_SYMBOL.items()}
SYMBOL_TO_Z = {"N": 7, "O": 8, "Fe": 26, "Cl": 17, "S": 16, "Br": 35,
               "F": 9, "I": 53, "B": 5, "P": 15}


def gen(smiles, name, n):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=n, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol


def split(mol):
    out = {}
    for c in mol.GetConformers():
        mi = Chem.Mol(mol, confId=c.GetId())
        mi.SetProp("_Name", f"{mol.GetProp('_Name')}_{c.GetId()}")
        out[mi.GetProp("_Name")] = mi
    return out


# ---- 1. overlay: sample dopamine, keep the best-combo conformer -----------
dopamine_confs = split(gen("C1=CC(=C(C=C1CC[NH3+])O)O", "dopamine_H+", N_QUERY_CONFS))
dataset = [gen("CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O", "PD128907_H+", N_DATASET_CONFS)]

scan = Roshambo2(list(dopamine_confs.values()), dataset, color=True,
                 remove_Hs_before_color_assignment=False,
                 conformers_have_unique_names=True,
                 color_generator=ProjectedPointPharmacophoreGenerator())
scores = scan.compute(backend="cuda", reduce_over_conformers=True,
                      optim_mode="combination", combination_param=0.5, write_scores=False)
best_q, best_row = max(((q.rsplit("_", 1)[0], df.iloc[0]) for q, df in scores.items()),
                       key=lambda x: x[1]["tanimoto_combination"])
print(f"best dopamine conformer: {best_q}   "
      f"shape={best_row['tanimoto_shape']:.3f} color={best_row['tanimoto_color']:.3f} "
      f"combo={best_row['tanimoto_combination']:.3f}")

best_mol = dopamine_confs[best_q]
best_mol.SetProp("_Name", "dopamine_H+")
calc = Roshambo2(best_mol, dataset, color=True, remove_Hs_before_color_assignment=False,
                 color_generator=ProjectedPointPharmacophoreGenerator())
df = next(iter(calc.compute(backend="cuda", reduce_over_conformers=True,
                            optim_mode="combination", combination_param=0.5,
                            write_scores=False).values()))
overlap_color = float(df["overlap_color"].iloc[0])

calc.write_best_fit_structures(hits_sdf_prefix="overlap_hits",
                               write_color_pseudomols=True, append_query=True,
                               feature_to_symbol_map=FEATURE_TO_SYMBOL)

_hits = list(Chem.SDMolSupplier("overlap_hits_dopamine_H+_0.sdf",
                                removeHs=False, sanitize=True))
for _m, _fn in zip(_hits, ("ligand_dopamine_ovl.sdf", "ligand_PD128907_ovl.sdf")):
    polar = [a.GetIdx() for a in _m.GetAtoms() if a.GetAtomicNum() in (7, 8)]
    with Chem.SDWriter(_fn) as _w:
        _w.write(Chem.AddHs(_m, addCoords=True, onlyOnAtoms=polar))   # polar H only

# ---- 2. read the two feature clouds (same coordinate frame) --------------
mols = list(Chem.SDMolSupplier("overlap_hits_dopamine_H+_0_color_features.sdf",
                               removeHs=False, sanitize=False))
qmol, fmol = mols[0], mols[1]


def feats(m):
    conf = m.GetConformer()
    return [(a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx())))
            for a in m.GetAtoms()]


qf, ff = feats(qmol), feats(fmol)

# ---- 3. greedy nearest matching within each family ----------------------
pairs = []
for sym in set(s for s, _ in qf) & set(s for s, _ in ff):
    qi = [i for i, (s, _) in enumerate(qf) if s == sym]
    fi = [j for j, (s, _) in enumerate(ff) if s == sym]
    cand = sorted(((np.linalg.norm(qf[i][1] - ff[j][1]), i, j)
                   for i in qi for j in fi), key=lambda x: x[0])
    used_q, used_f = set(), set()
    for d, i, j in cand:
        if i in used_q or j in used_f:
            continue
        used_q.add(i); used_f.add(j)
        fr = np.exp(-d * d / 2.0)
        if fr >= MIN_OVERLAP:
            pairs.append(dict(feature=SYMBOL_TO_FEATURE[sym], symbol=sym, d=d,
                              frac=fr, v=V0 * fr, xyz=0.5 * (qf[i][1] + ff[j][1])))

pairs.sort(key=lambda p: -p["v"])
explained = sum(p["v"] for p in pairs)
print(f"\noverlapping features (f >= {MIN_OVERLAP}):")
print(f"{'feature':<14}{'dist(A)':>9}{'overlap_frac':>14}{'v_color':>10}")
for p in pairs:
    print(f"{p['feature']:<14}{p['d']:>9.2f}{p['frac']:>14.3f}{p['v']:>10.3f}")
print(f"\nmatched pairs: {len(pairs)}   sum v_color = {explained:.2f} / "
      f"overlap_color {overlap_color:.2f} "
      f"({100 * explained / overlap_color:.0f}% of the colour overlap)")

# ---- 4. write the overlap-only pharmacophore model ---------------------
rw = Chem.RWMol()
conf = Chem.Conformer(len(pairs))
for k, p in enumerate(pairs):
    conf.SetAtomPosition(rw.AddAtom(Chem.Atom(SYMBOL_TO_Z[p["symbol"]])), p["xyz"].tolist())
out = rw.GetMol()
out.AddConformer(conf)
out.SetProp("_Name", "dopamine_x_PD128907_overlap_pharmacophore")
out.SetProp("min_overlap_fraction", str(MIN_OVERLAP))
out.SetProp("features", ", ".join(f"{p['feature']}(d={p['d']:.2f})" for p in pairs))
with Chem.SDWriter("pharmacophore_overlap_model.sdf") as w:
    w.write(out)
print("\nwrote pharmacophore_overlap_model.sdf")

# ---- 5. PyMOL loader --------------------------------------------------
with open("overlap_pharmacophore.pml", "w") as fh:
    fh.write(
        "# roshambo2 overlap pharmacophore - only the directional features\n"
        "# shared by dopamine & PD128907 (base points + projected D/A/ring points)\n"
        "load ligand_dopamine_ovl.sdf, dopamine\n"
        "load ligand_PD128907_ovl.sdf, PD128907\n"
        "hide everything, dopamine or PD128907\n"
        "show sticks, dopamine or PD128907\n"
        "set stick_radius, 0.1, dopamine or PD128907\n"
        "set valence, 0\n"
        "util.cbac dopamine\nutil.cbag PD128907\n"
        "color grey30, (dopamine or PD128907) and elem H\n"
        "show spheres, (dopamine or PD128907) and elem H\n"
        "set sphere_scale, 0.11, (dopamine or PD128907) and elem H\n\n"
        "load pharmacophore_overlap_model.sdf, pharm\n"
        "unbond pharm, pharm\n"
        "show spheres, pharm\n"
        "set sphere_scale, 0.4, pharm\n"
        "set sphere_transparency, 0.15, pharm\n"
        "color slate,     pharm and elem N\n"       # Donor
        "color firebrick, pharm and elem O\n"       # Acceptor
        "color orange,    pharm and elem Fe\n"      # PosIonizable
        "color green,     pharm and elem Cl\n"      # NegIonizable
        "color yellow,    pharm and elem S\n"       # Aromatic
        "color grey60,    pharm and elem Br\n"      # Hydrophobe
        "color deepblue,  pharm and elem F\n"       # DonorProj
        "color hotpink,   pharm and elem I\n"       # AcceptorProj
        "color wheat,     pharm and elem B\n"       # AromaticProj
        "color purple,    pharm and elem P\n"       # PosIonizableProj
        "python\n"
        "seen=set()\n"
        "fams={'N':'Donor','O':'Acceptor','S':'Aromatic','BR':'Hydrophobe',"
        "'F':'D-proj','I':'A-proj','B':'ring-proj','P':'Pos-proj'}\n"
        "off={'Donor':(0,1.6,3),'Acceptor':(0,-2.4,3),'Aromatic':(0,0,3),"
        "'Hydrophobe':(0,2.2,3),'D-proj':(-2.4,0,3),'A-proj':(0,-2.0,3),"
        "'ring-proj':(2.2,1.2,3),'Pos-proj':(2.4,-1.4,3)}\n"
        "for at in cmd.get_model('pharm').atom:\n"
        "    fam=fams.get(at.symbol.upper())\n"
        "    if fam and fam not in seen:\n"
        "        seen.add(fam)\n"
        "        sel=f'pharm and index {at.index}'\n"
        "        cmd.label(sel, repr(fam)); cmd.set('label_position', off[fam], sel)\n"
        "python end\n"
        "set label_size, 15\nset label_color, black\n"
        "bg_color white\nset orthoscopic, 1\nset ray_shadows, 0\norient\nzoom all, 2.5\n"
    )
print("wrote overlap_pharmacophore.pml")
