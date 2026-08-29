"""Build a 3D pharmacophore model from ONLY the pharmacophore ("color") features
that overlap between the query (dopamine) and the roshambo2-aligned hit (PD128907).

roshambo2's color term is a sum of Gaussian overlaps between same-family feature
points (cpp_helper_functions.cpp::volume_color). For the default interaction map
every same-family pair has width a=1, height p=1, point weights w=1, so a single
query/fit feature pair separated by d (Angstrom) contributes

    v(d) = (pi/2)**1.5 * exp(-d**2 / 2)          # (pi/2)**1.5 = 1.9687 = v(0)

We call a pair "overlapping" when its overlap fraction f = exp(-d**2/2) exceeds
--min-overlap (default 0.10, i.e. d <= 2.15 A), match query<->fit features
greedily by nearest distance within each family, and write the matched features
(placed at the pair midpoint) as a standalone pharmacophore model.

Outputs (cwd):
  hits_for_query_dopamine_H+_0.sdf                 query + aligned PD128907
  hits_for_query_dopamine_H+_0_color_features.sdf  ALL feature points (per molecule)
  pharmacophore_overlap_model.sdf                  ONLY the overlapping features
  overlap_pharmacophore.pml                        PyMOL loader for the model
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2

MIN_OVERLAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.10
V0 = (np.pi / 2) ** 1.5

FEATURE_TO_SYMBOL = {
    "Donor": "N", "Acceptor": "O", "PosIonizable": "Fe",
    "NegIonizable": "Cl", "Aromatic": "S", "Hydrophobe": "Br",
}
SYMBOL_TO_FEATURE = {v: k for k, v in FEATURE_TO_SYMBOL.items()}
SYMBOL_TO_Z = {"N": 7, "O": 8, "Fe": 26, "Cl": 17, "S": 16, "Br": 35}


def gen(smiles, name, nconfs):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=nconfs, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol


# ---- 1. overlay -----------------------------------------------------------
query = gen("C1=CC(=C(C=C1CC[NH3+])O)O", "dopamine_H+", 1)
dataset = [gen("CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O", "PD128907_H+", 20)]

calc = Roshambo2(query, dataset, color=True, remove_Hs_before_color_assignment=False)
scores = calc.compute(backend="cuda", reduce_over_conformers=True,
                      optim_mode="combination", combination_param=0.5,
                      write_scores=False)
df = next(iter(scores.values()))
overlap_color = float(df["overlap_color"].iloc[0])
print(df[["name", "tanimoto_shape", "tanimoto_color",
          "tanimoto_combination", "overlap_color"]].to_string(index=False))

calc.write_best_fit_structures(write_color_pseudomols=True, append_query=True,
                               feature_to_symbol_map=FEATURE_TO_SYMBOL)

# split the combined hits SDF into separate single-molecule files so PyMOL
# loads dopamine and PD128907 as distinct objects
_hits = list(Chem.SDMolSupplier("hits_for_query_dopamine_H+_0.sdf",
                                removeHs=False, sanitize=False))
for _m, _fn in zip(_hits, ("ligand_dopamine.sdf", "ligand_PD128907.sdf")):
    with Chem.SDWriter(_fn) as _w:
        _w.write(_m)

# ---- 2. read the two feature clouds (same coordinate frame) ---------------
supp = Chem.SDMolSupplier("hits_for_query_dopamine_H+_0_color_features.sdf",
                          removeHs=False, sanitize=False)
mols = [m for m in supp]
qmol, fmol = mols[0], mols[1]          # entry 0 = query, entry 1 = aligned PD128907


def feats(m):
    conf = m.GetConformer()
    return [(a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx())))
            for a in m.GetAtoms()]


qf, ff = feats(qmol), feats(fmol)

# ---- 3. greedy nearest matching within each family -----------------------
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
        f = np.exp(-d * d / 2.0)
        if f >= MIN_OVERLAP:
            mid = 0.5 * (qf[i][1] + ff[j][1])
            pairs.append(dict(feature=SYMBOL_TO_FEATURE[sym], symbol=sym,
                              d=d, frac=f, v=V0 * f, xyz=mid))

pairs.sort(key=lambda p: -p["v"])
explained = sum(p["v"] for p in pairs)
print(f"\noverlapping features (f >= {MIN_OVERLAP}):")
print(f"{'feature':<14}{'dist(A)':>9}{'overlap_frac':>14}{'v_color':>10}")
for p in pairs:
    print(f"{p['feature']:<14}{p['d']:>9.2f}{p['frac']:>14.3f}{p['v']:>10.3f}")
print(f"\nmatched pairs: {len(pairs)}   "
      f"sum v_color = {explained:.2f} / overlap_color {overlap_color:.2f} "
      f"({100*explained/overlap_color:.0f}% of the color overlap)")

# ---- 4. write the overlap-only pharmacophore model -----------------------
rw = Chem.RWMol()
conf = Chem.Conformer(len(pairs))
for k, p in enumerate(pairs):
    idx = rw.AddAtom(Chem.Atom(SYMBOL_TO_Z[p["symbol"]]))
    conf.SetAtomPosition(idx, p["xyz"].tolist())
out = rw.GetMol()
out.AddConformer(conf)
out.SetProp("_Name", "dopamine_x_PD128907_overlap_pharmacophore")
out.SetProp("min_overlap_fraction", str(MIN_OVERLAP))
out.SetProp("features", ", ".join(f"{p['feature']}(d={p['d']:.2f})" for p in pairs))
with Chem.SDWriter("pharmacophore_overlap_model.sdf") as w:
    w.write(out)
print("\nwrote pharmacophore_overlap_model.sdf")

# ---- 5. PyMOL loader ----------------------------------------------------
with open("overlap_pharmacophore.pml", "w") as fh:
    fh.write(
        "# roshambo2 overlap pharmacophore (only features shared by dopamine & PD128907)\n"
        "load ligand_dopamine.sdf, dopamine\n"
        "load ligand_PD128907.sdf, PD128907\n"
        "hide everything, dopamine or PD128907\n"
        "show sticks, dopamine or PD128907\n"
        "util.cbac dopamine\n"
        "util.cbag PD128907\n\n"
        "load pharmacophore_overlap_model.sdf, pharm\n"
        "show spheres, pharm\n"
        "set sphere_scale, 0.5, pharm\n"
        "set sphere_transparency, 0.3, pharm\n"
        "color slate,     pharm and elem N\n"     # Donor
        "color firebrick, pharm and elem O\n"     # Acceptor
        "color orange,    pharm and elem Fe\n"    # PosIonizable
        "color green,     pharm and elem Cl\n"    # NegIonizable
        "color yellow,    pharm and elem S\n"     # Aromatic
        "color grey60,    pharm and elem Br\n"    # Hydrophobe
        "label pharm, {'N':'Donor','O':'Acceptor','FE':'Pos','CL':'Neg',"
        "'S':'Aromatic','BR':'Hydrophobe'}.get(elem.upper(), elem)\n"
        "set label_size, 18\nbg_color white\nset orthoscopic, 1\nzoom all, 3\n"
    )
print("wrote overlap_pharmacophore.pml")
