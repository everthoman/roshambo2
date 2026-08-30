"""Shape/color overlay: dopamine (query) vs PD128907 (dataset), amines protonated.

Outputs (in cwd):
  hits_for_query_dopamine_H+_0.sdf                -> query + aligned PD128907
  hits_for_query_dopamine_H+_0_color_features.sdf -> pseudo-atoms at pharmacophore
                                                     feature centres, same order
  ligand_dopamine.sdf / ligand_PD128907.sdf      -> the two ligands as separate files
"""
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from ligtools import pose_with_H, polar_only

def gen(smiles, name, nconfs):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xf00d
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=nconfs, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol

# amino groups protonated (cationic at physiological pH)
query = gen("C1=CC(=C(C=C1CC[NH3+])O)O", "dopamine_H+", 1)
dataset = [gen("CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O", "PD128907_H+", 20)]

calc = Roshambo2(query, dataset, color=True, remove_Hs_before_color_assignment=False)
scores = calc.compute(backend="cuda", reduce_over_conformers=False,
                      optim_mode="combination", combination_param=0.5,
                      write_scores=False)

for qname, df in scores.items():
    cols = [c for c in ("name", "tanimoto_shape", "tanimoto_color",
                        "tanimoto_combination", "tanimoto_combo_legacy",
                        "overlap_volume", "overlap_color") if c in df.columns]
    print(f"\nquery = {qname}   (best of {len(df)} PD128907 conformers)")
    print(df[cols].head(3).to_string(index=False))

# map each pharmacophore family to a real element so PyMOL colours them distinctly
feature_to_symbol = {
    "Donor":        "N",   # blue
    "Acceptor":     "O",   # red
    "PosIonizable": "Fe",  # orange
    "NegIonizable": "Cl",  # green
    "Aromatic":     "S",   # yellow
    "Hydrophobe":   "Br",  # brown
}
calc.write_best_fit_structures(top_n=1, write_color_pseudomols=True, append_query=True,
                               feature_to_symbol_map=feature_to_symbol)
print("\nfeature -> element:", feature_to_symbol)

# one file per ligand (separate objects in PyMOL), carrying the exact hydrogens
# roshambo2 scored (a fresh AddHs guesses its own O-H rotamer - see ligtools)
hit_noH = calc.get_best_fit_structures(top_n=1)["dopamine_H+_0"][0]
k = int(hit_noH.GetProp("name").rsplit("_", 1)[-1])
with Chem.SDWriter("ligand_dopamine.sdf") as _w:
    _w.write(polar_only(query))
with Chem.SDWriter("ligand_PD128907.sdf") as _w:
    _w.write(polar_only(pose_with_H(hit_noH, Chem.Mol(dataset[0], confId=k))))
