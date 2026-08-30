"""Minimal roshambo2 overlay: protonated dopamine (query) vs PD128907 (dataset).

Just runs the shape+colour overlay and prints the score. For the pharmacophore
models and PyMOL renders see overlap_pharmacophore.py.
"""
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2


def gen(smiles, name, nconfs):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=nconfs, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol


query = gen("C1=CC(=C(C=C1CC[NH3+])O)O", "dopamine_H+", 1)                 # amines protonated
dataset = [gen("CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O", "PD128907_H+", 20)]

calc = Roshambo2(query, dataset, color=True, remove_Hs_before_color_assignment=False)
scores = calc.compute(backend="cuda", reduce_over_conformers=True,
                      optim_mode="combination", combination_param=0.5, write_scores=False)

for qname, df in scores.items():
    print(f"\nquery = {qname}")
    print(df[["name", "tanimoto_shape", "tanimoto_color",
              "tanimoto_combination", "overlap_volume", "overlap_color"]].to_string(index=False))
