"""How do dopamine/PD128907 overlay scores respond to conformer count?

Runs the base+projected colour model for a grid of
(dopamine query confs, PD128907 dataset confs). For >1 query conf each
conformer is a separately named query; we report the best-scoring one.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator

DOPAMINE = "C1=CC(=C(C=C1CC[NH3+])O)O"
PD128907 = "CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O"

QUERY_CONFS = [1, 5, 25]
DATA_CONFS = [20, 100, 400]


def confs(smiles, name, n):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(m, numConfs=n, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(m)
    m.SetProp("_Name", name)
    return m


def split(mol):
    out = []
    for c in mol.GetConformers():
        mi = Chem.Mol(mol, confId=c.GetId())
        mi.SetProp("_Name", f"{mol.GetProp('_Name')}_{c.GetId()}")
        out.append(mi)
    return out


print(f"{'q_confs':>8}{'d_confs':>9}{'shape':>9}{'color':>9}{'combo':>9}")
for nq in QUERY_CONFS:
    qmol = confs(DOPAMINE, "dopamine_H+", nq)
    queries = split(qmol) if nq > 1 else qmol
    for nd in DATA_CONFS:
        dset = [confs(PD128907, "PD128907_H+", nd)]
        calc = Roshambo2(queries, dset, color=True,
                         remove_Hs_before_color_assignment=False,
                         conformers_have_unique_names=(nq > 1),
                         color_generator=ProjectedPointPharmacophoreGenerator())
        scores = calc.compute(backend="cuda", reduce_over_conformers=True,
                              optim_mode="combination", combination_param=0.5,
                              write_scores=False)
        # best query conformer
        best = max((df.iloc[0] for df in scores.values()),
                   key=lambda r: r["tanimoto_combination"])
        print(f"{nq:>8}{nd:>9}{best['tanimoto_shape']:>9.3f}"
              f"{best['tanimoto_color']:>9.3f}{best['tanimoto_combination']:>9.3f}")
