"""dopamine vs PD128907: overlap-only pharmacophore models.

Overlays protonated dopamine (query, N_QUERY_CONFS conformers, best kept) on the
D3 agonist PD128907 (dataset, N_DATASET_CONFS conformers) with the CUDA backend,
then keeps ONLY the "color" feature points that actually overlap between the two
(overlaptools.match_overlap) and renders those.

Two models are produced:
  stock      - roshambo2's default 1-point-per-feature colour model
  projected  - ProjectedPointPharmacophoreGenerator (adds directional donor /
               acceptor lone-pair / ring-normal / cation projected points)

Outputs (cwd):
  overlap_model_stock.sdf      overlap_model_projected.sdf        the models
  overlap_stock.pml  / .png    overlap_projected.pml  / .png      renders
  ligand_{dopamine,PD128907}_{stock,proj}.sdf                     context ligands
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator
from ligtools import pose_with_H, polar_only
import overlaptools as ot

MIN_OVERLAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15   # d <= ~1.95 A
N_QUERY_CONFS = 25
N_DATASET_CONFS = 100
DOPAMINE = "C1=CC(=C(C=C1CC[NH3+])O)O"
PD128907 = "CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O"


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


dopamine_confs = split(gen(DOPAMINE, "dopamine_H+", N_QUERY_CONFS))
dataset = [gen(PD128907, "PD128907_H+", N_DATASET_CONFS)]


def best_conformer(color_generator):
    """25-conformer scan -> (best dopamine conf name, best score row)."""
    calc = Roshambo2(list(dopamine_confs.values()), dataset, color=True,
                     remove_Hs_before_color_assignment=False,
                     conformers_have_unique_names=True,
                     color_generator=color_generator)
    scores = calc.compute(backend="cuda", reduce_over_conformers=True,
                          optim_mode="combination", combination_param=0.5,
                          write_scores=False)
    return max(((q.rsplit("_", 1)[0], df.iloc[0]) for q, df in scores.items()),
               key=lambda x: x[1]["tanimoto_combination"])


def build(tag, color_generator):
    best_q, row = best_conformer(color_generator)
    print(f"{tag:<11}{row['tanimoto_shape']:>8.3f}{row['tanimoto_color']:>8.3f}"
          f"{row['tanimoto_combination']:>8.3f}   (dopamine {best_q})")

    qmol = Chem.Mol(dopamine_confs[best_q])       # copy - don't rename the shared conf
    qmol.SetProp("_Name", "dopamine_H+")
    calc = Roshambo2(qmol, dataset, color=True, remove_Hs_before_color_assignment=False,
                     color_generator=color_generator)
    df = next(iter(calc.compute(backend="cuda", reduce_over_conformers=False,
                                optim_mode="combination", combination_param=0.5,
                                write_scores=False).values()))
    calc.write_best_fit_structures(hits_sdf_prefix=f"_hits_{tag}", top_n=1,
                                   write_color_pseudomols=True, append_query=True,
                                   feature_to_symbol_map=ot.FEATURE_TO_SYMBOL)

    pairs = ot.match_overlap(f"_hits_{tag}_dopamine_H+_0_color_features.sdf", MIN_OVERLAP)
    print(f"           model: {tag}")
    ot.print_table(pairs, float(df["overlap_color"].iloc[0]))
    ot.write_model(pairs, f"overlap_model_{tag}.sdf", min_overlap=MIN_OVERLAP)

    # context ligands with the exact hydrogens roshambo2 scored
    hit_noH = calc.get_best_fit_structures(top_n=1)["dopamine_H+_0"][0]
    k = int(hit_noH.GetProp("name").rsplit("_", 1)[-1])
    lq, lh = f"ligand_dopamine_{tag}.sdf", f"ligand_PD128907_{tag}.sdf"
    with Chem.SDWriter(lq) as w:
        w.write(polar_only(qmol))
    with Chem.SDWriter(lh) as w:
        w.write(polar_only(pose_with_H(hit_noH, Chem.Mol(dataset[0], confId=k))))

    ot.write_pml(f"overlap_{tag}.pml", lq, lh, f"overlap_model_{tag}.sdf", pairs,
                 title=f"dopamine x PD128907 overlapping features - {tag} model")
    print(f"           wrote overlap_model_{tag}.sdf and overlap_{tag}.pml\n")


print(f"dopamine query confs = {N_QUERY_CONFS}, PD128907 dataset confs = {N_DATASET_CONFS}, "
      f"min overlap fraction = {MIN_OVERLAP}\n")
print(f"{'model':<11}{'shape':>8}{'color':>8}{'combo':>8}")
build("stock", None)
build("projected", ProjectedPointPharmacophoreGenerator())
