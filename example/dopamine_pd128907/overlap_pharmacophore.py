"""dopamine vs PD128907: overlap-only pharmacophore models.

Overlays protonated dopamine (query, N_QUERY_CONFS conformers, best kept) on the
D3 agonist PD128907 (dataset, N_DATASET_CONFS conformers) with the CUDA backend,
keeps ONLY the "color" feature points that overlap between the two, renders them.

    python overlap_pharmacophore.py [min_overlap] [tversky_alpha]

min_overlap   overlap-fraction cutoff for a shared feature (default 0.15)
tversky_alpha if given, rank/select by reference-Tversky ("fit") instead of the
              symmetric combo Tanimoto - alpha->1 rewards covering the query and
              ignores the hit's extra bulk (partial / sub-region matching)

Two colour models are built: `stock` (roshambo2 default, 1 pt/feature) and
`projected` (ProjectedPointPharmacophoreGenerator - directional points).
Outputs: overlap_model_{stock,projected}.sdf, overlap_{stock,projected}.{pml,png},
ligand_{dopamine,PD128907}_{stock,projected}.sdf
"""
import sys
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator
from ligtools import pose_with_H, polar_only
import overlaptools as ot

MIN_OVERLAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
TVERSKY = float(sys.argv[2]) if len(sys.argv) > 2 else None
RANK = "tversky_combo" if TVERSKY else "tanimoto_combination"
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


def rank(df):
    return ot.add_tversky(df, alpha=TVERSKY) if TVERSKY else df


def best_conformer(color_generator):
    """25-conformer scan -> best dopamine conf name by RANK."""
    calc = Roshambo2(list(dopamine_confs.values()), dataset, color=True,
                     remove_Hs_before_color_assignment=False,
                     conformers_have_unique_names=True,
                     color_generator=color_generator)
    scores = calc.compute(backend="cuda", reduce_over_conformers=True,
                          optim_mode="combination", combination_param=0.5,
                          write_scores=False)
    rows = [(q.rsplit("_", 1)[0], rank(df).iloc[0]) for q, df in scores.items()]
    return max(rows, key=lambda x: x[1][RANK])


def build(tag, color_generator):
    best_q, row = best_conformer(color_generator)
    extra = f" tversky={row['tversky_combo']:.3f}" if TVERSKY else ""
    print(f"{tag:<11}{row['tanimoto_shape']:>8.3f}{row['tanimoto_color']:>8.3f}"
          f"{row['tanimoto_combination']:>8.3f}{extra}   (dopamine {best_q})")

    qmol = Chem.Mol(dopamine_confs[best_q])
    qmol.SetProp("_Name", "dopamine_H+")
    calc = Roshambo2(qmol, dataset, color=True, remove_Hs_before_color_assignment=False,
                     color_generator=color_generator)
    df = rank(next(iter(calc.compute(backend="cuda", reduce_over_conformers=False,
                                     optim_mode="combination", combination_param=0.5,
                                     write_scores=False).values())))
    df = df.sort_values(RANK, ascending=False).reset_index(drop=True)
    best_hit = df["name"].iloc[0]

    poses = {m.GetProp("name"): m for m in calc.get_best_fit_structures()["dopamine_H+_0"]}
    k = int(best_hit.rsplit("_", 1)[-1])
    hit = pose_with_H(poses[best_hit], Chem.Mol(dataset[0], confId=k))

    qpts = ot.feature_points(qmol, color_generator)
    hpts = ot.feature_points(hit, color_generator)
    pairs = ot.match_overlap(qpts, hpts, MIN_OVERLAP)
    print(f"           model: {tag}")
    ot.print_table(pairs, float(df["overlap_color"].iloc[0]))
    ot.write_model(pairs, f"overlap_model_{tag}.sdf", min_overlap=MIN_OVERLAP)

    lq, lh = f"ligand_dopamine_{tag}.sdf", f"ligand_PD128907_{tag}.sdf"
    with Chem.SDWriter(lq) as w:
        w.write(polar_only(qmol))
    with Chem.SDWriter(lh) as w:
        w.write(polar_only(hit))

    ot.write_pml(f"overlap_{tag}.pml", lq, lh, f"overlap_model_{tag}.sdf", pairs,
                 title=f"dopamine x PD128907 overlapping features - {tag} model")
    print(f"           wrote overlap_model_{tag}.sdf and overlap_{tag}.pml\n")


mode = f"Tversky alpha={TVERSKY}" if TVERSKY else "combo Tanimoto"
print(f"dopamine confs {N_QUERY_CONFS}, PD128907 confs {N_DATASET_CONFS}, "
      f"min overlap {MIN_OVERLAP}, ranked by {mode}\n")
print(f"{'model':<11}{'shape':>8}{'color':>8}{'combo':>8}")
build("stock", None)
build("projected", ProjectedPointPharmacophoreGenerator())
