"""PD128907 (rigid template) vs dopamine: overlap-only pharmacophore models.

The rigid D3 agonist PD128907 is the query/template (its lowest-MMFF-energy
conformer); flexible dopamine is the hit - roshambo2 samples N_HIT_CONFS dopamine
conformers and finds the best-fitting pose. Only the "color" feature points that
overlap between the two are kept and rendered.

    python overlap_pharmacophore.py [min_overlap] [tversky_alpha]

min_overlap    overlap-fraction cutoff for a shared feature (default 0.15)
tversky_alpha  if given, rank/select by reference-Tversky ("fit") instead of the
               symmetric combo Tanimoto:
                 Tv = O_AB / (alpha*O_query + (1-alpha)*O_hit)
               alpha -> 1 : how well is the PD128907 template covered
               alpha -> 0 : how well does dopamine fit *inside* the template
                            (the "is the small molecule a sub-shape" question)

Two colour models: `stock` (roshambo2 default, 1 pt/feature) and `projected`
(ProjectedPointPharmacophoreGenerator - directional points).
Outputs: overlap_model_{stock,projected}.sdf, overlap_{stock,projected}.{pml,png},
ligand_{PD128907,dopamine}_{stock,projected}.sdf
"""
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator
from ligtools import pose_with_H, polar_only
import overlaptools as ot

MIN_OVERLAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
TVERSKY = float(sys.argv[2]) if len(sys.argv) > 2 else None
RANK = "tversky_combo" if TVERSKY else "tanimoto_combination"
N_TEMPLATE_CONFS = 10       # PD128907 is rigid - just pick the lowest-energy one
N_HIT_CONFS = 200           # dopamine is flexible
PD128907 = "CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O"
DOPAMINE = "C1=CC(=C(C=C1CC[NH3+])O)O"


def gen(smiles, name, n):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=n, params=p)
    ff = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol, [e for _, e in ff]


_pd, _pd_e = gen(PD128907, "PD128907_H+", N_TEMPLATE_CONFS)
template = Chem.Mol(_pd, confId=int(np.argmin(_pd_e)))     # rigid query
template.SetProp("_Name", "PD128907_H+")
dopamine = gen(DOPAMINE, "dopamine_H+", N_HIT_CONFS)[0]    # flexible hit
QKEY = "PD128907_H+_0"


def build(tag, color_generator):
    calc = Roshambo2(template, [dopamine], color=True,
                     remove_Hs_before_color_assignment=False,
                     color_generator=color_generator)
    df = next(iter(calc.compute(backend="cuda", reduce_over_conformers=False,
                                optim_mode="combination", combination_param=0.5,
                                write_scores=False).values()))
    if TVERSKY:
        df = ot.add_tversky(df, alpha=TVERSKY)
    df = df.sort_values(RANK, ascending=False).reset_index(drop=True)
    best = df.iloc[0]
    extra = f" tversky={best['tversky_combo']:.3f}" if TVERSKY else ""
    print(f"{tag:<11}{best['tanimoto_shape']:>8.3f}{best['tanimoto_color']:>8.3f}"
          f"{best['tanimoto_combination']:>8.3f}{extra}   (dopamine conf {best['name'].rsplit('_', 1)[-1]})")

    poses = {m.GetProp("name"): m for m in calc.get_best_fit_structures()[QKEY]}
    k = int(best["name"].rsplit("_", 1)[-1])
    hit = pose_with_H(poses[best["name"]], Chem.Mol(dopamine, confId=k))

    qpts = ot.feature_points(template, color_generator)     # PD128907 (query/template)
    hpts = ot.feature_points(hit, color_generator)          # dopamine (aligned hit)
    pairs = ot.match_overlap(qpts, hpts, MIN_OVERLAP)
    print(f"           model: {tag}")
    ot.print_table(pairs, float(best["overlap_color"]))
    ot.write_model(pairs, f"overlap_model_{tag}.sdf", min_overlap=MIN_OVERLAP)

    lq, lh = f"ligand_PD128907_{tag}.sdf", f"ligand_dopamine_{tag}.sdf"
    with Chem.SDWriter(lq) as w:
        w.write(polar_only(template))
    with Chem.SDWriter(lh) as w:
        w.write(polar_only(hit))
    ot.write_pml(f"overlap_{tag}.pml", lq, lh, f"overlap_model_{tag}.sdf", pairs,
                 title=f"PD128907 (template) x dopamine overlapping features - {tag} model")
    print(f"           wrote overlap_model_{tag}.sdf and overlap_{tag}.pml\n")


mode = f"Tversky alpha={TVERSKY}" if TVERSKY else "combo Tanimoto"
print(f"template = PD128907 (rigid, lowest-energy conf);  hit = dopamine "
      f"({N_HIT_CONFS} confs);  min overlap {MIN_OVERLAP};  ranked by {mode}\n")
print(f"{'model':<11}{'shape':>8}{'color':>8}{'combo':>8}")
build("stock", None)
build("projected", ProjectedPointPharmacophoreGenerator())
