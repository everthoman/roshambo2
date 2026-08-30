"""PD128907 vs dopamine: overlap-only pharmacophore models, conformer + strain aware.

Both molecules contribute every MMFF conformer within E_WINDOW kcal/mol of their
own global minimum (RMS-pruned). roshambo2 aligns every dopamine conformer to
every PD128907 conformer; each alignment is then ranked by

    adjusted = score - ENERGY_LAMBDA * (dE_template + dE_hit)

i.e. the raw shape/colour score minus a linear penalty on the conformational
strain both partners have to pay. The best-adjusted alignment builds the model.

    python overlap_pharmacophore.py [min_overlap] [tversky_alpha] [energy_lambda]

min_overlap    overlap-fraction cutoff for a shared feature (default 0.15)
tversky_alpha  rank by reference-Tversky ("fit") instead of combo Tanimoto
               (alpha->0 = "does dopamine fit inside PD128907")
energy_lambda  score units per kcal/mol of total strain (default 0.05)

Two colour models: `stock` and `projected`. Outputs:
overlap_model_{stock,projected}.sdf, overlap_{stock,projected}.{pml,png},
ligand_{PD128907,dopamine}_{stock,projected}.sdf
"""
import sys
from rdkit import Chem
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator
from ligtools import pose_with_H, polar_only
import overlaptools as ot

MIN_OVERLAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15
TVERSKY = float(sys.argv[2]) if len(sys.argv) > 2 else None
ENERGY_LAMBDA = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
E_WINDOW = 5.0
RANK = "tversky_combo" if TVERSKY else "tanimoto_combination"
PD128907 = "CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O"
DOPAMINE = "C1=CC(=C(C=C1CC[NH3+])O)O"

tmpl_mol, tmpl_dE = ot.conformer_ensemble(PD128907, "PD128907_H+", E_WINDOW)
hit_mol, hit_dE = ot.conformer_ensemble(DOPAMINE, "dopamine_H+", E_WINDOW)
print(f"conformers within {E_WINDOW} kcal/mol:  PD128907 {len(tmpl_dE)}  "
      f"(dE {min(tmpl_dE):.1f}-{max(tmpl_dE):.1f}),  dopamine {len(hit_dE)} "
      f"(dE {min(hit_dE):.1f}-{max(hit_dE):.1f})")

# one roshambo2 query per template conformer
tmpl_confs = {}
for i in range(tmpl_mol.GetNumConformers()):
    mi = Chem.Mol(tmpl_mol, confId=tmpl_mol.GetConformer(i).GetId())
    mi.SetProp("_Name", f"tmpl{i}")
    tmpl_confs[f"tmpl{i}_0"] = (mi, tmpl_dE[i])


def build(tag, color_generator):
    calc = Roshambo2([m for m, _ in tmpl_confs.values()], [hit_mol], color=True,
                     remove_Hs_before_color_assignment=False,
                     conformers_have_unique_names=True,
                     color_generator=color_generator)
    scores = calc.compute(backend="cuda", reduce_over_conformers=False,
                          optim_mode="combination", combination_param=0.5,
                          write_scores=False)
    poses = calc.get_best_fit_structures()

    best = None
    for qkey, df in scores.items():
        dEq = tmpl_confs[qkey][1]
        if TVERSKY:
            df = ot.add_tversky(df, alpha=TVERSKY)
        for _, r in df.iterrows():
            hc = int(r["name"].rsplit("_", 1)[-1])
            adj = r[RANK] - ENERGY_LAMBDA * (dEq + hit_dE[hc])
            if best is None or adj > best["adj"]:
                best = dict(adj=adj, raw=r[RANK], dEq=dEq, dEh=hit_dE[hc],
                            qkey=qkey, hname=r["name"], hc=hc, row=r)

    r = best["row"]
    print(f"{tag:<11}{r['tanimoto_shape']:>7.3f}{r['tanimoto_color']:>7.3f}"
          f"{r['tanimoto_combination']:>7.3f}   strain dE {best['dEq']:.1f}+{best['dEh']:.1f}"
          f"   adjusted {best['adj']:.3f}")

    hit_noH = {m.GetProp("name"): m for m in poses[best["qkey"]]}[best["hname"]]
    hc_id = hit_mol.GetConformer(best["hc"]).GetId()
    hit = pose_with_H(hit_noH, Chem.Mol(hit_mol, confId=hc_id))
    template = tmpl_confs[best["qkey"]][0]

    qpts = ot.feature_points(template, color_generator)     # PD128907 (template)
    hpts = ot.feature_points(hit, color_generator)          # dopamine (aligned hit)
    pairs = ot.match_overlap(qpts, hpts, MIN_OVERLAP)
    print(f"           model: {tag}")
    ot.print_table(pairs, float(r["overlap_color"]))
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
print(f"min overlap {MIN_OVERLAP};  ranked by {mode} - {ENERGY_LAMBDA}/kcal * total strain\n")
print(f"{'model':<11}{'shape':>7}{'color':>7}{'combo':>7}")
build("stock", None)
build("projected", ProjectedPointPharmacophoreGenerator())
