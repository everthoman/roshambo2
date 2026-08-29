"""dopamine vs PD128907 overlay with ROCS-style projected points.

Compares the stock single-point colour model against
ProjectedPointPharmacophoreGenerator (base points + projected donor / acceptor /
ring-normal points), then writes the projected feature cloud for PyMOL.

Outputs (cwd, all prefixed so the stock example's files are untouched):
  ligand_dopamine_proj.sdf / ligand_PD128907_proj.sdf   the two ligands
  projected_hits_dopamine_dopamine_H+_0.sdf             query + aligned hit
  projected_hits_dopamine_dopamine_H+_0_color_features.sdf   all feature points
  overlay_projected.pml                                 PyMOL loader
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator

SYMBOL = {                       # feature family -> element (for PyMOL colouring)
    "Donor": "N", "Acceptor": "O", "PosIonizable": "Fe",
    "NegIonizable": "Cl", "Aromatic": "S", "Hydrophobe": "Br",
    "DonorProj": "F", "AcceptorProj": "I", "AromaticProj": "B",
}


def gen(smiles, name, nconfs):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=nconfs, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol


query = gen("C1=CC(=C(C=C1CC[NH3+])O)O", "dopamine_H+", 1)
dataset = [gen("CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O", "PD128907_H+", 20)]


def run(color_generator, label):
    calc = Roshambo2(query, dataset, color=True,
                     remove_Hs_before_color_assignment=False,
                     color_generator=color_generator)
    df = next(iter(calc.compute(backend="cuda", reduce_over_conformers=True,
                                optim_mode="combination", combination_param=0.5,
                                write_scores=False).values()))
    row = df.iloc[0]
    print(f"{label:<22}{row['tanimoto_shape']:>10.3f}{row['tanimoto_color']:>10.3f}"
          f"{row['tanimoto_combination']:>12.3f}{row['overlap_color']:>12.2f}")
    return calc


print(f"{'model':<22}{'shape':>10}{'color':>10}{'combo':>12}{'ovlp_color':>12}")
run(None, "stock (1 pt/feature)")
calc = run(ProjectedPointPharmacophoreGenerator(), "base + projected pts")
run(ProjectedPointPharmacophoreGenerator(keep_base_points=False),
    "projected pts only")

# ---- write the projected feature cloud from the base+projected run ----------
# use a dedicated prefix so this example never touches the stock example's files
PFX = "projected_hits_dopamine"
calc.write_best_fit_structures(hits_sdf_prefix=PFX,
                               write_color_pseudomols=True, append_query=True,
                               feature_to_symbol_map=SYMBOL)
hits_sdf = f"{PFX}_dopamine_H+_0.sdf"
feats_sdf = f"{PFX}_dopamine_H+_0_color_features.sdf"

_hits = list(Chem.SDMolSupplier(hits_sdf, removeHs=False, sanitize=False))
for _m, _fn in zip(_hits, ("ligand_dopamine_proj.sdf", "ligand_PD128907_proj.sdf")):
    with Chem.SDWriter(_fn) as _w:
        _w.write(_m)

with open("overlay_projected.pml", "w") as fh:
    fh.write(
        "# roshambo2 overlay with projected pharmacophore points\n"
        "load ligand_dopamine_proj.sdf, dopamine\n"
        "load ligand_PD128907_proj.sdf, PD128907\n"
        "hide everything, dopamine or PD128907\n"
        "show sticks, dopamine or PD128907\n"
        "util.cbac dopamine\nutil.cbag PD128907\n\n"
        f"load {feats_sdf}, feats\n"
        "set all_states, on\n"
        "hide everything, feats\nshow spheres, feats\n"
        "set sphere_scale, 0.35, feats\nset sphere_transparency, 0.25, feats\n"
        "color slate,     feats and elem N\n"       # Donor
        "color firebrick, feats and elem O\n"       # Acceptor
        "color orange,    feats and elem Fe\n"      # PosIonizable
        "color green,     feats and elem Cl\n"      # NegIonizable
        "color yellow,    feats and elem S\n"       # Aromatic
        "color grey60,    feats and elem Br\n"      # Hydrophobe
        "color deepblue,  feats and elem F\n"       # DonorProj
        "color red,       feats and elem I\n"       # AcceptorProj
        "color wheat,     feats and elem B\n"       # AromaticProj
        "set sphere_scale, 0.5, feats and elem Fe\n"
        "bg_color white\nset orthoscopic, 1\nzoom all, 3\n"
    )
print(f"\nwrote {hits_sdf}, {feats_sdf} and overlay_projected.pml")
