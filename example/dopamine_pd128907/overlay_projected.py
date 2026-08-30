"""dopamine vs PD128907 overlay with ROCS-style projected points.

Compares the stock single-point colour model against
ProjectedPointPharmacophoreGenerator (base points + projected donor / acceptor /
ring-normal points), then writes the projected feature cloud for PyMOL.

dopamine (flexible) is sampled with N_QUERY_CONFS conformers, each a separate
named query; the best-scoring one is reported and used for the feature cloud.
PD128907 (rigid) is the dataset with N_DATASET_CONFS conformers, reduced over.

Outputs (cwd, all prefixed so the stock example's files are untouched):
  ligand_dopamine_proj.sdf / ligand_PD128907_proj.sdf   the two ligands
  projected_hits_dopamine_dopamine_H+_0.sdf             query + aligned hit
  projected_hits_dopamine_dopamine_H+_0_color_features.sdf   all feature points
  overlay_projected.pml                                 PyMOL loader
"""
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from roshambo2 import Roshambo2
from projected_pharmacophore import ProjectedPointPharmacophoreGenerator

N_QUERY_CONFS = 25       # dopamine (flexible) - sampled, best pose reported
N_DATASET_CONFS = 100    # PD128907 (rigid)   - reduced over

SYMBOL = {                       # feature family -> element (for PyMOL colouring)
    "Donor": "N", "Acceptor": "O", "PosIonizable": "Fe",
    "NegIonizable": "Cl", "Aromatic": "S", "Hydrophobe": "Br",
    "DonorProj": "F", "AcceptorProj": "I", "AromaticProj": "B",
    "PosIonizableProj": "P",
}


def gen(smiles, name, nconfs):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 0xF00D
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=nconfs, params=p)
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    mol.SetProp("_Name", name)
    return mol


def split(mol):
    """multi-conf mol -> {name_i: single-conf mol}"""
    out = {}
    for c in mol.GetConformers():
        mi = Chem.Mol(mol, confId=c.GetId())
        nm = f"{mol.GetProp('_Name')}_{c.GetId()}"
        mi.SetProp("_Name", nm)
        out[nm] = mi
    return out


dopamine_confs = split(gen("C1=CC(=C(C=C1CC[NH3+])O)O", "dopamine_H+", N_QUERY_CONFS))
queries = list(dopamine_confs.values())
dataset = [gen("CCC[NH+]1CCO[C@@H]2[C@@H]1COC3=C2C=C(C=C3)O", "PD128907_H+", N_DATASET_CONFS)]


def run(color_generator, label):
    calc = Roshambo2(queries, dataset, color=True,
                     remove_Hs_before_color_assignment=False,
                     conformers_have_unique_names=True,
                     color_generator=color_generator)
    scores = calc.compute(backend="cuda", reduce_over_conformers=True,
                          optim_mode="combination", combination_param=0.5,
                          write_scores=False)
    # best dopamine conformer (roshambo2 appends a trailing _<confId>; strip it)
    qname, row = max(((q.rsplit("_", 1)[0], df.iloc[0]) for q, df in scores.items()),
                     key=lambda x: x[1]["tanimoto_combination"])
    print(f"{label:<22}{row['tanimoto_shape']:>10.3f}{row['tanimoto_color']:>10.3f}"
          f"{row['tanimoto_combination']:>12.3f}{row['overlap_color']:>12.2f}"
          f"   best={qname}")
    return qname


print(f"dopamine query confs = {N_QUERY_CONFS}, PD128907 dataset confs = {N_DATASET_CONFS}\n")
print(f"{'model':<22}{'shape':>10}{'color':>10}{'combo':>12}{'ovlp_color':>12}")
run(None, "stock (1 pt/feature)")
best_q = run(ProjectedPointPharmacophoreGenerator(), "base + projected pts")
run(ProjectedPointPharmacophoreGenerator(keep_base_points=False), "projected pts only")

# ---- feature cloud from the best dopamine conformer (base+projected model) ----
best_mol = dopamine_confs[best_q]
best_mol.SetProp("_Name", "dopamine_H+")          # clean output filenames
calc = Roshambo2(best_mol, dataset, color=True,
                 remove_Hs_before_color_assignment=False,
                 color_generator=ProjectedPointPharmacophoreGenerator())
calc.compute(backend="cuda", reduce_over_conformers=True,
             optim_mode="combination", combination_param=0.5, write_scores=False)

PFX = "projected_hits_dopamine"                   # dedicated prefix, never clobbers stock files
calc.write_best_fit_structures(hits_sdf_prefix=PFX,
                               write_color_pseudomols=True, append_query=True,
                               feature_to_symbol_map=SYMBOL)
hits_sdf = f"{PFX}_dopamine_H+_0.sdf"
feats_sdf = f"{PFX}_dopamine_H+_0_color_features.sdf"

_hits = list(Chem.SDMolSupplier(hits_sdf, removeHs=False, sanitize=True))
for _m, _fn in zip(_hits, ("ligand_dopamine_proj.sdf", "ligand_PD128907_proj.sdf")):
    polar = [a.GetIdx() for a in _m.GetAtoms() if a.GetAtomicNum() in (7, 8)]
    with Chem.SDWriter(_fn) as _w:
        _w.write(Chem.AddHs(_m, addCoords=True, onlyOnAtoms=polar))   # polar H only

with open("overlay_projected.pml", "w") as fh:
    fh.write(
        "# roshambo2 overlay with projected pharmacophore points\n"
        "load ligand_dopamine_proj.sdf, dopamine\n"
        "load ligand_PD128907_proj.sdf, PD128907\n"
        "hide everything, dopamine or PD128907\n"
        "show sticks, dopamine or PD128907\n"
        "set stick_radius, 0.12, dopamine or PD128907\n"
        "set valence, 0\n"
        "util.cbac dopamine\nutil.cbag PD128907\n"
        "color grey30, (dopamine or PD128907) and elem H\n"
        "show spheres, (dopamine or PD128907) and elem H\n"
        "set sphere_scale, 0.11, (dopamine or PD128907) and elem H\n\n"
        f"load {feats_sdf}, feats\n"
        "set all_states, on\n"
        "unbond feats, feats\n"                     # kill spurious pseudo-atom bonds
        "hide everything, feats\nshow spheres, feats\n"
        "set sphere_scale, 0.3, feats\nset sphere_transparency, 0.4, feats\n"
        "color slate,     feats and elem N\n"       # Donor
        "color firebrick, feats and elem O\n"       # Acceptor
        "color orange,    feats and elem Fe\n"      # PosIonizable
        "color green,     feats and elem Cl\n"      # NegIonizable
        "color yellow,    feats and elem S\n"       # Aromatic
        "color grey70,    feats and elem Br\n"      # Hydrophobe
        "color deepblue,  feats and elem F\n"       # DonorProj
        "color hotpink,   feats and elem I\n"       # AcceptorProj
        "color wheat,     feats and elem B\n"       # AromaticProj
        "color purple,    feats and elem P\n"       # PosIonizableProj
        "set sphere_scale, 0.45, feats and elem Fe\n"
        "set sphere_transparency, 0.1, feats and (elem Fe or elem F or elem I or elem B or elem P)\n"
        "bg_color white\nset orthoscopic, 1\nset ray_shadows, 0\norient\nzoom all, 2\n"
    )
print(f"\nwrote {hits_sdf}, {feats_sdf} and overlay_projected.pml")
