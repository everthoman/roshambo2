# PyMOL: roshambo2 overlay + full pharmacophore ("color") feature clouds
# run from this directory:   pymol view.pml

# the two ligands as separate objects
load ligand_dopamine.sdf, dopamine
load ligand_PD128907.sdf, PD128907
hide everything, dopamine or PD128907
show sticks, dopamine or PD128907
util.cbac dopamine
util.cbag PD128907

# all per-molecule feature points (entry 1 = dopamine, entry 2 = PD128907)
load hits_for_query_dopamine_H+_0_color_features.sdf, feats
set all_states, on
hide everything, feats
show spheres, feats
set sphere_scale, 0.4, feats
set sphere_transparency, 0.25, feats

# element -> feature family colour key
color slate,     feats and elem N      # Donor
color firebrick, feats and elem O      # Acceptor
color orange,    feats and elem Fe     # PosIonizable  (protonated amine)
color green,     feats and elem Cl     # NegIonizable
color yellow,    feats and elem S      # Aromatic
color grey60,    feats and elem Br     # Hydrophobe

set sphere_scale, 0.5, feats and elem Fe
bg_color white
set orthoscopic, 1
zoom all, 2
