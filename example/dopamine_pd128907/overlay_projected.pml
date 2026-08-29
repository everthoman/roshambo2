# roshambo2 overlay with projected pharmacophore points
load ligand_dopamine_proj.sdf, dopamine
load ligand_PD128907_proj.sdf, PD128907
hide everything, dopamine or PD128907
show sticks, dopamine or PD128907
util.cbac dopamine
util.cbag PD128907

load projected_hits_dopamine_dopamine_H+_0_color_features.sdf, feats
set all_states, on
hide everything, feats
show spheres, feats
set sphere_scale, 0.35, feats
set sphere_transparency, 0.25, feats
color slate,     feats and elem N
color firebrick, feats and elem O
color orange,    feats and elem Fe
color green,     feats and elem Cl
color yellow,    feats and elem S
color grey60,    feats and elem Br
color deepblue,  feats and elem F
color red,       feats and elem I
color wheat,     feats and elem B
set sphere_scale, 0.5, feats and elem Fe
bg_color white
set orthoscopic, 1
zoom all, 3
