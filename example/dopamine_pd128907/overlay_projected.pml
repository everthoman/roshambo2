# roshambo2 overlay with projected pharmacophore points
load ligand_dopamine_proj.sdf, dopamine
load ligand_PD128907_proj.sdf, PD128907
hide everything, dopamine or PD128907
show sticks, dopamine or PD128907
set stick_radius, 0.12, dopamine or PD128907
set valence, 0
util.cbac dopamine
util.cbag PD128907
color grey30, (dopamine or PD128907) and elem H
show spheres, (dopamine or PD128907) and elem H
set sphere_scale, 0.11, (dopamine or PD128907) and elem H

load projected_hits_dopamine_dopamine_H+_0_color_features.sdf, feats
set all_states, on
unbond feats, feats
hide everything, feats
show spheres, feats
set sphere_scale, 0.3, feats
set sphere_transparency, 0.4, feats
color slate,     feats and elem N
color firebrick, feats and elem O
color orange,    feats and elem Fe
color green,     feats and elem Cl
color yellow,    feats and elem S
color grey70,    feats and elem Br
color deepblue,  feats and elem F
color hotpink,   feats and elem I
color wheat,     feats and elem B
color purple,    feats and elem P
set sphere_scale, 0.45, feats and elem Fe
set sphere_transparency, 0.1, feats and (elem Fe or elem F or elem I or elem B or elem P)
bg_color white
set orthoscopic, 1
set ray_shadows, 0
orient
zoom all, 2
