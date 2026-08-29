# roshambo2 overlap pharmacophore (only features shared by dopamine & PD128907)
load ligand_dopamine.sdf, dopamine
load ligand_PD128907.sdf, PD128907
hide everything, dopamine or PD128907
show sticks, dopamine or PD128907
util.cbac dopamine
util.cbag PD128907

load pharmacophore_overlap_model.sdf, pharm
show spheres, pharm
set sphere_scale, 0.5, pharm
set sphere_transparency, 0.3, pharm
color slate,     pharm and elem N
color firebrick, pharm and elem O
color orange,    pharm and elem Fe
color green,     pharm and elem Cl
color yellow,    pharm and elem S
color grey60,    pharm and elem Br
label pharm, {'N':'Donor','O':'Acceptor','FE':'Pos','CL':'Neg','S':'Aromatic','BR':'Hydrophobe'}.get(elem.upper(), elem)
set label_size, 18
bg_color white
set orthoscopic, 1
zoom all, 3
