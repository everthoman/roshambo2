# roshambo2 overlap pharmacophore (only features shared by dopamine & PD128907)
load ligand_dopamine.sdf, dopamine
load ligand_PD128907.sdf, PD128907
hide everything, dopamine or PD128907
show sticks, dopamine or PD128907
util.cbac dopamine
util.cbag PD128907

set stick_radius, 0.1, dopamine or PD128907

load pharmacophore_overlap_model.sdf, pharm
unbond pharm, pharm
show spheres, pharm
set sphere_scale, 0.4, pharm
set sphere_transparency, 0.2, pharm
color slate,     pharm and elem N
color firebrick, pharm and elem O
color orange,    pharm and elem Fe
color green,     pharm and elem Cl
color yellow,    pharm and elem S
color grey60,    pharm and elem Br
python
seen=set()
fams={'N':'Donor','O':'Acceptor','S':'Aromatic','BR':'Hydrophobe'}
for at in cmd.get_model('pharm').atom:
    fam=fams.get(at.symbol.upper())
    if fam and fam not in seen:
        seen.add(fam); cmd.label(f'pharm and index {at.index}', repr(fam))
python end
set label_size, 16
set label_color, black
set label_position, (0,0,3)
bg_color white
set orthoscopic, 1
set ray_shadows, 0
orient
zoom all, 2.5
