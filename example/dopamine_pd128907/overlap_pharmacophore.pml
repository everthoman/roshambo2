# roshambo2 overlap pharmacophore - only the directional features
# shared by dopamine & PD128907 (base points + projected D/A/ring points)
load ligand_dopamine_ovl.sdf, dopamine
load ligand_PD128907_ovl.sdf, PD128907
hide everything, dopamine or PD128907
show sticks, dopamine or PD128907
set stick_radius, 0.1, dopamine or PD128907
set valence, 0
util.cbac dopamine
util.cbag PD128907
color grey30, (dopamine or PD128907) and elem H
show spheres, (dopamine or PD128907) and elem H
set sphere_scale, 0.11, (dopamine or PD128907) and elem H

load pharmacophore_overlap_model.sdf, pharm
unbond pharm, pharm
show spheres, pharm
set sphere_scale, 0.4, pharm
set sphere_transparency, 0.15, pharm
color slate,     pharm and elem N
color firebrick, pharm and elem O
color orange,    pharm and elem Fe
color green,     pharm and elem Cl
color yellow,    pharm and elem S
color grey60,    pharm and elem Br
color deepblue,  pharm and elem F
color hotpink,   pharm and elem I
color wheat,     pharm and elem B
color purple,    pharm and elem P
python
seen=set()
fams={'N':'Donor','O':'Acceptor','S':'Aromatic','BR':'Hydrophobe','F':'D-proj','I':'A-proj','B':'ring-proj','P':'Pos-proj'}
off={'Donor':(0,1.6,3),'Acceptor':(0,-2.4,3),'Aromatic':(0,0,3),'Hydrophobe':(0,2.2,3),'D-proj':(-2.4,0,3),'A-proj':(0,-2.0,3),'ring-proj':(2.2,1.2,3),'Pos-proj':(2.4,-1.4,3)}
for at in cmd.get_model('pharm').atom:
    fam=fams.get(at.symbol.upper())
    if fam and fam not in seen:
        seen.add(fam)
        sel=f'pharm and index {at.index}'
        cmd.label(sel, repr(fam)); cmd.set('label_position', off[fam], sel)
python end
set label_size, 15
set label_color, black
bg_color white
set orthoscopic, 1
set ray_shadows, 0
orient
zoom all, 2.5
