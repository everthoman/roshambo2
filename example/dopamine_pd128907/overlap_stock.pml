# dopamine x PD128907 overlapping features - stock model
load ligand_dopamine_stock.sdf, query
load ligand_PD128907_stock.sdf, hit
hide everything, query or hit
show sticks, query or hit
set stick_radius, 0.10, query or hit
set valence, 0
util.cbac query
util.cbag hit
color grey30, (query or hit) and elem H
show spheres, (query or hit) and elem H
set sphere_scale, 0.11, (query or hit) and elem H

load overlap_model_stock.sdf, pharm
unbond pharm, pharm
show spheres, pharm
set sphere_scale, 0.40, pharm
set sphere_transparency, 0.15, pharm
color yellow, pharm and elem S
color orange, pharm and elem Fe
color slate, pharm and elem N
color grey60, pharm and elem Br
color firebrick, pharm and elem O
python
lab = {'N': 'Donor', 'O': 'Acceptor', 'FE': 'Pos', 'CL': 'Neg', 'S': 'Aromatic', 'BR': 'Hydrophobe', 'F': 'D-proj', 'I': 'A-proj', 'B': 'ring-proj', 'P': 'Pos-proj'}
off = {'Donor': (0, 1.8, 3), 'Acceptor': (0, -2.4, 3), 'Aromatic': (0, 0, 3), 'Hydrophobe': (0, 2.4, 3), 'Pos': (2.6, -1.6, 3), 'Neg': (-2.6, -1.6, 3), 'D-proj': (-2.6, 0.4, 3), 'A-proj': (0, -2.0, 3), 'ring-proj': (2.4, 1.4, 3), 'Pos-proj': (2.8, -0.2, 3)}
seen = set()
for at in cmd.get_model('pharm').atom:
    name = lab.get(at.symbol.upper())
    if name and name not in seen:
        seen.add(name)
        sel = f'pharm and index {at.index}'
        cmd.label(sel, repr(name))
        cmd.set('label_position', off.get(name, (0, 0, 3)), sel)
python end
set label_size, 15
set label_color, black
bg_color white
set orthoscopic, 1
set ray_shadows, 0
orient
zoom all, 2.5
