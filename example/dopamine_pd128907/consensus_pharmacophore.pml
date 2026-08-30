# dopaminergic consensus pharmacophore (roshambo2 pairwise alignment to apomorphine, then feature clustering)
load aligned_apomorphine.sdf, apomorphine
hide everything, apomorphine
show sticks, apomorphine
set stick_radius, 0.07, apomorphine
color cyan, apomorphine and elem C
load aligned_dopamine.sdf, dopamine
hide everything, dopamine
show sticks, dopamine
set stick_radius, 0.07, dopamine
color salmon, dopamine and elem C
load aligned_PD128907.sdf, PD128907
hide everything, PD128907
show sticks, PD128907
set stick_radius, 0.07, PD128907
color palegreen, PD128907 and elem C
load aligned_quinpirole.sdf, quinpirole
hide everything, quinpirole
show sticks, quinpirole
set stick_radius, 0.07, quinpirole
color lightblue, quinpirole and elem C
load aligned_rotigotine.sdf, rotigotine
hide everything, rotigotine
show sticks, rotigotine
set stick_radius, 0.07, rotigotine
color wheat, rotigotine and elem C
load aligned_R-7-OH-DPAT.sdf, R_7_OH_DPAT
hide everything, R_7_OH_DPAT
show sticks, R_7_OH_DPAT
set stick_radius, 0.07, R_7_OH_DPAT
color violet, R_7_OH_DPAT and elem C
set valence, 0
show spheres, (elem H and neighbor (elem N+O))
set sphere_scale, 0.10, (elem H)
color grey40, elem H

load consensus_model.sdf, pharm
unbond pharm, pharm
show spheres, pharm
set sphere_scale, 0.45, pharm
set sphere_transparency, 0.1, pharm
color yellow, pharm and elem S
color wheat, pharm and elem B
color grey60, pharm and elem Br
color firebrick, pharm and elem O
color slate, pharm and elem N
color hotpink, pharm and elem I
color deepblue, pharm and elem F
python
lab = {'N': 'Donor', 'O': 'Acceptor', 'FE': 'Pos', 'CL': 'Neg', 'S': 'Aromatic', 'BR': 'Hydrophobe', 'F': 'D-proj', 'I': 'A-proj', 'B': 'ring-proj', 'P': 'Pos-proj'}
seen = set()
for at in cmd.get_model('pharm').atom:
    nm = lab.get(at.symbol.upper())
    if nm and nm not in seen:
        seen.add(nm)
        cmd.label(f'pharm and index {at.index}', repr(nm))
        cmd.set('label_position', (0, 0, 3), f'pharm and index {at.index}')
python end
set label_size, 15
set label_color, black
bg_color white
set orthoscopic, 1
set ray_shadows, 0
orient
zoom all, 2.5
