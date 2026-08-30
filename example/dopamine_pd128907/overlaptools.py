"""Build a 3D pharmacophore model from ONLY the roshambo2 "color" feature points
that actually overlap between a query and its aligned hit.

roshambo2's colour term is a sum of Gaussian overlaps between same-family feature
points (cpp_helper_functions.cpp::volume_color). For the default interaction map
every same-family pair has width a=1, height p=1, point weights w=1, so one
query/hit feature pair separated by d (Angstrom) contributes

    v(d) = (pi/2)**1.5 * exp(-d**2 / 2)          # (pi/2)**1.5 = 1.9687 = v(0)

A pair is "overlapping" when its overlap fraction f = exp(-d**2/2) >= min_overlap.
Features are matched query<->hit greedily by nearest distance within each family
and written at the pair midpoint.
"""
import numpy as np
from rdkit import Chem

V0 = (np.pi / 2) ** 1.5

# feature family -> (element used in the SDF, PyMOL colour, short label)
FAMILY = {
    "Donor":            ("N",  "slate",     "Donor"),
    "Acceptor":         ("O",  "firebrick", "Acceptor"),
    "PosIonizable":     ("Fe", "orange",    "Pos"),
    "NegIonizable":     ("Cl", "green",     "Neg"),
    "Aromatic":         ("S",  "yellow",    "Aromatic"),
    "Hydrophobe":       ("Br", "grey60",    "Hydrophobe"),
    "DonorProj":        ("F",  "deepblue",  "D-proj"),
    "AcceptorProj":     ("I",  "hotpink",   "A-proj"),
    "AromaticProj":     ("B",  "wheat",     "ring-proj"),
    "PosIonizableProj": ("P",  "purple",    "Pos-proj"),
}
FEATURE_TO_SYMBOL = {f: e for f, (e, _, _) in FAMILY.items()}
SYMBOL_TO_FEATURE = {e: f for f, e in FEATURE_TO_SYMBOL.items()}
_Z = {"N": 7, "O": 8, "Fe": 26, "Cl": 17, "S": 16, "Br": 35,
      "F": 9, "I": 53, "B": 5, "P": 15}


def _points(mol):
    conf = mol.GetConformer()
    return [(a.GetSymbol(), np.array(conf.GetAtomPosition(a.GetIdx())))
            for a in mol.GetAtoms()]


def match_overlap(color_features_sdf, min_overlap=0.15):
    """color_features_sdf: the *_color_features.sdf written by roshambo2 with
    write_color_pseudomols=True + append_query=True (entry 0 = query, 1 = hit).
    Returns a list of pair dicts sorted by contribution, most first."""
    q, h = list(Chem.SDMolSupplier(color_features_sdf, removeHs=False, sanitize=False))[:2]
    qf, hf = _points(q), _points(h)
    pairs = []
    for sym in {s for s, _ in qf} & {s for s, _ in hf}:
        qi = [i for i, (s, _) in enumerate(qf) if s == sym]
        hi = [j for j, (s, _) in enumerate(hf) if s == sym]
        cand = sorted(((float(np.linalg.norm(qf[i][1] - hf[j][1])), i, j)
                       for i in qi for j in hi), key=lambda x: x[0])
        uq, uh = set(), set()
        for d, i, j in cand:
            if i in uq or j in uh:
                continue
            uq.add(i); uh.add(j)
            fr = float(np.exp(-d * d / 2.0))
            if fr >= min_overlap:
                pairs.append(dict(feature=SYMBOL_TO_FEATURE[sym], symbol=sym, d=d,
                                  frac=fr, v=V0 * fr, xyz=0.5 * (qf[i][1] + hf[j][1])))
    pairs.sort(key=lambda p: -p["v"])
    return pairs


def write_model(pairs, sdf_path, name="overlap_pharmacophore", min_overlap=None):
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(pairs))
    for p in pairs:
        conf.SetAtomPosition(rw.AddAtom(Chem.Atom(_Z[p["symbol"]])), p["xyz"].tolist())
    m = rw.GetMol()
    m.AddConformer(conf)
    m.SetProp("_Name", name)
    if min_overlap is not None:
        m.SetProp("min_overlap_fraction", str(min_overlap))
    m.SetProp("features", ", ".join(f"{p['feature']}(d={p['d']:.2f})" for p in pairs))
    with Chem.SDWriter(sdf_path) as w:
        w.write(m)


def print_table(pairs, overlap_color=None):
    print(f"\n{'feature':<16}{'dist(A)':>9}{'overlap_frac':>14}{'v_color':>10}")
    for p in pairs:
        print(f"{p['feature']:<16}{p['d']:>9.2f}{p['frac']:>14.3f}{p['v']:>10.3f}")
    tot = sum(p["v"] for p in pairs)
    tail = ""
    if overlap_color:
        tail = (f" / overlap_color {overlap_color:.2f} "
                f"({100 * tot / overlap_color:.0f}% of the colour overlap)")
    print(f"\n{len(pairs)} overlapping features   sum v_color = {tot:.2f}{tail}")


def write_pml(pml_path, ligand_q, ligand_h, model_sdf, pairs, title=""):
    present = []
    seen = set()
    for p in pairs:                       # families actually in the model, kept order
        if p["feature"] not in seen:
            seen.add(p["feature"]); present.append(p["feature"])
    lines = [
        f"# {title}" if title else "# roshambo2 overlapping-feature pharmacophore",
        f"load {ligand_q}, query",
        f"load {ligand_h}, hit",
        "hide everything, query or hit",
        "show sticks, query or hit",
        "set stick_radius, 0.10, query or hit",
        "set valence, 0",
        "util.cbac query",
        "util.cbag hit",
        "color grey30, (query or hit) and elem H",
        "show spheres, (query or hit) and elem H",
        "set sphere_scale, 0.11, (query or hit) and elem H",
        "",
        f"load {model_sdf}, pharm",
        "unbond pharm, pharm",
        "show spheres, pharm",
        "set sphere_scale, 0.40, pharm",
        "set sphere_transparency, 0.15, pharm",
    ]
    for fam in present:
        e, col, _ = FAMILY[fam]
        lines.append(f"color {col}, pharm and elem {e}")
    # one label per family, spread around so they do not pile up
    off = {"Donor": (0, 1.8, 3), "Acceptor": (0, -2.4, 3), "Aromatic": (0, 0, 3),
           "Hydrophobe": (0, 2.4, 3), "Pos": (2.6, -1.6, 3), "Neg": (-2.6, -1.6, 3),
           "D-proj": (-2.6, 0.4, 3), "A-proj": (0, -2.0, 3),
           "ring-proj": (2.4, 1.4, 3), "Pos-proj": (2.8, -0.2, 3)}
    lab = {e.upper(): FAMILY[f][2] for f, (e, _, _) in FAMILY.items()}
    lines += [
        "python",
        f"lab = {lab!r}",
        f"off = {off!r}",
        "seen = set()",
        "for at in cmd.get_model('pharm').atom:",
        "    name = lab.get(at.symbol.upper())",
        "    if name and name not in seen:",
        "        seen.add(name)",
        "        sel = f'pharm and index {at.index}'",
        "        cmd.label(sel, repr(name))",
        "        cmd.set('label_position', off.get(name, (0, 0, 3)), sel)",
        "python end",
        "set label_size, 15",
        "set label_color, black",
        "bg_color white",
        "set orthoscopic, 1",
        "set ray_shadows, 0",
        "orient",
        "zoom all, 2.5",
    ]
    open(pml_path, "w").write("\n".join(lines) + "\n")
