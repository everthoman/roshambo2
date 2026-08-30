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


def add_tversky(df, alpha=0.95, mixing=0.5):
    """Add reference-Tversky ("fit"/partial-match) columns to a roshambo2 score
    frame. alpha weights the *query* self-overlap: alpha->1 = "how much of the
    query is covered", ignoring the hit's extra bulk (so a small query can match
    a sub-region of a big molecule without a Tanimoto penalty). alpha=0.5 with
    equal terms would reduce to Tanimoto.

        Tv = O_AB / (alpha * O_AA + (1 - alpha) * O_BB)
    """
    df = df.copy()
    df["tversky_shape"] = df["overlap_volume"] / (
        alpha * df["self_overlap_volume_query"]
        + (1 - alpha) * df["self_overlap_volume_fit"])
    if "overlap_color" in df.columns:
        qc = df["self_overlap_color_query"].where(df["self_overlap_color_query"] > 0)
        df["tversky_color"] = (df["overlap_color"] /
                               (alpha * qc + (1 - alpha) * df["self_overlap_color_fit"])
                               ).fillna(0.0)
        df["tversky_combo"] = (1 - mixing) * df["tversky_shape"] + mixing * df["tversky_color"]
    return df

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


def feature_points(mol, color_generator=None):
    """[(family, xyz), ...] for an RDKit mol (needs its polar H for the donor /
    cation projected points). color_generator=None -> roshambo2's default."""
    if color_generator is None:
        from roshambo2.pharmacophore import PharmacophoreGenerator
        color_generator = PharmacophoreGenerator()
    coords, types = color_generator.generate_color_atoms(mol)
    i2f = color_generator.get_index_to_feature()
    return [(i2f[int(t)], np.asarray(xyz, float)) for xyz, t in zip(coords, types)]


def match_overlap(qpts, hpts, min_overlap=0.15):
    """qpts, hpts: [(family, xyz), ...] for query and aligned hit (same frame).
    Greedy nearest 1:1 match within each family; keep pairs with overlap
    fraction exp(-d²/2) >= min_overlap. Returns pair dicts, biggest first."""
    pairs = []
    fams = {f for f, _ in qpts} & {f for f, _ in hpts}
    for fam in fams:
        qi = [i for i, (f, _) in enumerate(qpts) if f == fam]
        hi = [j for j, (f, _) in enumerate(hpts) if f == fam]
        cand = sorted(((float(np.linalg.norm(qpts[i][1] - hpts[j][1])), i, j)
                       for i in qi for j in hi), key=lambda x: x[0])
        uq, uh = set(), set()
        for d, i, j in cand:
            if i in uq or j in uh:
                continue
            uq.add(i); uh.add(j)
            fr = float(np.exp(-d * d / 2.0))
            if fr >= min_overlap:
                pairs.append(dict(feature=fam, symbol=FEATURE_TO_SYMBOL[fam], d=d,
                                  frac=fr, v=V0 * fr,
                                  xyz=0.5 * (qpts[i][1] + hpts[j][1])))
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
