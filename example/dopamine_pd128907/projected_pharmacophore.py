"""ProjectedPointPharmacophoreGenerator for roshambo2.

roshambo2's stock colour model puts a single isotropic Gaussian on each
pharmacophore feature atom, so H-bond / ring-stacking *direction* is ignored.
This generator adds ROCS-style **projected points**:

  * Donor    -> one point `donor_dist` A out along each D-H bond
               (roughly where the acceptor would sit)
  * Acceptor -> points `acceptor_dist` A out along the lone-pair directions
               (see _acceptor_lone_pairs): two tetrahedral lone pairs for
               hydroxyl / ether O (out of the C-O-C plane, NOT coplanar with
               the substituents), two sp2 lone pairs at +/-120 deg for a
               carbonyl O / nitrile N
  * Aromatic -> two points `ring_dist` A above / below the ring plane
               along the ring normal
  * PosIonizable -> one point `pos_dist` A out along each N-H of a protonated
               amine (the salt-bridge vector), like a donor; straight out from
               the mean bond direction if the cation carries no H

The projected points get their own feature families (`DonorProj`,
`AcceptorProj`, `AromaticProj`) that only overlap their own kind, exactly like
the base families.  So the colour term now rewards two molecules for pointing
their donors / acceptors / ring faces the *same way*, not just for having the
heteroatom in the same place.

Usage
-----
    from projected_pharmacophore import ProjectedPointPharmacophoreGenerator
    from roshambo2 import Roshambo2

    cg = ProjectedPointPharmacophoreGenerator()          # base + projected
    calc = Roshambo2(query, dataset, color=True, color_generator=cg)
    calc.compute(backend="cuda", optim_mode="combination")

Knobs: donor_dist, acceptor_dist, ring_dist, pos_dist (Angstrom); proj_weight
(Gaussian height for the projected families, <1 down-weights them);
keep_base_points (set False for a projected-points-only model).  The molecule
handed to the generator keeps its explicit H atoms as long as Roshambo2 is
called with the default remove_Hs_before_color_assignment=False.
"""
import numpy as np
from roshambo2.pharmacophore import RDKitPharmacophoreGenerator

_BASE = ("Donor", "Acceptor", "PosIonizable", "NegIonizable", "Aromatic", "Hydrophobe")
_PROJ = ("DonorProj", "AcceptorProj", "AromaticProj", "PosIonizableProj")

_TET_HALF = np.deg2rad(54.75)   # half the tetrahedral angle (109.47/2)


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-6 else None


def _rotate(v, axis, ang):
    """Rodrigues rotation of v about unit `axis` by `ang` radians."""
    c, s = np.cos(ang), np.sin(ang)
    return v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1 - c)


def _acceptor_lone_pairs(a, xyz):
    """Unit vectors pointing from acceptor atom `a` along its lone pairs.

    Uses *all* neighbours (incl. explicit H) so the geometry is right for
    hydroxyl / ether / carbonyl oxygens and imine / aromatic nitrogens:
      - 2+ neighbours (sp3 O, ether, hydroxyl): two lone pairs straddling the
        plane that bisects the X-A-Y angle, tilted out of the X-A-Y plane by
        the tetrahedral half-angle -> NOT coplanar with the substituents.
      - 1 neighbour (carbonyl O, nitrile N): two lone pairs at +/-120 deg from
        the A=X bond, in the plane defined by X and its other substituents.
    """
    pa = xyz[a.GetIdx()]
    us = [u for u in (_unit(xyz[nb.GetIdx()] - pa) for nb in a.GetNeighbors())
          if u is not None]

    if len(us) >= 2:
        u1, u2 = us[0], us[1]
        bis = _unit(-(u1 + u2))          # in-plane, opposite the bonds
        perp = _unit(np.cross(u1, u2))   # normal to the X-A-Y plane
        if bis is None or perp is None:
            return []
        return [_unit(bis * np.cos(_TET_HALF) + s * perp * np.sin(_TET_HALF))
                for s in (+1.0, -1.0)]

    if len(us) == 1:
        u1 = us[0]
        nb0 = a.GetNeighbors()[0]
        refs = [xyz[x.GetIdx()] for x in nb0.GetNeighbors() if x.GetIdx() != a.GetIdx()]
        if not refs:
            return [(-u1)]              # nothing to define a plane -> straight out
        w = _unit(np.mean(refs, axis=0) - xyz[nb0.GetIdx()])
        axis = _unit(np.cross(u1, w if w is not None else u1 + np.array([1e-3, 0, 0])))
        if axis is None:
            return [(-u1)]
        return [_rotate(-u1, axis, ang) for ang in (np.deg2rad(120), np.deg2rad(-120))]

    return []


class ProjectedPointPharmacophoreGenerator(RDKitPharmacophoreGenerator):

    def __init__(self, donor_dist=3.0, acceptor_dist=3.0, ring_dist=3.5,
                 pos_dist=3.0, proj_weight=1.0, keep_base_points=True,
                 base_families=_BASE, fdefName=None):

        families = list(base_families) + list(_PROJ)
        feats = {f: "rdkit" for f in families}
        interactions = ([(f, f, 1.0, 1.0) for f in base_families]
                        + [(f, f, 1.0, float(proj_weight)) for f in _PROJ])

        super().__init__(features=feats, interactions=interactions, fdefName=fdefName)

        self.donor_dist = donor_dist
        self.acceptor_dist = acceptor_dist
        self.ring_dist = ring_dist
        self.pos_dist = pos_dist
        self.keep_base_points = keep_base_points

    # roshambo2 calls this with an RDKit mol that has a single conformer
    def generate_color_atoms(self, mol):
        base = self.assign_features(mol)                 # [ [family, atom_ids, [x,y,z]], ... ]
        conf = mol.GetConformer()
        xyz = {a.GetIdx(): np.asarray(conf.GetAtomPosition(a.GetIdx()), float)
               for a in mol.GetAtoms()}

        coords, types = [], []

        def add(p, fam):
            coords.append(np.asarray(p, np.float32))
            types.append(self.FEATURES_ENUM[fam])

        for family, atom_ids, pos in base:
            pos = np.asarray(pos, float)
            if self.keep_base_points:
                add(pos, family)

            if family == "Donor":
                for ai in atom_ids:
                    a = mol.GetAtomWithIdx(int(ai))
                    for nb in a.GetNeighbors():
                        if nb.GetAtomicNum() != 1:
                            continue
                        v = xyz[nb.GetIdx()] - xyz[ai]
                        n = np.linalg.norm(v)
                        if n > 1e-3:
                            add(xyz[ai] + v / n * self.donor_dist, "DonorProj")

            elif family == "PosIonizable":
                # cation directionality: along each N-H of a protonated amine
                # (the salt-bridge vector), or straight out if there is no H
                for ai in atom_ids:
                    a = mol.GetAtomWithIdx(int(ai))
                    hs = [nb for nb in a.GetNeighbors() if nb.GetAtomicNum() == 1]
                    if hs:
                        for nb in hs:
                            u = _unit(xyz[nb.GetIdx()] - xyz[ai])
                            if u is not None:
                                add(xyz[ai] + u * self.pos_dist, "PosIonizableProj")
                    else:
                        heavy = [xyz[nb.GetIdx()] for nb in a.GetNeighbors()]
                        u = _unit(xyz[ai] - np.mean(heavy, axis=0)) if heavy else None
                        if u is not None:
                            add(xyz[ai] + u * self.pos_dist, "PosIonizableProj")

            elif family == "Acceptor":
                for ai in atom_ids:
                    a = mol.GetAtomWithIdx(int(ai))
                    for lp in _acceptor_lone_pairs(a, xyz):
                        lp = _unit(lp)
                        if lp is not None:
                            add(xyz[ai] + lp * self.acceptor_dist, "AcceptorProj")

            elif family == "Aromatic":
                ring = np.array([xyz[int(ai)] for ai in atom_ids])
                c = ring.mean(axis=0)
                # ring-plane normal = smallest-variance direction of the ring atoms
                _, _, vt = np.linalg.svd(ring - c)
                normal = vt[2] / np.linalg.norm(vt[2])
                add(c + normal * self.ring_dist, "AromaticProj")
                add(c - normal * self.ring_dist, "AromaticProj")

        if not coords:
            return np.zeros((0, 3), np.float32), np.zeros(0, int)
        return np.asarray(coords, np.float32), np.asarray(types, int)


if __name__ == "__main__":
    # tiny self-test: methanol donor+acceptor, benzene ring
    from rdkit import Chem
    from rdkit.Chem import AllChem

    for smi in ("CO", "c1ccccc1"):
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=1)
        cg = ProjectedPointPharmacophoreGenerator()
        c, t = cg.generate_color_atoms(m)
        idx2fam = cg.get_index_to_feature()
        print(smi, "->", [idx2fam[i] for i in t])
