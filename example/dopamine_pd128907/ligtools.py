"""Keep the rendered ligand poses consistent with what roshambo2 scored.

roshambo2 rebuilds its aligned output from SMILES with NO hydrogens
(`_smiles_to_3d_noH`). If you just `Chem.AddHs(..., addCoords=True)` that
back, RDKit guesses hydrogen positions from local geometry - and for a
free rotor like a phenol/hydroxyl O-H it will often pick a *different*
rotamer than the ETKDG/MMFF conformer roshambo2 actually used to place the
donor feature. The rendered O-H then disagrees with the DonorProj point.

`pose_with_H` fixes this: it rigid-body-fits the original H-bearing
conformer onto roshambo2's aligned heavy-atom skeleton, so the hydrogens
shown are exactly the ones that were scored.
"""
from rdkit import Chem
from rdkit.Chem import rdMolAlign


def polar_only(mol):
    """Drop C-H hydrogens, keep the ones on N / O (for visualisation)."""
    ed = Chem.RWMol(mol)
    drop = sorted((a.GetIdx() for a in mol.GetAtoms()
                   if a.GetAtomicNum() == 1
                   and a.GetNeighbors()[0].GetAtomicNum() not in (7, 8)), reverse=True)
    for i in drop:
        ed.RemoveAtom(i)
    out = ed.GetMol()
    Chem.SanitizeMol(out)          # fix implicit-H counts left by the removals
    return out


def pose_with_H(aligned_noH, original_withH):
    """Return a copy of `original_withH` (all H kept) rigid-body moved so its
    heavy atoms sit on `aligned_noH`'s heavy atoms. Same molecule, any atom
    order; picks the substructure match with the lowest RMSD (handles ring
    symmetry)."""
    prb = Chem.Mol(original_withH)
    ref_heavy = Chem.RemoveHs(Chem.Mol(original_withH))
    matches = aligned_noH.GetSubstructMatches(ref_heavy, uniquify=False)
    if not matches:                       # fall back to a plain re-hydrogenate
        return Chem.AddHs(Chem.Mol(aligned_noH), addCoords=True)

    best_rms, best = None, None
    for m in matches:
        cand = Chem.Mol(prb)
        atom_map = [(int(i), int(m[i])) for i in range(len(m))]   # (prb heavy, ref)
        rms = rdMolAlign.AlignMol(cand, aligned_noH, atomMap=atom_map, maxIters=200)
        if best_rms is None or rms < best_rms:
            best_rms, best = rms, cand
    return best
