# dopamine vs PD128907 — overlap-only pharmacophore models

Protonated dopamine (query) overlaid on the D3 agonist PD128907 (dataset) with
the CUDA backend, `optim_mode="combination"`. Only the roshambo2 "color" feature
points that **actually overlap** between the two molecules are output/rendered.

| file | what it does |
|---|---|
| `roshambo_dopamine.py` | minimal overlay — prints the score, no files |
| `projected_pharmacophore.py` | `ProjectedPointPharmacophoreGenerator`: adds ROCS-style projected donor / acceptor lone-pair / ring-normal / cation points so colour scoring rewards H-bond & stacking **direction** (importable) |
| `overlaptools.py` | feature-point extraction, overlap matching (`exp(-d²/2) ≥ min_overlap`), model SDF + PyMOL loader, and `add_tversky` for partial / "fit" scoring (importable) |
| `ligtools.py` | `pose_with_H` (put the *scored* hydrogens on the aligned pose) + `polar_only` (importable) |
| `overlap_pharmacophore.py` | **the script** — samples 25 dopamine confs, keeps the best, builds the overlap-only model for the **stock** and the **projected** colour model, renders both |
| `consensus_pharmacophore.py` | multi-ligand: aligns 6 dopaminergic agonists to the most rigid one (apomorphine), clusters feature points, keeps those seen in ≥ N ligands |
| `conf_scan.py` | conformer-count sweep for the projected model |

## Run

```
python overlap_pharmacophore.py            # min_overlap 0.15 (d ≤ ~1.95 Å)
python overlap_pharmacophore.py 0.30       # stricter
python overlap_pharmacophore.py 0.15 0.95  # rank/select by Tversky (partial match)

pymol overlap_stock.pml                    # 8 shared features, stock model
pymol overlap_projected.pml                # 15 shared features, directional model
```

**Partial matching (Tversky / "fit")** — `overlaptools.add_tversky(df, alpha)`
computes `O_AB / (alpha·O_AA + (1-alpha)·O_BB)` from the overlap + self-overlap
columns roshambo2 already returns. `alpha → 1` scores *how well the query is
covered*, ignoring the hit's extra bulk — so a small query / fragment can match a
sub-region of a bigger molecule without a Tanimoto penalty. Pass it as the 2nd
arg to `overlap_pharmacophore.py`; `consensus_pharmacophore.py` uses `alpha=0.95`
by default so each ligand snaps onto the relevant part of the (rigid) reference.
roshambo2's own optimiser already maximises raw overlap, so it aligns
partial-friendly regardless — Tversky only changes which pose/conformer is
*picked* and *reported*.

Outputs: `overlap_model_{stock,projected}.sdf` (the models),
`overlap_{stock,projected}.{pml,png}` (renders),
`ligand_{dopamine,PD128907}_{stock,projected}.sdf` (context ligands, carrying
their polar N/O hydrogens).

Feature → element → colour: Donor `N` slate · Acceptor `O` firebrick ·
PosIonizable `Fe` orange · NegIonizable `Cl` green · Aromatic `S` yellow ·
Hydrophobe `Br` grey · DonorProj `F` deepblue · AcceptorProj `I` hotpink ·
AromaticProj `B` wheat · PosIonizableProj `P` purple.

## Scores (this machine, ETKDGv3 + MMFF, seed 0xF00D; 25 dopamine / 100 PD128907 confs)

| colour model | shape | colour | combo | shared features (f ≥ 0.15) | % of colour overlap |
|---|---|---|---|---|---|
| stock (1 pt/feature) | 0.633 | 0.500 | 0.566 | 8  | 74 % |
| projected            | 0.629 | 0.362 | 0.496 | 15 | 84 % |

Projected points lower the *combined* colour Tanimoto — dopamine and PD128907
only partly agree on feature *direction*, not just position. The `PosIonizableProj`
salt-bridge vector and one `DonorProj` are the weak links (d ≈ 1.6 Å).

### Consensus of 6 agonists (`consensus_pharmacophore.py`)

Ligands (neutral SMILES) come from `dopamine_ligands.txt`; the basic amine of
each is protonated. roshambo2 aligns only pairwise, so: pick the most rigid
ligand (apomorphine), roshambo2-align the other five to it, then cluster
same-family feature points and keep clusters present in ≥ 4/6 ligands.

13 consensus points — the classic aminergic layout:

- **protonated amine** (`PosIonizable`, 4/6)
- **aromatic ring** — `Aromatic` + ring-normal projected + `Hydrophobe`, 6/6
- **ring hydroxyl** — `Acceptor` / `Donor` and their projected points, 4–5/6

Only `PosIonizableProj` stays below threshold: the amine *position* is shared
but its N–H points different ways across the scaffolds, so the salt-bridge
*vector* does not cluster. (With less accurate input structures even
`PosIonizable` misses — shape+colour overlay pins the large ring systems more
tightly than the amine.)

### Conformer count (`conf_scan.py`, projected model)

| dopamine confs | PD128907 confs | shape | colour | combo |
|---|---|---|---|---|
| 1  | 20  | 0.623 | 0.318 | 0.471 |
| 1  | 100 | 0.629 | 0.333 | 0.481 |
| 25 | 100 | 0.629 | 0.362 | 0.496 |
| 25 | 400 | 0.637 | 0.378 | 0.507 |

PD128907 is rigid — its dataset confs saturate by ~100. Dopamine is the flexible
partner; sampling its query conformers (→25) is what moves the score.

## Notes

- Needs the `roshambo2/roshambo2.py` fix that stops `compute()` from discarding a
  user-supplied `color_generator` (on the `local/megingjord-gpu` branch).
- `pose_with_H` rigid-body-fits the original H-bearing conformer onto roshambo2's
  aligned (H-stripped) skeleton — a plain `Chem.AddHs` guesses its own phenol
  O–H rotamer and the drawn O–H then disagrees with `DonorProj`.
- A phenol / hydroxyl O–H is a low-barrier rotor, so its `DonorProj` direction is
  only as good as the conformer picked. N–H donors (amines, `PosIonizableProj`)
  are fixed by the heavy-atom framework and are reliable.
