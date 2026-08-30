# PD128907 vs dopamine — overlap-only pharmacophore models

Both molecules contribute **every MMFF conformer within 5 kcal/mol** of their own
global minimum (`overlaptools.conformer_ensemble`, RMS-pruned including the polar
N/O hydrogens so O–H / N–H rotamers survive). roshambo2 aligns every dopamine
conformer to every PD128907 conformer, and each alignment is ranked by

    adjusted = score − ENERGY_LAMBDA · (dE_PD128907 + dE_dopamine)

— the raw shape/colour score minus a linear penalty on the conformational strain
both partners pay. The best-adjusted alignment builds the model. CUDA backend,
`optim_mode="combination"`. Only the "color" feature points that **actually
overlap** are output/rendered.

| file | what it does |
|---|---|
| `roshambo_dopamine.py` | minimal overlay — prints the score, no files |
| `projected_pharmacophore.py` | `ProjectedPointPharmacophoreGenerator`: adds ROCS-style projected donor / acceptor lone-pair / ring-normal / cation points so colour scoring rewards H-bond & stacking **direction** (importable) |
| `overlaptools.py` | `conformer_ensemble` (≤ E-window MMFF confs, polar-H RMS prune), feature-point extraction, overlap matching, model SDF + PyMOL loader, `add_tversky` (importable) |
| `ligtools.py` | `pose_with_H` (put the *scored* hydrogens on the aligned pose) + `polar_only` (importable) |
| `overlap_pharmacophore.py` | **the script** — every ≤ 5 kcal/mol conformer of both molecules, all pairwise alignments, ranked by score − strain penalty; builds the overlap-only model for the **stock** and the **projected** colour model, renders both |
| `consensus_pharmacophore.py` | multi-ligand: aligns 6 dopaminergic agonists to the most rigid one (apomorphine), clusters feature points, keeps those seen in ≥ N ligands |
| `conf_scan.py` | conformer-count sweep for the projected model |

## Run

```
python overlap_pharmacophore.py                 # min_overlap 0.15, strain penalty 0.05/kcal
python overlap_pharmacophore.py 0.30            # stricter feature cutoff
python overlap_pharmacophore.py 0.15 0.05 0.10  # Tversky alpha 0.05, strain penalty 0.10/kcal

pymol overlap_stock.pml                         # 7 shared features, stock model
pymol overlap_projected.pml                     # 11 shared features, directional model
```

Args: `[min_overlap] [tversky_alpha] [energy_lambda]`. `energy_lambda` is score
units per kcal/mol of total strain (default 0.05, so a 2 kcal/mol strained
conformer must score 0.1 higher to be chosen).

**Partial matching (Tversky / "fit")** — `overlaptools.add_tversky(df, alpha)`
computes `O_AB / (alpha·O_query + (1-alpha)·O_hit)` from the overlap + self-overlap
columns roshambo2 already returns.

- `alpha → 1`: *how well is the PD128907 template covered* (rewards a hit that spans it)
- `alpha → 0`: *how well does dopamine fit inside the template* — the "is the small
  molecule a sub-shape of the big one" question; the useful direction here since
  dopamine is the smaller partner

`consensus_pharmacophore.py` uses `alpha=0.95` (the reference apomorphine is the
larger, rigid partner, so "cover the reference" is right there). roshambo2's own
optimiser already maximises raw overlap, so it aligns partial-friendly regardless
— Tversky only changes which pose/conformer is *picked* and *reported*.

Outputs: `overlap_model_{stock,projected}.sdf` (the models),
`overlap_{stock,projected}.{pml,png}` (renders),
`ligand_{PD128907,dopamine}_{stock,projected}.sdf` (context ligands, carrying
their polar N/O hydrogens).

Feature → element → colour: Donor `N` slate · Acceptor `O` firebrick ·
PosIonizable `Fe` orange · NegIonizable `Cl` green · Aromatic `S` yellow ·
Hydrophobe `Br` grey · DonorProj `F` deepblue · AcceptorProj `I` hotpink ·
AromaticProj `B` wheat · PosIonizableProj `P` purple.

## Scores (this machine, ETKDGv3 + MMFF, seed 0xF00D)

≤ 5 kcal/mol ensembles: PD128907 10 conformers (dE 0.0–3.9), dopamine 10
(dE 0.0–4.7). Best-adjusted alignment for both colour models uses each
molecule's **global minimum** (0.0 + 0.0 kcal/mol strain).

| colour model | shape | colour | combo | shared features (f ≥ 0.15) |
|---|---|---|---|---|
| stock (1 pt/feature) | 0.629 | 0.293 | 0.461 | 7  |
| projected            | 0.527 | 0.279 | 0.403 | 11 |

Lower than an unconstrained-conformer search would give: restricted to
physically reasonable geometries, low-energy dopamine cannot span PD128907's
ring **and** its amine at once — the winning pose stacks the aromatic rings and
matches the ring hydroxyl, and the amine (`PosIonizable`) drops out of the
overlap. The strain penalty never bites here because no strained conformer
scores enough better to pay for itself.

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
