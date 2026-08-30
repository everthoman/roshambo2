# dopamine vs PD128907 — overlap-only pharmacophore models

Protonated dopamine (query) overlaid on the D3 agonist PD128907 (dataset) with
the CUDA backend, `optim_mode="combination"`. Only the roshambo2 "color" feature
points that **actually overlap** between the two molecules are output/rendered.

| file | what it does |
|---|---|
| `roshambo_dopamine.py` | minimal overlay — prints the score, no files |
| `projected_pharmacophore.py` | `ProjectedPointPharmacophoreGenerator`: adds ROCS-style projected donor / acceptor lone-pair / ring-normal / cation points so colour scoring rewards H-bond & stacking **direction** (importable) |
| `overlaptools.py` | keep only the query↔hit feature pairs whose Gaussian overlap fraction `exp(-d²/2) ≥ min_overlap`; write the model SDF + PyMOL loader (importable) |
| `ligtools.py` | `pose_with_H` (put the *scored* hydrogens on the aligned pose) + `polar_only` (importable) |
| `overlap_pharmacophore.py` | **the script** — samples 25 dopamine confs, keeps the best, builds the overlap-only model for the **stock** and the **projected** colour model, renders both |
| `conf_scan.py` | conformer-count sweep for the projected model |

## Run

```
python overlap_pharmacophore.py            # min_overlap 0.15 (d ≤ ~1.95 Å)
python overlap_pharmacophore.py 0.30       # stricter

pymol overlap_stock.pml                    # 8 shared features, stock model
pymol overlap_projected.pml                # 15 shared features, directional model
```

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
