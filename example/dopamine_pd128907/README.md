# dopamine vs PD128907 — shape overlay + pharmacophore models

Protonated dopamine (query) overlaid on the D3 agonist PD128907 with the CUDA
backend, `optim_mode="combination"`.

| script | what it does | key outputs |
|---|---|---|
| `roshambo_dopamine.py` | shape/colour overlay, stock 1-point-per-feature colour model | `hits_for_query_dopamine_H+_0.sdf`, `..._color_features.sdf`, `ligand_{dopamine,PD128907}.sdf` |
| `overlap_pharmacophore.py` | keeps only the query/hit feature pairs roshambo2's colour term actually overlaps (Gaussian overlap-fraction cutoff, arg = min fraction, default 0.10) and writes a standalone model | `pharmacophore_overlap_model.sdf` |
| `projected_pharmacophore.py` | reusable `ProjectedPointPharmacophoreGenerator` — adds ROCS-style projected donor / acceptor / ring-normal points so colour scoring rewards H-bond / stacking **direction** | (importable module) |
| `overlay_projected.py` | runs the overlay with stock vs base+projected vs projected-only colour models and writes the projected feature cloud | `projected_hits_dopamine_*.sdf`, `ligand_*_proj.sdf` |
| `conf_scan.py` | scans (dopamine query confs) x (PD128907 dataset confs) for the projected model | prints a table |

## PyMOL

```
pymol view.pml                  # overlay + full per-molecule feature clouds
pymol overlap_pharmacophore.pml # overlay + overlapping-only model (labelled)
pymol overlay_projected.pml     # overlay + projected feature cloud
```

Feature → element → colour: Donor `N` slate · Acceptor `O` firebrick ·
PosIonizable `Fe` orange · NegIonizable `Cl` green · Aromatic `S` yellow ·
Hydrophobe `Br` grey · DonorProj `F` · AcceptorProj `I` · AromaticProj `B`.

## Scores (this machine, ETKDGv3 + MMFF, seed 0xF00D)

| colour model | shape | colour | combo |
|---|---|---|---|
| stock (1 pt/feature) | 0.627 | 0.478 | 0.552 |
| base + projected pts | 0.634 | 0.355 | 0.494 |
| projected pts only   | 0.618 | 0.280 | 0.449 |

Projected points lower the colour Tanimoto because the two molecules only
partly agree on feature *direction*, not just position.

### Conformer count (`conf_scan.py`, base+projected model)

| dopamine confs | PD128907 confs | shape | colour | combo |
|---|---|---|---|---|
| 1  | 20  | 0.634 | 0.355 | 0.494 |
| 1  | 100 | 0.620 | 0.397 | 0.508 |
| 25 | 100 | 0.639 | 0.411 | 0.525 |
| 25 | 400 | 0.639 | 0.419 | 0.529 |

PD128907 is rigid, so its dataset confs saturate by ~100. Dopamine is the
flexible partner — sampling its query conformers (→25) is what actually moves
the score. Plateau ~(25, 100); combo 0.49 → 0.53.

## Note

`overlay_projected.py` needs the fix in `roshambo2/roshambo2.py` that stops
`compute()` from discarding a user-supplied `color_generator` (committed on the
`local/megingjord-gpu` branch).
