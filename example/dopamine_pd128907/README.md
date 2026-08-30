# dopamine vs PD128907 — shape overlay + pharmacophore models

Protonated dopamine (query) overlaid on the D3 agonist PD128907 with the CUDA
backend, `optim_mode="combination"`.

| script | what it does | key outputs |
|---|---|---|
| `roshambo_dopamine.py` | shape/colour overlay, stock 1-point-per-feature colour model | `hits_for_query_dopamine_H+_0.sdf`, `..._color_features.sdf`, `ligand_{dopamine,PD128907}.sdf` |
| `projected_pharmacophore.py` | reusable `ProjectedPointPharmacophoreGenerator` — adds ROCS-style projected donor / acceptor lone-pair / ring-normal points so colour scoring rewards H-bond / stacking **direction** | (importable module) |
| `overlay_projected.py` | runs the overlay with stock vs base+projected vs projected-only colour models and writes the projected feature cloud | `projected_hits_dopamine_*.sdf`, `ligand_*_proj.sdf` |
| `overlap_pharmacophore.py` | using the projected model, keeps only the query/hit feature pairs roshambo2's colour term actually overlaps (Gaussian overlap-fraction cutoff, arg = min fraction, default 0.15) and writes a standalone model | `pharmacophore_overlap_model.sdf`, `ligand_*_ovl.sdf` |
| `conf_scan.py` | scans (dopamine query confs) x (PD128907 dataset confs) for the projected model | prints a table |

## PyMOL

```
pymol view.pml                  # overlay + full per-molecule feature clouds (stock model)
pymol overlay_projected.pml     # overlay + full projected feature cloud
pymol overlap_pharmacophore.pml # overlay + directional overlapping-only model (labelled)
```

Feature → element → colour: Donor `N` slate · Acceptor `O` firebrick ·
PosIonizable `Fe` orange · NegIonizable `Cl` green · Aromatic `S` yellow ·
Hydrophobe `Br` grey · DonorProj `F` · AcceptorProj `I` · AromaticProj `B`.

## Scores (this machine, ETKDGv3 + MMFF, seed 0xF00D)

`overlay_projected.py` — 25 dopamine query confs (best reported), 100 PD128907 dataset confs:

| colour model | shape | colour | combo |
|---|---|---|---|
| stock (1 pt/feature) | 0.633 | 0.500 | 0.566 |
| base + projected pts | 0.656 | 0.360 | 0.508 |
| projected pts only   | 0.503 | 0.512 | 0.507 |

Projected points lower the *combined* colour Tanimoto because the two molecules
only partly agree on feature *direction*, not just position — and "projected
only" will trade shape for a much better directional colour match.

### Conformer count (`conf_scan.py`, base+projected model)

| dopamine confs | PD128907 confs | shape | colour | combo |
|---|---|---|---|---|
| 1  | 20  | 0.625 | 0.337 | 0.481 |
| 1  | 100 | 0.644 | 0.346 | 0.495 |
| 25 | 100 | 0.656 | 0.360 | 0.508 |
| 25 | 400 | 0.642 | 0.391 | 0.516 |

PD128907 is rigid, so its dataset confs saturate by ~100. Dopamine is the
flexible partner — sampling its query conformers (→25) is what actually moves
the score. Plateau ~(25, 100).

## Note

`overlay_projected.py` and `overlap_pharmacophore.py` need the fix in
`roshambo2/roshambo2.py` that stops `compute()` from discarding a user-supplied
`color_generator` (committed on the `local/megingjord-gpu` branch).
