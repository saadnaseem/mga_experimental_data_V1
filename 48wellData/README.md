# 48-well OD600 kinetics — top-5 media follow-up (DBTL round 1)

Follow-up growth curves for the top-5 media designs identified in the SN1
Biolog PM08 screen (`../figures/condition_summary_clean.csv` in the parent
repository), run in a 48-well plate at OD600 to feed the omics (proteomics /
transcriptomics) round.

Everything in this folder is self-contained — drop it anywhere, run
`python clean_omics.py`, and the cleaned CSVs regenerate from the xlsx.

## Files

| File | What it is |
|------|------------|
| `05222026_top5_media_dbtl1.xlsx` | Raw plate-reader export (unchanged from the instrument). |
| `clean_omics.py` | Script that produces the four CSVs below from the xlsx. |
| `plate_map.csv` | Well → condition/replicate/role mapping for all 48 wells. |
| `omics_growth_long.csv` | Per-well tidy long-format kinetics for the 21 sample wells (raw + blank-subtracted OD600). |
| `omics_growth_by_condition.csv` | Aggregated per-condition per-timepoint: mean ± std of blank-subtracted OD600 across replicates. |
| `omics_growth_summary_per_condition.csv` | One row per condition: OD_final, OD_max, µ_max, doubling time (mean ± std across 3 replicates). |

## Experiment summary

- **Organism:** _Methylobacterium extorquens_ AM1 (or equivalent — see manuscript).
- **Read:** OD600 in a plate reader, 20-min cycle, 30 °C.
- **Time course:** 145 timepoints, ~48 h total.
- **Plate:** single 48-well plate (6 rows A–F × 8 cols 1–8).

## Plate layout (all 48 wells)

```
     1     2      3      4      5      6      7      8
A    -     -      -      -      -      -      -      -       (unused edge row)
B  blank  008    008    008    058    058    058   blank
C  blank  045    045    045    041    041    041   blank
D  blank  068    068    068  ctrl_m  ctrl_m ctrl_m blank
E  blank   MP     MP     MP     -      -      -    blank
F    -     -      -      -      -      -      -      -       (unused edge row)
```

- **7 conditions × 3 replicates = 21 sample wells** + **8 blank wells** (media
  only, no cells) + **19 unused wells** (rows A/F and E5–E7 left empty to reduce
  evaporation artefacts).
- `MPOB_008 / 058 / 045 / 041 / 068` are the top-5 media designs from the SN1
  Biolog PM08 screen (see `../figures/condition_summary_clean.csv` for their
  ranking and lift over the in-plate control at 590 / 740 nm).
- `ctrl_media` is the same base-media control used across the SN1 Biolog screen.
- `MP` is the standard _M. extorquens_ minimal medium — an external reference
  for the follow-up.

## How the data was processed

1. **Load raw xlsx** — sheet `Sheet1` has a header row (`Time`, `T° abs600:600`,
   then 48 well names A1..F8), followed by 145 timepoint rows. `Time` is a
   wall-clock value (mixed `datetime.time` for the first day and
   `datetime.timedelta` for later readings) — converted to hours since assay
   start.
2. **Blank series** — at each timepoint, mean of the 8 blank wells
   (B1, C1, D1, E1, B8, C8, D8, E8). Captures both media background and
   plate-reader drift.
3. **Blank-subtract** per well:
   `OD600_blank_subtracted(t) = OD600_raw(well, t) − blank_mean(t)`.
4. **Aggregate replicates** to (mean, std, n) per condition per timepoint
   (`omics_growth_by_condition.csv`).
5. **Per-condition kinetics** (`omics_growth_summary_per_condition.csv`):
   - `OD_final` — blank-subtracted OD600 at the last timepoint.
   - `OD_max` — blank-subtracted OD600 peak across the whole trace.
   - `µ_max` — specific growth rate (h⁻¹), fit by the SAME routine as the SN1
     Biolog analysis: Savitzky–Golay smooth (window 7, poly 2) of raw OD minus
     its own t=0 baseline, then sliding 4-h window log-linear regression on
     `ln(Δabs)`, accepting windows where every point has Δabs ≥ 0.015, ≥ 6
     points in window, R² ≥ 0.95. Max accepted slope = µ_max.
   - `doubling_time_h = ln(2) / µ_max`.

## Column reference — `omics_growth_long.csv`

| column | description |
|--------|-------------|
| `time_h` | Hours since assay start. |
| `well` | 48-well plate coordinate (e.g. `B3`). |
| `row`, `col` | Row letter (A–F) and column number (1–8). |
| `condition` | One of `MPOB_008 / 058 / 045 / 041 / 068 / ctrl_media / MP`. |
| `replicate` | 1, 2, or 3 within a condition. |
| `OD600_raw` | Raw reader value (unitless absorbance at 600 nm). |
| `blank_mean` | Mean OD600 across the 8 blank wells at the same timepoint. |
| `OD600_blank_subtracted` | `OD600_raw − blank_mean`. Use this for growth curves. |
| `temperature_C` | Recorded plate temperature. |

## Per-condition summary (reproduced from `omics_growth_summary_per_condition.csv`)

| Condition | n | Wells | OD_final | OD_max | µ_max (h⁻¹) | Doubling time (h) |
|-----------|---|-------|----------|--------|-------------|-------------------|
| MPOB_008 | 3 | B2–B4 | 0.694 ± 0.005 | 0.695 ± 0.004 | 0.491 ± 0.025 | 1.41 |
| MPOB_068 | 3 | D2–D4 | 0.689 ± 0.006 | 0.706 ± 0.013 | 0.487 ± 0.085 | 1.42 |
| MPOB_045 | 3 | C2–C4 | 0.654 ± 0.018 | 0.665 ± 0.017 | 0.360 ± 0.086 | 1.93 |
| ctrl_media | 3 | D5–D7 | 0.497 ± 0.021 | 0.671 ± 0.172 | 0.434 ± 0.034 | 1.60 |
| MPOB_058 | 3 | B5–B7 | 0.320 ± 0.002 | 0.617 ± 0.018 | 0.356 ± 0.096 | 1.95 |
| MP | 3 | E2–E4 | 0.293 ± 0.098 | 0.409 ± 0.087 | 0.364 ± 0.013 | 1.90 |
| MPOB_041 | 3 | C5–C7 | 0.286 ± 0.002 | 0.449 ± 0.005 | 0.356 ± 0.126 | 1.95 |

Note the OD_max vs OD_final gap for MPOB_058 and MPOB_041 (peak > final): the
culture hit a maximum, then declined (nutrient depletion / lysis). MPOB_008
and MPOB_068 stayed flat at their peak.

## Reproducing the aggregates

```bash
cd 48wellData
python clean_omics.py
```

Regenerates `plate_map.csv`, `omics_growth_long.csv`, `omics_growth_by_condition.csv`,
and `omics_growth_summary_per_condition.csv` from the xlsx. Idempotent.

Dependencies (same as the parent repo): `numpy`, `pandas`, `scipy`, `openpyxl`
— see `../requirements.txt`.

## Provenance

- Instrument export: `05222026_top5_media_dbtl1.xlsx`
  (also present in the parent repo at `../data/omics_growth/` because
  `../top5_vs_control_590.ipynb` references it there).
- Cleaning code: `clean_omics.py` in this folder.
