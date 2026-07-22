# SN1 Biolog PM08 — media-optimization screen

Raw plate-reader data and analysis code accompanying the manuscript:

> **_[Manuscript title placeholder]_**
> _[Author list placeholder]_ · _[Journal / preprint placeholder]_ · _[DOI placeholder]_

This repository contains everything needed to reproduce the media-optimization
figures and tables from raw absorbance readings.

---

## What's in this repository

```
.
├── README.md                                      this file
├── LICENSE                                        MIT (code)
├── LICENSE-DATA                                   CC-BY-4.0 (data)
├── CITATION.cff                                   how to cite this dataset/code
├── requirements.txt                               Python dependencies
├── .gitignore
│
├── plot_growth.py                                 raw CSVs -> tidy long table + per-plate grid + interactive viewer
├── analyze_conditions.py                          condition-level ranking, growth kinetics, plots, tables
├── top5_vs_control_590.ipynb                      Jupyter notebook for the top-5-vs-control figure at 590 nm
│                                                  (also contains a companion OD600 cell for the omics follow-up)
│
├── SN1_PlateInfo.csv                              plate UUID -> sample name (SN1_1c, SN1_2c, SN1_3c, SN1_4c)
├── well_assignment.xlsx                           well -> Condition_ID + Replicate + Plate mapping
├── design_concentrations.csv                      per-condition media composition (Phosphates, NH4SO4, CoCl2, Succinate, Methanol, PQQ)
│
├── SN1_RawReads/                                  raw Biolog PM08 reads (4 plates x 96 wells x 2 wavelengths x ~145 timepoints)
│   ├── SN1_PM08_custom_SN1_1c_..._RawReads.csv    Plate 1 raw reads (SN1_1c)
│   ├── SN1_PM08_custom_SN1_2c_..._RawReads.csv    Plate 2 raw reads (SN1_2c)
│   ├── SN1_PM08_custom_SN1_3c_..._RawReads.csv    Plate 3 raw reads (SN1_3c)
│   └── SN1_PM08_custom_SN1_4c_..._RawReads.csv    Plate 4 raw reads (SN1_4c)
│
├── data/omics_growth/                             follow-up OD600 kinetics for the top-5 designs used in the proteomics/RNAseq arm
│   └── 05222026_top5_media_dbtl1.xlsx
│
└── figures/                                       representative outputs (regeneratable — kept for reviewers)
    ├── report_top5_vs_control_740.png             top-5 media vs ctrl_media at 740 nm (turbidity)
    ├── report_top5_vs_control_590.png             top-5 media vs ctrl_media at 590 nm (dye reduction)
    ├── condition_summary_clean.csv                per-condition summary at both wavelengths (AUC, mu_max, doubling time, ...)
    ├── improvement_vs_control_740.csv             top-5 % lift over plate-matched control at 740 nm
    └── improvement_vs_control_590.csv             top-5 % lift over plate-matched control at 590 nm
```

---

## Experimental design (very briefly)

- **Organism.** _Pseudomonas putida_ KT2440 (or equivalent isolate — see manuscript for strain details).
- **Assay.** Biolog PM08 plates (96 wells each), 4 plates in parallel (`SN1_1c`, `SN1_2c`, `SN1_3c`, `SN1_4c`).
- **Reads.** Kinetic absorbance at two wavelengths, **590 nm** (tetrazolium-dye reduction — a metabolic-activity proxy) and **740 nm** (light scattering — a turbidity/biomass proxy), every ~20 min for 48 h at 30 °C.
- **Media designs.** 66 test media (`MPOB_005` … `MPOB_068`) plus the `ctrl_media` control (n = 4 replicate wells per plate, giving n = 16 across all four plates) and calibrant wells (`Nd*`, excluded from analysis).
- **Composition.** Each `MPOB_xxx` design specifies concentrations of six ingredients: Phosphates, NH₄SO₄, CoCl₂, Succinate, Methanol, PQQ. Full compositions in `design_concentrations.csv`.

---

## How to reproduce the analysis

### 1. Environment

```bash
python -m venv .venv
source .venv/bin/activate          # macOS/Linux; Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Python 3.9+ is required. Analysis was originally run on macOS with Python 3.9.

### 2. Regenerate everything from raw CSVs

```bash
python plot_growth.py         # -> outputs/growth_long.csv, outputs/grid_SN1_*.png, outputs/viewer.html
python analyze_conditions.py  # -> outputs/conditions/*.png, *.csv (all condition-level results)
```

The scripts create an `outputs/` folder alongside them. They are idempotent — running
them again overwrites the same files. See the `figures/` folder for reference copies
of the two main manuscript figures and the accompanying CSVs, so reviewers can
check against them without executing anything.

### 3. Notebook (top 5 media at 590 nm)

```bash
jupyter lab top5_vs_control_590.ipynb
```

The notebook is self-contained: it loads `outputs/growth_long.csv` (produced by
`plot_growth.py`) and `well_assignment.xlsx`, ranks conditions by AUC of
Δabs₅₉₀, computes µ_max via a Savitzky–Golay-smoothed sliding-window log-linear
regression, and writes:

- `outputs/conditions/top5_vs_control_590_nb.png` + `.pdf`
- `outputs/conditions/improvement_vs_control_590_nb.csv`

The notebook also has an OD600 companion section that loads
`data/omics_growth/05222026_top5_media_dbtl1.xlsx` (the follow-up growth curves
for the top-5 designs that fed the proteomics/RNAseq measurements) and plots
them alongside `ctrl_media` and the `MP` standard.

---

## Methods (what the code does)

### `plot_growth.py`

- Reads each raw CSV (`Read At`, `Wavelength`, 96 well columns), melts to long form.
- Computes `t_h` = hours since the first read on each plate.
- Joins `SN1_PlateInfo.csv` to map plate UUIDs → sample names (`SN1_1c` etc.).
- Writes `outputs/growth_long.csv` (tidy: `plate, well, row, col, wavelength, t_h, abs`).
- Renders an 8×12 static grid per plate and an interactive Plotly viewer for QC.

### `analyze_conditions.py`

- Joins `growth_long.csv` with `well_assignment.xlsx` on `(plate, well)`.
- Excludes `Nd*` calibrant wells.
- **Per-well kinetic features** (function `features_one_well`):
  - **AUC** — trapezoidal integral of Δabs (baseline-subtracted vs t=0) over the full time course. Units: abs·h.
  - **max Δabs** — highest Δabs reached.
  - **µ_max** — specific growth rate, computed by `calculate_mu_max`:
    1. Baseline-subtract each well against its own t=0.
    2. Savitzky–Golay smooth (window = 7 points, polyorder = 2).
    3. Slide a **4-h window** across the trace; linear-regress `ln(Δabs_smoothed)` on `t`.
    4. Reject any window unless: every point has Δabs ≥ 0.015 (well above the noise floor), ≥ 6 points in the window, and R² ≥ 0.95.
    5. Take the maximum accepted slope. Doubling time = `ln(2)/µ_max`.
  - **lag** — time at which Δabs first crosses 10 % of its maximum.
- **Aggregation** — replicate wells (typically n = 4) → condition mean ± std, at each wavelength.
- **Plots** — ranked bar charts, growth-rate rankings, per-plate kinetics, top-N vs control overlay, plate-plate concordance, z-scored metric heatmap, standalone HTML viewer with checkbox condition picker.
- **Statistics** — one-way ANOVA across conditions per metric per wavelength; % lift computed vs **plate-matched** control (each condition compared to `ctrl_media` on the same plate, absorbing plate effects).

### `top5_vs_control_590.ipynb`

Focused notebook that reproduces the 590 nm top-5-vs-control comparison and
computes both AUC and µ_max for the top 5 designs with an accompanying
improvement table. The µ_max implementation is identical to the routine in
`analyze_conditions.py` and its output matches to 3 decimal places (verified
against `figures/condition_summary_clean.csv`).

---

## Key results (top 5 media at 590 nm)

Reproduced from `figures/improvement_vs_control_590.csv`:

| Condition | n | Plate | AUC₅₉₀ (mean ± std) | AUC vs ctrl | µ_max (h⁻¹) | µ_max vs ctrl | t_d (h) |
|-----------|---|-------|---------------------|-------------|-------------|----------------|---------|
| MPOB_008 | 4 | SN1_1c | 7.066 ± 0.252 | +70.8 % | 0.279 ± 0.007 | +29.6 % | 2.48 |
| MPOB_058 | 4 | SN1_4c | 6.948 ± 0.294 | +76.4 % | 0.281 ± 0.008 | +38.9 % | 2.47 |
| MPOB_041 | 4 | SN1_2c | 5.748 ± 0.255 | +47.5 % | 0.296 ± 0.012 | +49.9 % | 2.34 |
| MPOB_045 | 4 | SN1_2c | 5.620 ± 0.241 | +44.3 % | 0.229 ± 0.014 | +16.2 % | 3.02 |
| MPOB_068 | 4 | SN1_2c | 5.403 ± 0.144 | +38.7 % | 0.239 ± 0.007 | +21.4 % | 2.90 |

All comparisons are against **plate-matched** `ctrl_media` (mean of the 4
control wells on the same plate as the test condition), which absorbs
plate-to-plate variation.

---

## Data-format reference

### Raw CSV (one file per plate, in `SN1_RawReads/`)

Columns: `PlateId, Wavelength, Read At, A01, A02, …, H12` (96 well columns).
Each row is a single read of the whole plate at one wavelength. `PlateId` is a
UUID that maps to a sample name via `SN1_PlateInfo.csv`.

### `well_assignment.xlsx`

One row per (Plate, Well). Columns: `Plate` (1–4), `Well` (`A1` … `H12`, no
zero-pad), `Condition_ID` (e.g. `MPOB_008`, `ctrl_media`, `Nd*`), `Replicate`
(1–4). The scripts zero-pad wells (`A1` → `A01`) and map `Plate` → sample name
via `PLATE_MAP` (`1: SN1_1c, 2: SN1_2c, 3: SN1_3c, 4: SN1_4c`).

### `design_concentrations.csv`

One row per (Condition_ID, Plate, Well, Replicate). Columns:
`Phosphates, NH4SO4, CoCl2, Succinate, Methanol, PQQ` (concentrations in the
units used in the manuscript — see paper for details).

### `SN1_PlateInfo.csv`

Plate-level metadata: UUID, sample name, incubation temperature (30 °C),
incubation time (48 h), cycle time (20 min), etc. Written by the Biolog reader.

---

## Software versions

Reference environment:

- Python 3.9
- pandas ≥ 1.5
- numpy ≥ 1.23 (either NumPy 1.x with `np.trapz` or NumPy 2.x with `np.trapezoid` works)
- scipy ≥ 1.9
- matplotlib ≥ 3.6
- plotly ≥ 5.14
- openpyxl (for reading `.xlsx`)
- jupyterlab (optional, for the notebook)

See `requirements.txt` for pinnable versions.

---

## License

- **Code** (`*.py`, `*.ipynb`) — MIT License, see [`LICENSE`](LICENSE).
- **Data** (`SN1_RawReads/*`, `SN1_PlateInfo.csv`, `well_assignment.xlsx`, `design_concentrations.csv`, `data/omics_growth/*`, all files in `figures/`) — Creative Commons Attribution 4.0 International (CC-BY-4.0), see [`LICENSE-DATA`](LICENSE-DATA).

If you use this dataset or code, please cite the manuscript and this
repository (see `CITATION.cff`).

---

## Contact

Questions about the data or analysis: _[Author name]_ · _[email placeholder]_ · Lawrence Berkeley National Laboratory
