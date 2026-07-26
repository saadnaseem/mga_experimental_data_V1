"""Pairwise response-surface grids for the six designed factors of the SN1
Biolog PM08 screen, at 590 nm.

Produces two combined 6×6 grids in the reference-figure style — one metric in
the upper triangle, another in the lower triangle, with independent colorbars:

    response_surface_grid_Amax_upper_mumax_lower_590.png    Amax  (upper) + µ_max (lower)
    response_surface_grid_AUC_upper_mumax_lower_590.png     AUC   (upper) + µ_max (lower)

Also emits single-metric 6×6 grids (symmetric — for detailed per-metric
inspection) and the flat per-condition table used to fit all surfaces.

Method:
1. Joins raw kinetics (outputs/growth_long.csv) with design_concentrations.csv
   and computes per-well responses (n = 4 replicates per condition).
2. Aggregates to one row per Condition_ID (66 MPOB designs after excluding
   ctrl_media and Nd calibrants) with 6 factors + 3 responses (AUC, Amax, µ_max).
3. Fits a global thin-plate-spline RBF interpolator in the 6-D factor space
   (factor values z-standardised for numerical stability), with a small
   smoothing term so measurement noise doesn't drive the wiggles.
4. For every ordered pair (row_factor, col_factor), predicts the metric on a
   40 × 40 grid of (col, row) values while holding the OTHER four factors at
   the reference condition's values (default: MPOB_008). This is the classical
   "response surface slice" convention (matches the reference figure's use of
   the best condition as anchor).
5. Draws contourf panels with observed dots (colour = measured value), a star
   at each panel's predicted argmax, and a filled square at the reference.

Outputs (release/figures/):
    response_surface_grid_Amax_upper_mumax_lower_590.{png,pdf}   combined figure
    response_surface_grid_AUC_upper_mumax_lower_590.{png,pdf}    combined figure
    response_surface_AUC_590.{png,pdf}       single-metric symmetric grid
    response_surface_Amax_590.{png,pdf}      single-metric symmetric grid
    response_surface_mu_max_590.{png,pdf}    single-metric symmetric grid
    per_condition_responses_590.csv          66 rows, 6 factors + 3 responses

Run:
    python response_surfaces_590.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.interpolate import RBFInterpolator
from scipy.signal import savgol_filter

# ---------- config ----------
HERE = Path(__file__).resolve().parent
GROWTH_CSV = HERE / 'outputs' / 'growth_long.csv'
DESIGN_CSV = HERE / 'design_concentrations.csv'
FIG_DIR = HERE / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

WAVELENGTH = 590
PLATE_MAP = {1: 'SN1_1c', 2: 'SN1_2c', 3: 'SN1_3c', 4: 'SN1_4c'}
REFERENCE = 'MPOB_008'      # anchor point for holding "other" factors fixed
GRID_N = 40                 # grid resolution per surface
EXCLUDE_PREFIXES = ('Nd', 'ctrl')

FACTORS = ['Phosphates', 'NH4SO4', 'CoCl2', 'Succinate', 'Methanol', 'PQQ']
FACTOR_UNITS = {'Phosphates': 'mM', 'NH4SO4': 'mM', 'CoCl2': 'µM',
                'Succinate': 'mM', 'Methanol': 'mM', 'PQQ': 'µM'}


# ---------- per-well kinetics (identical routine to analyze_conditions.py) ----------
def calculate_mu_max(t, y, win_h=4.0, min_r2=0.95, min_delta_for_fit=0.015,
                     min_window_points=6, smooth_window=7, smooth_polyorder=2,
                     epsilon=1e-4):
    order = np.argsort(t); t = t[order]; y = y[order]
    if len(t) < max(smooth_window, 4):
        return np.nan
    dy = y - y[0]
    win = min(smooth_window, len(dy) - (1 - len(dy) % 2))
    if win < smooth_polyorder + 2:
        dy_s = dy
    else:
        if win % 2 == 0: win -= 1
        dy_s = savgol_filter(dy, win, smooth_polyorder)
    if dy_s.max() < min_delta_for_fit:
        return np.nan
    log_y = np.log(np.clip(dy_s, epsilon, None))
    best = -np.inf
    for i in range(len(t)):
        mask = (t >= t[i]) & (t <= t[i] + win_h)
        if mask.sum() < min_window_points: continue
        if dy_s[mask].min() < min_delta_for_fit: continue
        tt = t[mask]; yy = log_y[mask]
        slope, intercept = np.polyfit(tt, yy, 1)
        ss_res = float(((yy - (slope*tt+intercept))**2).sum())
        ss_tot = float(((yy - yy.mean())**2).sum())
        if ss_tot < 1e-12: continue
        r2 = 1.0 - ss_res / ss_tot
        if r2 < min_r2 or slope <= 0: continue
        if slope > best: best = float(slope)
    return best if np.isfinite(best) else np.nan


def features_one_well(t, y):
    order = np.argsort(t); t = t[order]; y = y[order]
    dy = y - y[0]
    trapz = getattr(np, 'trapezoid', np.trapz)
    return dict(AUC=float(trapz(dy, t)), Amax=float(dy.max()),
                mu_max=calculate_mu_max(t, y))


# ---------- load + join + per-condition table ----------
def build_per_condition_table():
    growth = pd.read_csv(GROWTH_CSV)
    growth = growth[growth['wavelength'] == WAVELENGTH].copy()

    design = pd.read_csv(DESIGN_CSV)
    design['plate'] = design['Plate'].map(PLATE_MAP)
    design['well'] = design['Well'].map(lambda w: f'{w[0].upper()}{int(w[1:]):02d}')

    joined = growth.merge(
        design[['plate', 'well', 'Condition_ID', 'Replicate'] + FACTORS],
        on=['plate', 'well'], how='inner'
    )

    rows = []
    for (cond, plate, well, rep), g in joined.groupby(
            ['Condition_ID', 'plate', 'well', 'Replicate']):
        feats = features_one_well(g['t_h'].to_numpy(), g['abs'].to_numpy())
        row = {'Condition_ID': cond, 'plate': plate, 'well': well, 'rep': rep}
        for f in FACTORS:
            row[f] = float(g[f].iloc[0])
        row.update(feats)
        rows.append(row)
    per_well = pd.DataFrame(rows)

    per_cond = (per_well.groupby('Condition_ID')
                        .agg(**{f: (f, 'first') for f in FACTORS},
                             AUC=('AUC', 'mean'),
                             AUC_std=('AUC', 'std'),
                             Amax=('Amax', 'mean'),
                             Amax_std=('Amax', 'std'),
                             mu_max=('mu_max', 'mean'),
                             mu_max_std=('mu_max', 'std'),
                             n_reps=('AUC', 'count'))
                        .reset_index())
    # keep only real MPOB designs
    mask = ~per_cond['Condition_ID'].str.startswith(EXCLUDE_PREFIXES)
    per_cond = per_cond.loc[mask].dropna(subset=['AUC', 'Amax', 'mu_max']).reset_index(drop=True)
    return per_cond


# ---------- surface fitting (used by both single and combined plotters) ----------
def fit_surfaces(data: pd.DataFrame, metric: str, ref_vals: np.ndarray):
    """Fit a global thin-plate-spline RBF for `metric` on the 6-D factor space
    and predict a 40×40 grid for every ordered (row_factor, col_factor) pair,
    with the other 4 factors held at `ref_vals`.

    Returns:
      grids[(i, j)] -> (XI, XJ, Z, xi_grid, xj_grid)
      argmax_pt[(i, j)] -> (x_j_star, x_i_star)
      vmin, vmax    -> colour range covering both observed and predicted values
    """
    X = data[FACTORS].to_numpy(dtype=float)
    y = data[metric].to_numpy(dtype=float)
    X_mean = X.mean(axis=0); X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0
    Xn = (X - X_mean) / X_std
    smoothing = max(1e-6, 0.05 * float(np.nanstd(y)))
    interp = RBFInterpolator(Xn, y, kernel='thin_plate_spline', smoothing=smoothing)

    grids: dict = {}
    argmax_pt: dict = {}
    for i, fi in enumerate(FACTORS):
        for j, fj in enumerate(FACTORS):
            if i == j: continue
            xi_grid = np.linspace(data[fi].min(), data[fi].max(), GRID_N)
            xj_grid = np.linspace(data[fj].min(), data[fj].max(), GRID_N)
            XJ, XI = np.meshgrid(xj_grid, xi_grid)
            pts = np.tile(ref_vals, (XI.size, 1))
            pts[:, i] = XI.ravel()
            pts[:, j] = XJ.ravel()
            Z = interp((pts - X_mean) / X_std).reshape(XI.shape)
            grids[(i, j)] = (XI, XJ, Z, xi_grid, xj_grid)
            k = int(np.argmax(Z))
            argmax_pt[(i, j)] = (XJ.ravel()[k], XI.ravel()[k])

    z_all = np.concatenate([g[2].ravel() for g in grids.values()])
    obs_min, obs_max = float(np.nanmin(y)), float(np.nanmax(y))
    vmin = min(obs_min, float(np.nanpercentile(z_all, 2)))
    vmax = max(obs_max, float(np.nanpercentile(z_all, 98)))
    return grids, argmax_pt, vmin, vmax, obs_min, obs_max, float(np.nanmax(z_all))


def _draw_panel(ax, XI, XJ, Z, argmax, data, metric, ref_vals,
                i, j, vmin, vmax, cmap):
    from matplotlib.ticker import MaxNLocator
    levels = np.linspace(vmin, vmax, 20)
    ax.contourf(XJ, XI, Z, levels=levels, cmap=cmap,
                vmin=vmin, vmax=vmax, extend='both')
    ax.scatter(data[FACTORS[j]], data[FACTORS[i]], c=data[metric],
               cmap=cmap, vmin=vmin, vmax=vmax,
               s=14, edgecolor='black', lw=0.35, alpha=0.9)
    ax.plot(argmax[0], argmax[1], marker='*', color='yellow',
            markeredgecolor='black', markersize=14, markeredgewidth=0.7, zorder=6)
    ax.plot(ref_vals[j], ref_vals[i], marker='s', color='black',
            markersize=7, markeredgecolor='white', markeredgewidth=0.7, zorder=6)
    # concentration ticks on EVERY panel so the range is readable without
    # eye-tracking to the outer edge
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))
    ax.tick_params(axis='x', labelsize=6.5, rotation=0, pad=1)
    ax.tick_params(axis='y', labelsize=6.5, pad=1)


def _draw_diagonal(ax, name, unit):
    ax.set_facecolor('#f7f7f7')
    ax.text(0.5, 0.55, name, ha='center', va='center',
            fontsize=13, fontweight='bold', color='#333')
    ax.text(0.5, 0.35, f'({unit})', ha='center', va='center', fontsize=10, color='#666')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ('top', 'right', 'bottom', 'left'):
        ax.spines[s].set_visible(False)


def _get_ref_vals(data: pd.DataFrame, fallback_metric: str) -> np.ndarray:
    ref_row = data[data['Condition_ID'] == REFERENCE]
    if ref_row.empty:
        ref_row = data.iloc[[data[fallback_metric].idxmax()]]
    return ref_row[FACTORS].iloc[0].to_numpy(dtype=float)


# ---------- single-metric symmetric grid ----------
def plot_single_metric_grid(data: pd.DataFrame, metric: str, metric_label: str,
                             cbar_label: str, out_png: Path, out_pdf: Path,
                             cmap: str = 'viridis'):
    ref_vals = _get_ref_vals(data, metric)
    grids, argmax, vmin, vmax, obs_min, obs_max, argmax_pred = fit_surfaces(data, metric, ref_vals)

    fig, axes = plt.subplots(6, 6, figsize=(17, 17))
    fig.suptitle(f'Pairwise response-surface grid — {metric_label} @ {WAVELENGTH} nm\n'
                 f'"Other" factors held at {REFERENCE}; dots = observed conditions '
                 f'(colour = measured {metric_label}); ★ = surface argmax; ■ = {REFERENCE}',
                 fontsize=13, y=0.995)

    for i in range(6):
        for j in range(6):
            ax = axes[i, j]
            ax.tick_params(axis='both', labelsize=7)
            if i == j:
                _draw_diagonal(ax, FACTORS[i], FACTOR_UNITS[FACTORS[i]])
                continue
            XI, XJ, Z, _, _ = grids[(i, j)]
            _draw_panel(ax, XI, XJ, Z, argmax[(i, j)], data, metric, ref_vals,
                        i, j, vmin, vmax, cmap)
            # outer axis labels (bold factor name); inner panels keep numeric ticks only
            if i == 5:
                ax.set_xlabel(f'{FACTORS[j]} ({FACTOR_UNITS[FACTORS[j]]})', fontsize=9)
            if j == 0:
                ax.set_ylabel(f'{FACTORS[i]} ({FACTOR_UNITS[FACTORS[i]]})', fontsize=9)

    fig.subplots_adjust(left=0.06, right=0.90, top=0.94, bottom=0.06,
                        wspace=0.32, hspace=0.32)
    cax = fig.add_axes([0.915, 0.15, 0.014, 0.7])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(cbar_label, fontsize=11)
    cb.ax.tick_params(labelsize=9)

    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)
    return {'vmin': vmin, 'vmax': vmax,
            'obs_range': (obs_min, obs_max),
            'argmax_predicted': argmax_pred}


# ---------- combined-triangles grid (upper / lower show different metrics) ----------
def plot_combined_triangles_grid(data: pd.DataFrame,
                                  upper_metric: str, upper_label: str, upper_cbar: str,
                                  lower_metric: str, lower_label: str, lower_cbar: str,
                                  out_png: Path, out_pdf: Path,
                                  upper_cmap: str = 'viridis',
                                  lower_cmap: str = 'plasma'):
    """Reference-figure layout: one metric per triangle, two colorbars, no
    redundant panels."""
    ref_vals = _get_ref_vals(data, upper_metric)
    up_grids, up_argmax, up_vmin, up_vmax, *_ = fit_surfaces(data, upper_metric, ref_vals)
    lo_grids, lo_argmax, lo_vmin, lo_vmax, *_ = fit_surfaces(data, lower_metric, ref_vals)

    fig, axes = plt.subplots(6, 6, figsize=(17, 17))
    fig.suptitle(
        f'Pairwise response-surface grid @ {WAVELENGTH} nm — '
        f'{upper_label} (upper triangle) vs {lower_label} (lower triangle)\n'
        f'"Other" factors held at {REFERENCE}; dots = observed conditions '
        f'(colour = measured metric of that triangle); ★ = surface argmax; ■ = {REFERENCE}',
        fontsize=13, y=0.995
    )

    for i in range(6):
        for j in range(6):
            ax = axes[i, j]
            ax.tick_params(axis='both', labelsize=7)
            if i == j:
                _draw_diagonal(ax, FACTORS[i], FACTOR_UNITS[FACTORS[i]])
                continue
            if j > i:      # upper triangle
                XI, XJ, Z, _, _ = up_grids[(i, j)]
                _draw_panel(ax, XI, XJ, Z, up_argmax[(i, j)], data, upper_metric,
                            ref_vals, i, j, up_vmin, up_vmax, upper_cmap)
            else:          # lower triangle
                XI, XJ, Z, _, _ = lo_grids[(i, j)]
                _draw_panel(ax, XI, XJ, Z, lo_argmax[(i, j)], data, lower_metric,
                            ref_vals, i, j, lo_vmin, lo_vmax, lower_cmap)

            # outer axis labels (bold factor name); inner panels keep numeric ticks only
            if i == 5:
                ax.set_xlabel(f'{FACTORS[j]} ({FACTOR_UNITS[FACTORS[j]]})', fontsize=9)
            if j == 0:
                ax.set_ylabel(f'{FACTORS[i]} ({FACTOR_UNITS[FACTORS[i]]})', fontsize=9)

    fig.subplots_adjust(left=0.06, right=0.88, top=0.93, bottom=0.06,
                        wspace=0.32, hspace=0.32)
    # two colorbars: upper on the right (upper half), lower under it
    cax_up = fig.add_axes([0.90, 0.52, 0.014, 0.36])
    sm_up = plt.cm.ScalarMappable(cmap=upper_cmap, norm=plt.Normalize(vmin=up_vmin, vmax=up_vmax))
    sm_up.set_array([])
    cb_up = fig.colorbar(sm_up, cax=cax_up)
    cb_up.set_label(f'{upper_cbar} (upper △)', fontsize=10)
    cb_up.ax.tick_params(labelsize=8)

    cax_lo = fig.add_axes([0.90, 0.10, 0.014, 0.36])
    sm_lo = plt.cm.ScalarMappable(cmap=lower_cmap, norm=plt.Normalize(vmin=lo_vmin, vmax=lo_vmax))
    sm_lo.set_array([])
    cb_lo = fig.colorbar(sm_lo, cax=cax_lo)
    cb_lo.set_label(f'{lower_cbar} (lower △)', fontsize=10)
    cb_lo.ax.tick_params(labelsize=8)

    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)
    return {'upper': (up_vmin, up_vmax), 'lower': (lo_vmin, lo_vmax)}


# ---------- main ----------
def main():
    print(f'Loading kinetics @ {WAVELENGTH} nm and design factors …')
    per_cond = build_per_condition_table()
    per_cond.to_csv(FIG_DIR / f'per_condition_responses_{WAVELENGTH}.csv', index=False)
    print(f'  {len(per_cond)} conditions × 6 factors × 3 responses')
    print(f'  reference for slicing = {REFERENCE}')
    if REFERENCE not in per_cond['Condition_ID'].values:
        print(f'  ! {REFERENCE} not found — using argmax fallback')

    specs = [
        ('AUC',    'AUC (integrated Δabs)', f'AUC$_{{{WAVELENGTH}\\,\\mathrm{{nm}}}}$  (abs·h)'),
        ('Amax',   'A$_{max}$ (peak Δabs)',  f'A$_{{max,{WAVELENGTH}\\,\\mathrm{{nm}}}}$  (abs)'),
        ('mu_max', 'µ$_{max}$ (specific growth rate)',
                                             f'µ$_{{max,{WAVELENGTH}\\,\\mathrm{{nm}}}}$  (h$^{{-1}}$)'),
    ]

    # ---- single-metric symmetric grids (for detailed per-metric inspection) ----
    print('\n[1/2] single-metric grids (symmetric)')
    for metric, label, cb in specs:
        out_png = FIG_DIR / f'response_surface_{metric}_{WAVELENGTH}.png'
        out_pdf = FIG_DIR / f'response_surface_{metric}_{WAVELENGTH}.pdf'
        stats = plot_single_metric_grid(per_cond, metric, label, cb, out_png, out_pdf)
        print(f'  {metric:<7}  obs range = {stats["obs_range"][0]:.3f} … {stats["obs_range"][1]:.3f}, '
              f'predicted argmax across slices = {stats["argmax_predicted"]:.3f}')
        print(f'           -> {out_png.name} + .pdf')

    # ---- combined-triangle grids (one metric per triangle, no redundant panels) ----
    print('\n[2/2] combined-triangle grids (upper vs lower)')
    combined_specs = [
        ('Amax_upper_mumax_lower',
         ('Amax', specs[1][1], specs[1][2], 'viridis'),
         ('mu_max', specs[2][1], specs[2][2], 'plasma')),
        ('AUC_upper_mumax_lower',
         ('AUC', specs[0][1], specs[0][2], 'viridis'),
         ('mu_max', specs[2][1], specs[2][2], 'plasma')),
    ]
    for tag, (u_metric, u_lab, u_cb, u_cm), (l_metric, l_lab, l_cb, l_cm) in combined_specs:
        out_png = FIG_DIR / f'response_surface_grid_{tag}_{WAVELENGTH}.png'
        out_pdf = FIG_DIR / f'response_surface_grid_{tag}_{WAVELENGTH}.pdf'
        r = plot_combined_triangles_grid(
            per_cond,
            upper_metric=u_metric, upper_label=u_lab, upper_cbar=u_cb,
            lower_metric=l_metric, lower_label=l_lab, lower_cbar=l_cb,
            out_png=out_png, out_pdf=out_pdf,
            upper_cmap=u_cm, lower_cmap=l_cm,
        )
        print(f'  upper={u_metric} [{r["upper"][0]:.3f}, {r["upper"][1]:.3f}]  '
              f'lower={l_metric} [{r["lower"][0]:.3f}, {r["lower"][1]:.3f}]')
        print(f'           -> {out_png.name} + .pdf')


if __name__ == '__main__':
    main()
