"""
SN1 Biolog PM08 — condition-level analysis.

Joins per-well kinetics (outputs/growth_long.csv from plot_growth.py) with
well_assignment.xlsx, computes per-replicate kinetic features, aggregates to
(Condition_ID, wavelength) mean ± std, runs stats, and plots a ranked
"best media condition" figure.

Plate mapping
-------------
well_assignment.Plate == 1  ->  raw plate SN1_1c
well_assignment.Plate == 2  ->  raw plate SN1_2c
Plates 3 / 4 in the assignment have no raw-data counterpart (skipped).

Outputs (./outputs/conditions/)
-------------------------------
replicate_features.csv     one row per (plate, well, condition, replicate, wavelength)
condition_summary.csv      mean +/- std of every metric per (Condition_ID, wavelength)
rank_AUC740.png            primary "best condition" ranked bar with error bars
rank_max_delta_740.png     same ranking but by max delta abs740 (max OD-equivalent)
top10_kinetics.png         mean +/- std curves for top 10 conditions (740 nm)
plate_concordance.png      scatter of mean AUC740 on plate 1 vs plate 2 (5 overlapping conds)
heatmap_metrics.png        z-scored heatmap of all metrics across conditions

Run:
    python analyze_conditions.py
(Re-run plot_growth.py first if outputs/growth_long.csv is missing.)
"""

from __future__ import annotations

import colorsys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
from scipy.signal import savgol_filter
import plotly.graph_objects as go


HERE = Path(__file__).resolve().parent
GROWTH_CSV = HERE / "outputs" / "growth_long.csv"
ASSIGN_XLSX = HERE / "well_assignment.xlsx"
OUT_DIR = HERE / "outputs" / "conditions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLATE_MAP = {1: "SN1_1c", 2: "SN1_2c", 3: "SN1_3c", 4: "SN1_4c"}  # assignment Plate -> raw plate name
CONTROL_CONDITION = "ctrl_media"         # highlighted in per-plate kinetics plots
PRIMARY_METRIC = "AUC_740"               # the one we rank "best" on
TOP_N = 10                               # how many to show in mean-curve plot
# condition names matching any prefix here are dropped before analysis (calibrants, etc.)
EXCLUDE_PREFIXES = ("Nd",)


# ---------- load + join ----------

def normalize_well(w: str) -> str:
    """A8 -> A08, H11 -> H11."""
    w = str(w).strip()
    return f"{w[0].upper()}{int(w[1:]):02d}"


def load_assignment() -> pd.DataFrame:
    df = pd.read_excel(ASSIGN_XLSX)
    df["Well"] = df["Well"].map(normalize_well)
    df["plate"] = df["Plate"].map(PLATE_MAP)
    df = df.dropna(subset=["plate"]).copy()           # drop plates 3 / 4
    df = df.rename(columns={"Well": "well", "Condition_ID": "condition", "Replicate": "rep"})
    if EXCLUDE_PREFIXES:
        excl_mask = df["condition"].str.startswith(EXCLUDE_PREFIXES)
        if excl_mask.any():
            dropped = sorted(df.loc[excl_mask, "condition"].unique())
            print(f"excluding {len(dropped)} calibrant condition(s): {', '.join(dropped)}")
            df = df[~excl_mask].copy()
    return df[["plate", "well", "condition", "rep"]]


def load_growth() -> pd.DataFrame:
    if not GROWTH_CSV.exists():
        raise FileNotFoundError(f"{GROWTH_CSV} missing — run plot_growth.py first")
    return pd.read_csv(GROWTH_CSV)


# ---------- per-replicate kinetic features ----------

def calculate_mu_max(t: np.ndarray, y: np.ndarray,
                     win_h: float = 4.0,
                     min_r2: float = 0.95,
                     min_delta_for_fit: float = 0.015,
                     min_window_points: int = 6,
                     smooth_window: int = 7,
                     smooth_polyorder: int = 2,
                     epsilon: float = 1e-4) -> dict:
    """Specific growth rate µ_max from an absorbance trace, following standard practice.

    Recipe
    ------
    1. Sort by t, subtract per-well t=0 baseline (cell-specific signal Δabs).
    2. Savitzky-Golay smooth Δabs to suppress measurement noise.
    3. Sliding linear regression on ln(Δabs) over a `win_h`-hour window.
       (Subtracting the blank before log avoids underestimating µ when the
       background absorbance dominates raw OD.)
    4. Reject windows that don't look exponential: require
         (a) ALL points in the window have Δabs >= `min_delta_for_fit`
             (so we're well above the noise floor on every fitted point —
              this is what restricts the fit to genuine exponential phase
              rather than the noise→detection transition where log is unstable),
         (b) at least `min_window_points` points in the window
             (avoids spurious near-perfect R² from 3-point fits),
         (c) R² >= `min_r2` (truly log-linear, not e.g. plateauing).
    5. µ_max = max valid slope across all accepted windows.

    Returns: {mu_max, mu_max_t (window centre), mu_max_r2,
              mu_max_window_n (points in best window)}.  All NaN if no window passes.
    """
    order = np.argsort(t)
    t = np.asarray(t, dtype=float)[order]
    y = np.asarray(y, dtype=float)[order]
    if len(t) < max(smooth_window, 4):
        return {"mu_max": np.nan, "mu_max_t": np.nan,
                "mu_max_r2": np.nan, "mu_max_window_n": 0}

    dy = y - y[0]                                            # blank-subtracted (Δabs)
    win = min(smooth_window, len(dy) - (1 - len(dy) % 2))    # ensure odd & ≤ N
    if win < smooth_polyorder + 2:
        dy_smooth = dy
    else:
        if win % 2 == 0:
            win -= 1
        dy_smooth = savgol_filter(dy, win, smooth_polyorder)

    if dy_smooth.max() < min_delta_for_fit:                  # never grew
        return {"mu_max": np.nan, "mu_max_t": np.nan,
                "mu_max_r2": np.nan, "mu_max_window_n": 0}

    log_y = np.log(np.clip(dy_smooth, epsilon, None))

    best = {"slope": -np.inf, "t_centre": np.nan, "r2": np.nan, "n": 0}
    for i in range(len(t)):
        mask = (t >= t[i]) & (t <= t[i] + win_h)
        if mask.sum() < min_window_points:
            continue
        tt = t[mask]
        yy = log_y[mask]
        # require EVERY point in the window to be above the noise floor
        # (excludes the noise→detection transition, where log is dominated by ε)
        if dy_smooth[mask].min() < min_delta_for_fit:
            continue
        slope, intercept = np.polyfit(tt, yy, 1)
        y_fit = slope * tt + intercept
        ss_res = float(((yy - y_fit) ** 2).sum())
        ss_tot = float(((yy - yy.mean()) ** 2).sum())
        if ss_tot < 1e-12:
            continue
        r2 = 1.0 - ss_res / ss_tot
        if r2 < min_r2 or slope <= 0:
            continue
        if slope > best["slope"]:
            best = {"slope": float(slope),
                    "t_centre": float(t[i] + win_h / 2),
                    "r2": float(r2),
                    "n": int(mask.sum())}

    if not np.isfinite(best["slope"]):
        return {"mu_max": np.nan, "mu_max_t": np.nan,
                "mu_max_r2": np.nan, "mu_max_window_n": 0}
    return {"mu_max": best["slope"], "mu_max_t": best["t_centre"],
            "mu_max_r2": best["r2"], "mu_max_window_n": best["n"]}


def features_one_well(t: np.ndarray, y: np.ndarray) -> dict:
    """Compute summary metrics for a single (well, wavelength) curve."""
    order = np.argsort(t)
    t = t[order]
    y = y[order]
    y0 = y[0]
    dy = y - y0                                       # baseline-subtracted (delta abs)

    auc = float(np.trapezoid(dy, t))                 # trapezoid AUC of delta abs over hours
    max_d = float(dy.max())
    t_max = float(t[int(np.argmax(dy))])
    final = float(dy[-1])

    # specific growth rate via the best-practice routine above
    mu = calculate_mu_max(t, y)
    mu_max = mu["mu_max"]

    # lag = first time delta crosses 10% of max (only meaningful if there's growth)
    lag = np.nan
    if max_d > 1e-3:
        thr = 0.1 * max_d
        crossed = np.where(dy > thr)[0]
        if len(crossed):
            lag = float(t[crossed[0]])

    # doubling time = ln(2)/µ. Only meaningful when µ_max passed the R² gate.
    doubling_time = np.nan
    if np.isfinite(mu_max) and mu_max > 1e-3:
        doubling_time = float(np.log(2) / mu_max)

    return {"AUC": auc, "max_delta": max_d, "t_max": t_max, "final_delta": final,
            "mu_max": mu_max,
            "mu_max_t": mu["mu_max_t"],
            "mu_max_r2": mu["mu_max_r2"],
            "doubling_time": doubling_time,
            "lag": lag, "y0": float(y0)}


def replicate_features(growth: pd.DataFrame, assign: pd.DataFrame) -> pd.DataFrame:
    rows = []
    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    for (plate, well, cond, rep, wl), g in joined.groupby(
        ["plate", "well", "condition", "rep", "wavelength"]
    ):
        f = features_one_well(g["t_h"].to_numpy(), g["abs"].to_numpy())
        rows.append({"plate": plate, "well": well, "condition": cond, "rep": rep,
                     "wavelength": wl, **f})
    return pd.DataFrame(rows)


# ---------- condition-level aggregate ----------

def summarize(rep_df: pd.DataFrame) -> pd.DataFrame:
    """One row per (condition, wavelength) with mean / std / n / sem of every metric."""
    metrics = ["AUC", "max_delta", "mu_max", "doubling_time",
               "lag", "t_max", "final_delta", "y0"]
    agg = (rep_df
           .groupby(["condition", "wavelength"])[metrics]
           .agg(["mean", "std", "count"])
           .reset_index())
    # flatten MultiIndex columns
    agg.columns = ["_".join([str(x) for x in c if x]).rstrip("_") for c in agg.columns]
    # add sem
    for m in metrics:
        agg[f"{m}_sem"] = agg[f"{m}_std"] / np.sqrt(agg[f"{m}_count"].clip(lower=1))
    return agg


def clean_summary(wide: pd.DataFrame) -> pd.DataFrame:
    """Readable condition-level table: condition + (mean, std) for each headline metric,
    columns ordered logically, rows sorted by µ_max @ 740 nm desc, values rounded."""
    out_cols = [
        ("condition",                "condition"),
        ("n_reps",                   "AUC_count_740"),
        ("AUC_740_mean",             "AUC_mean_740"),
        ("AUC_740_std",              "AUC_std_740"),
        ("max_delta_740_mean",       "max_delta_mean_740"),
        ("max_delta_740_std",        "max_delta_std_740"),
        ("mu_max_740_mean_per_h",    "mu_max_mean_740"),
        ("mu_max_740_std_per_h",     "mu_max_std_740"),
        ("doubling_time_740_mean_h", "doubling_time_mean_740"),
        ("doubling_time_740_std_h",  "doubling_time_std_740"),
        ("AUC_590_mean",             "AUC_mean_590"),
        ("AUC_590_std",              "AUC_std_590"),
        ("max_delta_590_mean",       "max_delta_mean_590"),
        ("max_delta_590_std",        "max_delta_std_590"),
        ("mu_max_590_mean_per_h",    "mu_max_mean_590"),
        ("mu_max_590_std_per_h",     "mu_max_std_590"),
        ("doubling_time_590_mean_h", "doubling_time_mean_590"),
        ("doubling_time_590_std_h",  "doubling_time_std_590"),
    ]
    df = wide.copy()
    out = pd.DataFrame({new: df[old] for new, old in out_cols if old in df.columns})
    out = out.sort_values("mu_max_740_mean_per_h", ascending=False, na_position="last")
    # round numeric columns: integer for n_reps, 3 dp for everything else
    for c in out.columns:
        if c == "n_reps":
            out[c] = out[c].astype("Int64")
        elif pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].round(3)
    return out


def wide_by_wavelength(summary: pd.DataFrame) -> pd.DataFrame:
    """Pivot so each condition has both 740 and 590 columns (e.g. AUC_740, AUC_590)."""
    metrics = ["AUC_mean", "AUC_std", "AUC_sem",
               "max_delta_mean", "max_delta_std",
               "mu_max_mean", "mu_max_std",
               "doubling_time_mean", "doubling_time_std", "doubling_time_count",
               "lag_mean", "lag_std",
               "AUC_count"]
    out = summary.pivot(index="condition", columns="wavelength", values=metrics)
    out.columns = [f"{m}_{int(wl)}" for m, wl in out.columns]
    return out.reset_index()


# ---------- stats ----------

def anova_across_conditions(rep_df: pd.DataFrame, metric: str, wavelength: int) -> tuple[float, float]:
    sub = rep_df[rep_df["wavelength"] == wavelength].dropna(subset=[metric])
    groups = [g[metric].to_numpy() for _, g in sub.groupby("condition") if len(g) >= 2]
    if len(groups) < 2:
        return (np.nan, np.nan)
    F, p = stats.f_oneway(*groups)
    return float(F), float(p)


# ---------- plots ----------

def plot_ranked(wide: pd.DataFrame, mean_col: str, std_col: str, title: str, out: Path) -> None:
    df = wide.dropna(subset=[mean_col]).sort_values(mean_col, ascending=True).copy()
    n = len(df)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.22 * n)))
    y = np.arange(n)
    means = df[mean_col].to_numpy()
    stds = df[std_col].fillna(0).to_numpy()
    # color top 10% green, bottom 10% grey
    q90 = np.quantile(means, 0.90)
    q10 = np.quantile(means, 0.10)
    colors = np.where(means >= q90, "tab:green",
              np.where(means <= q10, "lightgrey", "tab:blue"))
    ax.barh(y, means, xerr=stds, color=colors, edgecolor="black", lw=0.4,
            error_kw=dict(lw=0.7, capsize=2))
    ax.set_yticks(y)
    ax.set_yticklabels(df["condition"].to_list(), fontsize=8)
    ax.set_xlabel(mean_col)
    ax.set_title(title)
    ax.axvline(0, color="black", lw=0.5)
    ax.grid(axis="x", lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_top_kinetics(growth: pd.DataFrame, assign: pd.DataFrame, wide: pd.DataFrame,
                       wavelength: int, n: int, out: Path) -> None:
    """For top-N conditions by AUC at this wavelength, plot mean+/-std of delta abs over time."""
    rank_col = f"AUC_mean_{wavelength}"
    top = wide.dropna(subset=[rank_col]).nlargest(n, rank_col)["condition"].to_list()

    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    joined = joined[joined["wavelength"] == wavelength]
    # baseline-subtract per well (delta vs that well's t=0)
    joined["abs_delta"] = joined["abs"] - joined.groupby(["plate", "well", "wavelength"])["abs"].transform("first")

    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.colormaps["tab10"]
    for i, cond in enumerate(top):
        sub = joined[joined["condition"] == cond]
        # bin time per plate (timestamps shared within plate); mean across replicates per t_h
        agg = sub.groupby(["t_h"])["abs_delta"].agg(["mean", "std"]).reset_index()
        agg = agg.sort_values("t_h")
        ax.plot(agg["t_h"], agg["mean"], color=cmap(i % 10), lw=1.5, label=cond)
        ax.fill_between(agg["t_h"],
                        agg["mean"] - agg["std"].fillna(0),
                        agg["mean"] + agg["std"].fillna(0),
                        color=cmap(i % 10), alpha=0.15, lw=0)
    ax.set_xlabel("time (h)")
    ax.set_ylabel(f"delta abs{wavelength}  (mean +/- std across reps)")
    ax.set_title(f"Top {n} conditions by mean AUC_{wavelength}")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_top_n_vs_control_report(growth: pd.DataFrame, assign: pd.DataFrame,
                                  wide: pd.DataFrame, n: int = 5,
                                  wavelength: int = 740,
                                  out_png: Path = None, out_pdf: Path = None) -> list[str]:
    """Publication-style figure: top N conditions by AUC + ctrl_media overlaid,
    mean ± shaded std across all replicate wells (pooled across plates).
    Returns the ordered list of condition names plotted."""
    rank_col = f"AUC_mean_{wavelength}"
    top = (wide.dropna(subset=[rank_col])
                .sort_values(rank_col, ascending=False)
                .loc[lambda d: d["condition"] != CONTROL_CONDITION]
                .head(n)["condition"].tolist())
    selected = top + [CONTROL_CONDITION]

    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    sub = joined[(joined["condition"].isin(selected))
                 & (joined["wavelength"] == wavelength)].copy()
    sub["abs_delta"] = (
        sub["abs"] - sub.groupby(["plate", "well", "wavelength"])["abs"].transform("first")
    )
    sub["t_bin"] = (sub["t_h"] * 10).round() / 10
    agg = (sub.groupby(["condition", "t_bin"])["abs_delta"]
              .agg(["mean", "std", "count"]).reset_index())
    agg["std"] = agg["std"].fillna(0.0)

    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.colormaps["viridis"](np.linspace(0.08, 0.85, n))

    # top conditions
    for i, cond in enumerate(top):
        df = agg[agg["condition"] == cond].sort_values("t_bin")
        nrep = int(df["count"].max())
        ax.plot(df["t_bin"], df["mean"], color=cmap[i], lw=2.2,
                label=f"{cond}  (n={nrep})", zorder=5 - 0.01 * i)
        ax.fill_between(df["t_bin"],
                        df["mean"] - df["std"], df["mean"] + df["std"],
                        color=cmap[i], alpha=0.18, lw=0)

    # control on top, dashed black so it's instantly readable
    ctrl = agg[agg["condition"] == CONTROL_CONDITION].sort_values("t_bin")
    nrep_c = int(ctrl["count"].max())
    ax.plot(ctrl["t_bin"], ctrl["mean"], color="black", lw=2.4, ls="--",
            label=f"{CONTROL_CONDITION}  (control, n={nrep_c})", zorder=10)
    ax.fill_between(ctrl["t_bin"],
                    ctrl["mean"] - ctrl["std"], ctrl["mean"] + ctrl["std"],
                    color="black", alpha=0.12, lw=0, zorder=9)

    ax.set_xlabel("Time (h)", fontsize=12)
    ax.set_ylabel(f"ΔAbs$_{{{wavelength}\\,\\mathrm{{nm}}}}$  "
                  f"(baseline-subtracted, mean ± std)", fontsize=12)
    ax.set_title(f"Top {n} media conditions vs ctrl_media — Biolog PM08 @ {wavelength} nm",
                 fontsize=13)
    ax.legend(loc="upper left", fontsize=10, frameon=False, handlelength=2.5)
    ax.grid(lw=0.3, alpha=0.35)
    ax.axhline(0, color="black", lw=0.4, alpha=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=10)
    fig.tight_layout()

    if out_png:
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
    if out_pdf:
        fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return top


def improvement_vs_control_table(rep_df: pd.DataFrame, top_conds: list[str],
                                  wavelength: int = 740,
                                  out_csv: Path = None) -> pd.DataFrame:
    """For each top condition, compare to PLATE-MATCHED ctrl_media (not pooled).
    This is the rigorous fold-change because per-plate controls already absorb
    plate-effect (e.g. Plate 3 controls ran ~28% lower than Plate 1 controls)."""
    sub = rep_df[rep_df["wavelength"] == wavelength]

    # per-plate control mean of AUC and µ_max
    ctrl_by_plate = (
        sub[sub["condition"] == CONTROL_CONDITION]
        .groupby("plate").agg(ctrl_AUC=("AUC", "mean"),
                              ctrl_AUC_std=("AUC", "std"),
                              ctrl_mu=("mu_max", "mean"),
                              ctrl_mu_std=("mu_max", "std"))
    )

    rows = []
    for cond in top_conds:
        cond_rows = sub[sub["condition"] == cond]
        if cond_rows.empty:
            continue
        plates_for_cond = sorted(cond_rows["plate"].unique())
        # control reference = mean across the plates this condition is on
        ctrl_AUC_ref = ctrl_by_plate.loc[plates_for_cond, "ctrl_AUC"].mean()
        ctrl_mu_ref  = ctrl_by_plate.loc[plates_for_cond, "ctrl_mu"].mean()

        auc_mean = cond_rows["AUC"].mean()
        auc_std  = cond_rows["AUC"].std(ddof=1)
        mu_mean  = cond_rows["mu_max"].mean()
        mu_std   = cond_rows["mu_max"].std(ddof=1)

        rows.append({
            "condition": cond,
            "n_reps": int(len(cond_rows)),
            "plates_present": ", ".join(plates_for_cond),
            f"AUC_{wavelength}_mean":             round(auc_mean, 3),
            f"AUC_{wavelength}_std":              round(auc_std, 3),
            f"ctrl_AUC_{wavelength}_matched":     round(ctrl_AUC_ref, 3),
            f"AUC_{wavelength}_pct_vs_ctrl":      round((auc_mean / ctrl_AUC_ref - 1) * 100, 1),
            f"mu_max_{wavelength}_mean_per_h":    round(mu_mean, 3),
            f"mu_max_{wavelength}_std_per_h":     round(mu_std, 3),
            f"ctrl_mu_max_{wavelength}_matched":  round(ctrl_mu_ref, 3),
            f"mu_max_{wavelength}_pct_vs_ctrl":   round((mu_mean / ctrl_mu_ref - 1) * 100, 1),
        })

    out = pd.DataFrame(rows)
    if out_csv:
        out.to_csv(out_csv, index=False)
    return out


def plot_control_across_plates(growth: pd.DataFrame, assign: pd.DataFrame,
                                out: Path) -> None:
    """Overlay ctrl_media (mean +/- std across reps) from each plate, both wavelengths."""
    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    ctrl = joined[joined["condition"] == CONTROL_CONDITION].copy()
    if ctrl.empty:
        print(f"  [skip] no {CONTROL_CONDITION} wells found")
        return
    ctrl["abs_delta"] = (
        ctrl["abs"] - ctrl.groupby(["plate", "well", "wavelength"])["abs"].transform("first")
    )
    ctrl["t_bin"] = (ctrl["t_h"] * 10).round() / 10

    plates = sorted(ctrl["plate"].unique())
    palette = _hsv_palette(max(4, len(plates)))
    plate_color = {p: palette[i] for i, p in enumerate(plates)}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharex=True)
    summary_rows = []

    for ax, wl in zip(axes, (740, 590)):
        for plate in plates:
            sub = ctrl[(ctrl["plate"] == plate) & (ctrl["wavelength"] == wl)]
            if sub.empty:
                continue
            agg = (sub.groupby("t_bin")["abs_delta"]
                      .agg(["mean", "std", "count"]).reset_index().sort_values("t_bin"))
            agg["std"] = agg["std"].fillna(0.0)
            n = int(agg["count"].max())
            color = plate_color[plate]
            ax.plot(agg["t_bin"], agg["mean"], color=color, lw=2.0,
                    label=f"{plate} (n={n})")
            ax.fill_between(agg["t_bin"],
                            agg["mean"] - agg["std"],
                            agg["mean"] + agg["std"],
                            color=color, alpha=0.18, lw=0)
            # per-plate AUC and final delta for the console summary
            mean_curve = agg.set_index("t_bin")["mean"]
            auc = float(np.trapezoid(mean_curve.values, mean_curve.index.values))
            summary_rows.append({"plate": plate, "wavelength": wl,
                                  "n_reps": n,
                                  "AUC": round(auc, 3),
                                  "final_delta": round(float(mean_curve.iloc[-1]), 3),
                                  "max_delta": round(float(mean_curve.max()), 3)})

        ax.axhline(0, color="black", lw=0.4, alpha=0.4)
        ax.set_xlabel("time (h)")
        ax.set_ylabel(f"Δabs{wl}  (mean ± std across reps)")
        ax.set_title(f"ctrl_media @ {wl} nm")
        ax.legend(loc="best", fontsize=9, frameon=False)
        ax.grid(lw=0.3, alpha=0.4)

    fig.suptitle("Plate-to-plate variation of ctrl_media (n=4 reps per plate)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=160)
    plt.close(fig)

    print("\nctrl_media per-plate summary:")
    print(pd.DataFrame(summary_rows).to_string(index=False))


def plot_per_plate_kinetics(growth: pd.DataFrame, assign: pd.DataFrame,
                             plate: str, wavelength: int, out: Path) -> None:
    """Mean +/- std curves for every condition on one plate, with the control highlighted."""
    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    sub = joined[(joined["plate"] == plate) & (joined["wavelength"] == wavelength)].copy()
    if sub.empty:
        print(f"  [skip] no data for plate {plate} at {wavelength} nm")
        return

    sub["abs_delta"] = (
        sub["abs"] - sub.groupby(["plate", "well", "wavelength"])["abs"].transform("first")
    )
    sub["t_bin"] = (sub["t_h"] * 10).round() / 10

    # rank conditions on THIS plate by AUC of the mean curve, descending → consistent legend
    auc_per_cond: dict[str, float] = {}
    for cond, g in sub.groupby("condition"):
        mean_curve = g.groupby("t_bin")["abs_delta"].mean().sort_index()
        auc_per_cond[cond] = float(np.trapezoid(mean_curve.values, mean_curve.index.values))
    cond_order = sorted(auc_per_cond, key=auc_per_cond.get, reverse=True)
    # put control last so it draws on top (highlighted)
    if CONTROL_CONDITION in cond_order:
        cond_order = [c for c in cond_order if c != CONTROL_CONDITION] + [CONTROL_CONDITION]
    palette = _hsv_palette(max(1, len(cond_order) - (1 if CONTROL_CONDITION in cond_order else 0)))

    fig, ax = plt.subplots(figsize=(11, 7))
    p_idx = 0
    for cond in cond_order:
        agg = (sub[sub["condition"] == cond]
               .groupby("t_bin")["abs_delta"].agg(["mean", "std"])
               .reset_index().sort_values("t_bin"))
        agg["std"] = agg["std"].fillna(0.0)
        if cond == CONTROL_CONDITION:
            ax.plot(agg["t_bin"], agg["mean"], color="black", lw=2.5,
                    label=f"{cond} (control)", zorder=10)
            ax.fill_between(agg["t_bin"],
                            agg["mean"] - agg["std"],
                            agg["mean"] + agg["std"],
                            color="black", alpha=0.18, lw=0, zorder=9)
        else:
            color = palette[p_idx]; p_idx += 1
            ax.plot(agg["t_bin"], agg["mean"], color=color, lw=1.0, label=cond)
            ax.fill_between(agg["t_bin"],
                            agg["mean"] - agg["std"],
                            agg["mean"] + agg["std"],
                            color=color, alpha=0.10, lw=0)
    ax.axhline(0, color="black", lw=0.4, alpha=0.4)
    ax.set_xlabel("time (h)")
    ax.set_ylabel(f"Δabs{wavelength}  (mean ± std across reps)")
    n_cond = len(cond_order)
    ax.set_title(f"Plate {plate} — kinetics at {wavelength} nm "
                 f"({n_cond} conditions, Nd calibrants excluded; control = ctrl_media in black)")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=7, ncol=1, frameon=False)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_growth_rate_top_bottom(wide: pd.DataFrame, wavelength: int,
                                 n: int = 10,
                                 control: str = "ctrl_media",
                                 out: Path = None) -> None:
    """Top N + bottom N conditions by µ_max + the control, sorted high→low.
    Two panels (740 & 590) when called per wavelength."""
    mu_col = f"mu_max_mean_{wavelength}"
    mu_std = f"mu_max_std_{wavelength}"
    auc_col = f"AUC_mean_{wavelength}"
    td_col = f"doubling_time_count_{wavelength}"

    df = wide.copy()
    # real growers only (positive µ, positive AUC, ≥2 reps with valid doubling)
    real = (df[mu_col] > 1e-3) & (df[auc_col] > 0) & (df[td_col].fillna(0) >= 2)
    df_real = df[real].dropna(subset=[mu_col]).sort_values(mu_col, ascending=False)

    top = df_real.head(n)
    bottom = df_real.tail(n)
    ctrl_row = df[df["condition"] == control]                # always include even if low

    # build display order: top (descending) → control → bottom (descending)
    pieces = [top]
    if not ctrl_row.empty and control not in top["condition"].values \
       and control not in bottom["condition"].values:
        pieces.append(ctrl_row)
    pieces.append(bottom)
    plot_df = pd.concat(pieces).drop_duplicates(subset=["condition"]).copy()

    # color/style per row
    def row_style(cond: str) -> tuple[str, str, float]:
        if cond == control:
            return ("#000000", "#FFA500", 1.4)        # face=black, edge=orange
        if cond in set(top["condition"]):
            return ("#2ca02c", "#222222", 0.5)        # green
        return ("#a6a6a6", "#222222", 0.5)            # grey

    n_rows = len(plot_df)
    fig, ax = plt.subplots(figsize=(9, max(5, 0.32 * n_rows)))
    y_pos = np.arange(n_rows)[::-1]                  # invert so highest is at top
    means = plot_df[mu_col].to_numpy()
    stds = plot_df[mu_std].fillna(0).to_numpy()
    faces, edges, ews = zip(*[row_style(c) for c in plot_df["condition"]])

    bars = ax.barh(y_pos, means, xerr=stds,
                   color=list(faces), edgecolor=list(edges),
                   linewidth=list(ews),
                   error_kw=dict(lw=0.7, capsize=2))

    ax.set_yticks(y_pos)
    labels = [f"{c} (control)" if c == control else c for c in plot_df["condition"]]
    ax.set_yticklabels(labels, fontsize=9)

    # bold the control label
    for tick, c in zip(ax.get_yticklabels(), plot_df["condition"]):
        if c == control:
            tick.set_fontweight("bold")
            tick.set_color("#cc6600")

    # visual separator between top group and (control / bottom)
    if not top.empty and not bottom.empty:
        sep_y = y_pos[len(top) - 1] - 0.5
        ax.axhline(sep_y, color="#aaa", lw=0.6, ls="--", alpha=0.7)
        if (control not in top["condition"].values
            and not ctrl_row.empty and control not in bottom["condition"].values):
            sep_y2 = y_pos[len(top)] - 0.5
            ax.axhline(sep_y2, color="#aaa", lw=0.6, ls="--", alpha=0.7)

    ax.set_xlabel(f"specific growth rate µ_max @ {wavelength} nm  [h⁻¹]")
    ax.set_title(
        f"µ_max — top {n} (green) + ctrl_media (black/orange) + bottom {n} (grey)\n"
        f"sorted high → low; error bars = std across reps"
    )
    ax.grid(axis="x", lw=0.3, alpha=0.4)
    ax.axvline(0, color="black", lw=0.5)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_growth_rate(wide: pd.DataFrame, wavelength: int, out: Path) -> None:
    """Two-panel ranked bar: specific growth rate µ_max and doubling time t_d.
    Conditions with no real growth (µ < 1e-3 or NaN doubling time) are dropped.
    """
    mu_col, mu_std = f"mu_max_mean_{wavelength}", f"mu_max_std_{wavelength}"
    td_col, td_std = f"doubling_time_mean_{wavelength}", f"doubling_time_std_{wavelength}"

    df = wide.dropna(subset=[mu_col]).copy()
    auc_col = f"AUC_mean_{wavelength}"
    td_count = f"doubling_time_count_{wavelength}"
    # require: positive µ, net-positive AUC (real growth, not decline), and a doubling
    # time supported by ≥ 2 replicates (drops noisy single-rep outliers like Nd wells)
    grew = (df[mu_col] > 1e-3) & (df[auc_col] > 0) & (df[td_count].fillna(0) >= 2)
    skipped = sorted(df.loc[~grew, "condition"].tolist())
    df = df[grew].copy()

    df_mu = df.sort_values(mu_col, ascending=True)
    df_td = df.dropna(subset=[td_col]).sort_values(td_col, ascending=False)  # fastest at top

    fig, axes = plt.subplots(1, 2, figsize=(13, max(5, 0.22 * len(df_mu))))

    # left panel: µ_max
    ax = axes[0]
    y = np.arange(len(df_mu))
    means = df_mu[mu_col].to_numpy()
    stds = df_mu[mu_std].fillna(0).to_numpy()
    q90 = np.quantile(means, 0.90)
    colors = np.where(means >= q90, "tab:green", "tab:blue")
    ax.barh(y, means, xerr=stds, color=colors, edgecolor="black", lw=0.4,
            error_kw=dict(lw=0.7, capsize=2))
    ax.set_yticks(y)
    ax.set_yticklabels(df_mu["condition"].to_list(), fontsize=8)
    ax.set_xlabel(f"specific growth rate µ_max ({wavelength} nm)  [h⁻¹]")
    ax.set_title("µ_max — error bars = std across reps")
    ax.grid(axis="x", lw=0.3, alpha=0.4)

    # right panel: doubling time (fastest at the top, i.e. shortest time)
    ax = axes[1]
    y = np.arange(len(df_td))
    means = df_td[td_col].to_numpy()
    stds = df_td[td_std].fillna(0).to_numpy()
    q10 = np.quantile(means, 0.10)
    colors = np.where(means <= q10, "tab:green", "tab:blue")
    ax.barh(y, means, xerr=stds, color=colors, edgecolor="black", lw=0.4,
            error_kw=dict(lw=0.7, capsize=2))
    ax.set_yticks(y)
    ax.set_yticklabels(df_td["condition"].to_list(), fontsize=8)
    ax.set_xlabel(f"doubling time t_d ({wavelength} nm)  [h]")
    ax.set_title("t_d = ln(2)/µ — error bars = std across reps")
    ax.grid(axis="x", lw=0.3, alpha=0.4)

    suptitle = f"Specific growth rate and doubling time at {wavelength} nm"
    if skipped:
        suptitle += f"\n({len(skipped)} non-growers excluded: {', '.join(skipped)})"
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_all_conditions_kinetics(growth: pd.DataFrame, assign: pd.DataFrame,
                                  wide: pd.DataFrame, wavelength: int,
                                  out: Path) -> None:
    """Static PNG: mean +/- std band per condition for all conditions, sorted by AUC_740."""
    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    joined = joined[joined["wavelength"] == wavelength].copy()
    joined["abs_delta"] = (
        joined["abs"] - joined.groupby(["plate", "well", "wavelength"])["abs"].transform("first")
    )
    joined["t_bin"] = (joined["t_h"] * 10).round() / 10

    rank_col = f"AUC_mean_{wavelength}"
    rank = (wide.dropna(subset=[rank_col])
                .sort_values(rank_col, ascending=False)["condition"].tolist())
    palette = _hsv_palette(len(rank))

    fig, ax = plt.subplots(figsize=(11, 7))
    for i, cond in enumerate(rank):
        sub = joined[joined["condition"] == cond]
        agg = (sub.groupby("t_bin")["abs_delta"]
                  .agg(["mean", "std"]).reset_index().sort_values("t_bin"))
        agg["std"] = agg["std"].fillna(0.0)
        ax.plot(agg["t_bin"], agg["mean"], color=palette[i], lw=1.0, label=cond)
        ax.fill_between(agg["t_bin"],
                        agg["mean"] - agg["std"],
                        agg["mean"] + agg["std"],
                        color=palette[i], alpha=0.10, lw=0)
    ax.axhline(0, color="black", lw=0.4, alpha=0.5)
    ax.set_xlabel("time (h)")
    ax.set_ylabel(f"Δabs{wavelength}  (mean ± std across reps)")
    ax.set_title(f"All {len(rank)} conditions — mean ± std at {wavelength} nm "
                 f"(sorted by AUC_{wavelength})")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=7, ncol=2, frameon=False)
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_plate_concordance(rep_df: pd.DataFrame, out: Path) -> None:
    """For conditions present on both plates, scatter mean AUC740 plate1 vs plate2."""
    sub = rep_df[rep_df["wavelength"] == 740]
    means = sub.groupby(["plate", "condition"])["AUC"].mean().unstack("plate")
    if means.shape[1] < 2:
        return
    means = means.dropna()
    if means.empty:
        return
    p1, p2 = means.columns[0], means.columns[1]
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(means[p1], means[p2], s=30, alpha=0.8, edgecolor="black", lw=0.5)
    for cond, row in means.iterrows():
        ax.annotate(cond, (row[p1], row[p2]), fontsize=7, alpha=0.7)
    lo = min(means[p1].min(), means[p2].min())
    hi = max(means[p1].max(), means[p2].max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.6, alpha=0.5)
    if len(means) >= 3:
        r, p = stats.pearsonr(means[p1], means[p2])
        ax.set_title(f"AUC_740 mean: {p1} vs {p2}\n"
                     f"n={len(means)} shared conditions, Pearson r={r:.2f} (p={p:.3g})")
    else:
        ax.set_title(f"AUC_740 mean: {p1} vs {p2} (n={len(means)})")
    ax.set_xlabel(f"{p1}  mean AUC_740")
    ax.set_ylabel(f"{p2}  mean AUC_740")
    ax.grid(lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def _hsv_palette(n: int) -> list[str]:
    """n evenly-spaced hex colors via HSV."""
    return [
        "#" + "".join(f"{int(c*255):02x}" for c in colorsys.hsv_to_rgb(i / n, 0.62, 0.85))
        for i in range(n)
    ]


def _hex_to_rgba(hex_str: str, alpha: float) -> str:
    h = hex_str.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _condition_group(name: str) -> str:
    if name.startswith("MPOB"):
        return "MPOB"
    if name.lower().startswith("nd"):
        return "Nd series"
    if name.lower().startswith("ctrl"):
        return "Control"
    return "Other"


def build_condition_viewer(growth: pd.DataFrame, assign: pd.DataFrame,
                           wide: pd.DataFrame, out_path: Path,
                           default_n_visible: int = 10) -> None:
    """Standalone HTML with a real checkbox sidebar driving a Plotly figure."""
    joined = growth.merge(assign, on=["plate", "well"], how="inner")
    joined["abs_delta"] = (
        joined["abs"] - joined.groupby(["plate", "well", "wavelength"])["abs"].transform("first")
    )
    joined["t_bin"] = (joined["t_h"] * 10).round() / 10

    agg = (joined.groupby(["condition", "wavelength", "t_bin"])["abs_delta"]
           .agg(["mean", "std", "count"])
           .reset_index()
           .sort_values(["condition", "wavelength", "t_bin"]))
    agg["std"] = agg["std"].fillna(0.0)

    rank = (wide.dropna(subset=["AUC_mean_740"])
                .sort_values("AUC_mean_740", ascending=False)["condition"].tolist())
    extras = [c for c in agg["condition"].unique() if c not in rank]
    ordered = rank + extras
    visible_default = set(rank[:default_n_visible])
    palette = _hsv_palette(len(ordered))
    cond_color = {c: palette[i] for i, c in enumerate(ordered)}

    fig = go.Figure()
    trace_index: dict[str, dict[str, dict[str, int]]] = {}  # {cond: {wl: {role: idx}}}
    cur_idx = 0

    for cond in ordered:
        color = cond_color[cond]
        for wl in (740, 590):
            sub = agg[(agg["condition"] == cond) & (agg["wavelength"] == wl)]
            if sub.empty:
                continue
            t = sub["t_bin"].to_numpy()
            mean = sub["mean"].to_numpy()
            std = sub["std"].to_numpy()
            on = (cond in visible_default) and (wl == 740)

            fig.add_trace(go.Scatter(
                x=t, y=mean + std, line=dict(width=0), mode="lines",
                visible=on, hoverinfo="skip", showlegend=False,
                name=f"{cond} {wl} +std",
            ))
            trace_index.setdefault(cond, {}).setdefault(str(wl), {})["upper"] = cur_idx
            cur_idx += 1

            fig.add_trace(go.Scatter(
                x=t, y=mean - std, line=dict(width=0), mode="lines",
                fill="tonexty", fillcolor=_hex_to_rgba(color, 0.18),
                visible=on, hoverinfo="skip", showlegend=False,
                name=f"{cond} {wl} -std",
            ))
            trace_index[cond][str(wl)]["lower"] = cur_idx
            cur_idx += 1

            fig.add_trace(go.Scatter(
                x=t, y=mean, name=f"{cond} {wl}",
                line=dict(color=color, width=1.7,
                          dash=("dash" if wl == 590 else "solid")),
                mode="lines", visible=on, showlegend=False,
                hovertemplate=f"{cond} {wl}nm<br>t=%{{x:.1f}} h<br>"
                              "Δabs=%{y:.3f}<extra></extra>",
            ))
            trace_index[cond][str(wl)]["mean"] = cur_idx
            cur_idx += 1

    fig.update_layout(
        height=720,
        margin=dict(l=70, r=20, t=20, b=60),
        showlegend=False,
        xaxis_title="time (h)",
        yaxis_title="Δabs (baseline-subtracted)",
        plot_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor="#eee"),
        yaxis=dict(showgrid=True, gridcolor="#eee", zeroline=True, zerolinecolor="#999"),
    )
    # include_plotlyjs=True embeds plotly.js (~3.5 MB) so the file works offline /
    # behind firewalls that block CDNs. Switch to "cdn" to get a smaller file
    # if collaborators always have internet access.
    plot_div = fig.to_html(include_plotlyjs=True, full_html=False, div_id="plot",
                           config={"responsive": True, "displaylogo": False})

    # metadata for the sidebar JS
    auc740 = wide.set_index("condition")["AUC_mean_740"].to_dict()
    auc740_std = wide.set_index("condition")["AUC_std_740"].to_dict()
    auc590 = wide.set_index("condition")["AUC_mean_590"].to_dict()
    cond_info = {
        c: {
            "color": cond_color[c],
            "auc740": None if pd.isna(auc740.get(c, np.nan)) else round(float(auc740[c]), 2),
            "auc740_std": None if pd.isna(auc740_std.get(c, np.nan)) else round(float(auc740_std[c]), 2),
            "auc590": None if pd.isna(auc590.get(c, np.nan)) else round(float(auc590.get(c, np.nan)), 2),
            "group": _condition_group(c),
        }
        for c in ordered
    }
    groups: dict[str, list[str]] = {}
    for c in ordered:
        groups.setdefault(_condition_group(c), []).append(c)
    group_order = ["MPOB", "Nd series", "Control", "Other"]
    groups = {g: groups[g] for g in group_order if g in groups}

    top10 = rank[:10]
    bottom10 = rank[-10:]

    payload = {
        "trace_index": trace_index,
        "cond_info": cond_info,
        "groups": groups,
        "top10": top10,
        "bottom10": bottom10,
        "all_conditions": ordered,
        "default_visible": list(visible_default),
    }

    out_path.write_text(_VIEWER_TEMPLATE.replace("__PLOT_DIV__", plot_div)
                                         .replace("__PAYLOAD__", json.dumps(payload)))


_VIEWER_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SN1 Biolog PM08 — condition viewer</title>
<style>
  * { box-sizing: border-box; }
  body { margin: 0; font-family: -apple-system, system-ui, "Helvetica Neue", Arial, sans-serif;
         color: #222; background: #fff; }
  .app { display: flex; height: 100vh; min-height: 600px; }
  .sidebar {
    width: 340px; flex-shrink: 0; border-right: 1px solid #ddd; background: #fafafa;
    display: flex; flex-direction: column; overflow: hidden;
  }
  .sidebar-header { padding: 12px 14px 8px; border-bottom: 1px solid #e6e6e6; }
  .sidebar-header h1 { font-size: 14px; margin: 0 0 6px; font-weight: 600; }
  .sidebar-header .sub { font-size: 11px; color: #666; }
  .sidebar-controls { padding: 10px 14px; border-bottom: 1px solid #e6e6e6; }
  .sidebar-list { flex: 1; overflow-y: auto; padding: 4px 8px 12px; }
  .plot-area { flex: 1; min-width: 0; padding: 12px; overflow: auto; }

  input[type="search"] {
    width: 100%; padding: 6px 8px; font-size: 12px;
    border: 1px solid #ccc; border-radius: 4px; margin-bottom: 8px;
  }
  .toolbar { display: flex; flex-wrap: wrap; gap: 4px; margin-bottom: 6px; }
  .toolbar button {
    font-size: 11px; padding: 4px 8px; border: 1px solid #c0c0c0;
    background: #fff; border-radius: 3px; cursor: pointer;
  }
  .toolbar button:hover { background: #eef; border-color: #88a; }

  .opts { display: flex; flex-wrap: wrap; gap: 10px; font-size: 12px; margin-top: 4px; }
  .opts label { display: inline-flex; align-items: center; gap: 4px; cursor: pointer; }

  .group-header {
    font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em;
    color: #555; margin: 10px 4px 2px; padding-bottom: 2px;
    border-bottom: 1px solid #ddd; display: flex; justify-content: space-between;
  }
  .group-header .gcount { font-weight: 400; color: #999; }
  .cond-row {
    display: flex; align-items: center; padding: 3px 4px; font-size: 12px;
    gap: 6px; border-radius: 3px; cursor: pointer;
  }
  .cond-row:hover { background: #eef3ff; }
  .cond-row .cb { margin: 0; cursor: pointer; }
  .swatch { width: 11px; height: 11px; border-radius: 2px; flex-shrink: 0;
            border: 1px solid rgba(0,0,0,0.15); }
  .cond-name { flex: 1; font-variant-numeric: tabular-nums; }
  .cond-auc { font-size: 11px; color: #666; font-variant-numeric: tabular-nums; }
  .hidden { display: none !important; }
  #plot { width: 100%; height: 100%; }
</style>
</head>
<body>
<div class="app">
  <aside class="sidebar">
    <div class="sidebar-header">
      <h1>SN1 Biolog PM08 — condition picker</h1>
      <div class="sub">mean &plusmn; std across replicates &middot; sorted by AUC<sub>740</sub></div>
    </div>

    <div class="sidebar-controls">
      <input type="search" id="search" placeholder="filter conditions… (e.g. MPOB_008, Nd)">
      <div class="toolbar">
        <button data-action="all">All</button>
        <button data-action="none">None</button>
        <button data-action="invert">Invert</button>
        <button data-action="top10">Top 10</button>
        <button data-action="bottom10">Bottom 10</button>
        <button data-action="reset">Reset</button>
      </div>
      <div class="opts">
        <label><input type="checkbox" id="wl740" checked> 740 nm</label>
        <label><input type="checkbox" id="wl590"> 590 nm</label>
        <label><input type="checkbox" id="bands" checked> std band</label>
        <label><input type="checkbox" id="clipy"> clip Y</label>
      </div>
    </div>

    <div class="sidebar-list" id="list"></div>
  </aside>

  <main class="plot-area">__PLOT_DIV__</main>
</div>

<script>
const PAYLOAD = __PAYLOAD__;

function buildList() {
  const list = document.getElementById('list');
  list.innerHTML = '';
  for (const [group, conds] of Object.entries(PAYLOAD.groups)) {
    const header = document.createElement('div');
    header.className = 'group-header';
    header.innerHTML = `<span>${group}</span><span class="gcount">${conds.length}</span>`;
    list.appendChild(header);

    for (const cond of conds) {
      const info = PAYLOAD.cond_info[cond];
      const row = document.createElement('label');
      row.className = 'cond-row';
      row.dataset.name = cond;
      row.dataset.group = group;

      const aucTxt = (info.auc740 === null) ? '—'
        : `${info.auc740.toFixed(2)}${info.auc740_std !== null ? ' ±' + info.auc740_std.toFixed(2) : ''}`;
      const isOn = PAYLOAD.default_visible.includes(cond);

      row.innerHTML = `
        <input type="checkbox" class="cb" value="${cond}" ${isOn ? 'checked' : ''}>
        <span class="swatch" style="background:${info.color}"></span>
        <span class="cond-name">${cond}</span>
        <span class="cond-auc">${aucTxt}</span>
      `;
      list.appendChild(row);
    }
  }

  list.addEventListener('change', (e) => {
    if (e.target.classList.contains('cb')) applyVisibility();
  });
}

function checkboxes() {
  return Array.from(document.querySelectorAll('#list .cb'));
}
function visibleRows() {
  return Array.from(document.querySelectorAll('#list .cond-row:not(.hidden)'));
}
function setChecked(condSet, mode = 'set') {
  for (const cb of checkboxes()) {
    if (mode === 'set')         cb.checked = condSet.has(cb.value);
    else if (mode === 'union')  cb.checked = cb.checked || condSet.has(cb.value);
    else if (mode === 'invert') cb.checked = !cb.checked;
  }
}

function applyVisibility() {
  const checked = new Set(checkboxes().filter(cb => cb.checked).map(cb => cb.value));
  const wl740 = document.getElementById('wl740').checked;
  const wl590 = document.getElementById('wl590').checked;
  const bands = document.getElementById('bands').checked;

  const idxs = [];
  const vis  = [];
  for (const [cond, byWl] of Object.entries(PAYLOAD.trace_index)) {
    const condOn = checked.has(cond);
    for (const [wl, byRole] of Object.entries(byWl)) {
      const wlOn = (wl === '740' && wl740) || (wl === '590' && wl590);
      for (const [role, i] of Object.entries(byRole)) {
        const roleOn = (role === 'mean') ? true : bands;
        idxs.push(i);
        vis.push(condOn && wlOn && roleOn);
      }
    }
  }
  Plotly.restyle('plot', {visible: vis}, idxs);
}

function applyYAxis() {
  const clip = document.getElementById('clipy').checked;
  Plotly.relayout('plot', {'yaxis.autorange': !clip,
                           'yaxis.range': clip ? [-0.05, 0.25] : null});
}

function attachToolbar() {
  document.querySelector('.toolbar').addEventListener('click', (e) => {
    const action = e.target.dataset.action;
    if (!action) return;
    const allConds = new Set(PAYLOAD.all_conditions);
    if (action === 'all')          setChecked(allConds);
    else if (action === 'none')    setChecked(new Set());
    else if (action === 'invert')  setChecked(new Set(), 'invert');
    else if (action === 'top10')   setChecked(new Set(PAYLOAD.top10));
    else if (action === 'bottom10')setChecked(new Set(PAYLOAD.bottom10));
    else if (action === 'reset') {
      setChecked(new Set(PAYLOAD.default_visible));
      document.getElementById('wl740').checked = true;
      document.getElementById('wl590').checked = false;
      document.getElementById('bands').checked = true;
      document.getElementById('clipy').checked = false;
      document.getElementById('search').value = '';
      filterList('');
      applyYAxis();
    }
    applyVisibility();
  });

  for (const id of ['wl740', 'wl590', 'bands']) {
    document.getElementById(id).addEventListener('change', applyVisibility);
  }
  document.getElementById('clipy').addEventListener('change', applyYAxis);

  document.getElementById('search').addEventListener('input', (e) => {
    filterList(e.target.value.trim().toLowerCase());
  });
}

function filterList(q) {
  for (const row of document.querySelectorAll('#list .cond-row')) {
    row.classList.toggle('hidden', q && !row.dataset.name.toLowerCase().includes(q));
  }
  // hide empty group headers
  let lastHeader = null, anyVisibleSinceHeader = false;
  for (const el of document.querySelectorAll('#list > *')) {
    if (el.classList.contains('group-header')) {
      if (lastHeader) lastHeader.classList.toggle('hidden', !anyVisibleSinceHeader);
      lastHeader = el;
      anyVisibleSinceHeader = false;
    } else if (!el.classList.contains('hidden')) {
      anyVisibleSinceHeader = true;
    }
  }
  if (lastHeader) lastHeader.classList.toggle('hidden', !anyVisibleSinceHeader);
}

document.addEventListener('DOMContentLoaded', () => {
  buildList();
  attachToolbar();
});
</script>
</body>
</html>
"""


def plot_metric_heatmap(wide: pd.DataFrame, out: Path) -> None:
    cols = ["AUC_mean_740", "max_delta_mean_740", "mu_max_mean_740",
            "AUC_mean_590", "max_delta_mean_590", "mu_max_mean_590"]
    cols = [c for c in cols if c in wide.columns]
    df = wide.set_index("condition")[cols].copy()
    df = df.sort_values("AUC_mean_740", ascending=False)
    # z-score each column independently (across conditions)
    z = (df - df.mean()) / df.std(ddof=0)
    fig, ax = plt.subplots(figsize=(7, max(4, 0.22 * len(z))))
    norm = TwoSlopeNorm(vcenter=0, vmin=z.values.min(), vmax=z.values.max())
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(z.index, fontsize=7)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=30, ha="right", fontsize=8)
    ax.set_title("Per-metric z-score across conditions (sorted by AUC_740)")
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------- main ----------

def main() -> None:
    growth = load_growth()
    assign = load_assignment()
    print(f"loaded {len(growth):,} growth rows, {len(assign):,} well assignments "
          f"({assign['condition'].nunique()} conditions on plates {sorted(assign['plate'].unique())})")

    rep_df = replicate_features(growth, assign)
    rep_csv = OUT_DIR / "replicate_features.csv"
    rep_df.to_csv(rep_csv, index=False)
    print(f"wrote {rep_csv}  ({len(rep_df):,} rows)")

    summary = summarize(rep_df)
    wide = wide_by_wavelength(summary)
    sum_csv = OUT_DIR / "condition_summary.csv"
    wide.to_csv(sum_csv, index=False)
    print(f"wrote {sum_csv}  ({len(wide)} conditions)")

    clean = clean_summary(wide)
    clean_csv = OUT_DIR / "condition_summary_clean.csv"
    clean.to_csv(clean_csv, index=False)
    print(f"wrote {clean_csv}  (sorted by µ_max @ 740 nm)")

    # one-way ANOVA across conditions (overall test of differential growth)
    F740, p740 = anova_across_conditions(rep_df, "AUC", 740)
    F590, p590 = anova_across_conditions(rep_df, "AUC", 590)
    print(f"ANOVA AUC across conditions:  740nm F={F740:.2f}, p={p740:.3g} | "
          f"590nm F={F590:.2f}, p={p590:.3g}")

    # rank: best condition by mean AUC_740 (with std error bars)
    best = wide.sort_values("AUC_mean_740", ascending=False).head(5)
    print("\ntop 5 by mean AUC_740 (biomass-equivalent):")
    print(best[["condition", "AUC_mean_740", "AUC_std_740", "AUC_count_740",
                "max_delta_mean_740", "mu_max_mean_740"]].to_string(index=False))

    # rank by specific growth rate (drop non-growers and single-rep outliers)
    real_growers = (wide["mu_max_mean_740"] > 1e-3) & \
                   (wide["AUC_mean_740"] > 0) & \
                   (wide["doubling_time_count_740"].fillna(0) >= 2)
    fastest = (wide[real_growers]
                   .sort_values("mu_max_mean_740", ascending=False).head(5))
    print("\ntop 5 by mean µ_max @ 740 nm (specific growth rate, h⁻¹):")
    print(fastest[["condition", "mu_max_mean_740", "mu_max_std_740",
                   "doubling_time_mean_740", "doubling_time_std_740"]]
          .round(3).to_string(index=False))

    # plots
    plot_ranked(wide, "AUC_mean_740", "AUC_std_740",
                "Conditions ranked by mean AUC_740 (biomass) — error bars = std across reps",
                OUT_DIR / "rank_AUC740.png")
    plot_ranked(wide, "AUC_mean_590", "AUC_std_590",
                "Conditions ranked by mean AUC_590 (dye reduction / respiration) — error bars = std across reps",
                OUT_DIR / "rank_AUC590.png")
    plot_ranked(wide, "max_delta_mean_740", "max_delta_std_740",
                "Conditions ranked by mean max delta abs740 (max OD reached above t0)",
                OUT_DIR / "rank_max_delta_740.png")
    plot_ranked(wide, "max_delta_mean_590", "max_delta_std_590",
                "Conditions ranked by mean max delta abs590 (max dye reduction above t0)",
                OUT_DIR / "rank_max_delta_590.png")
    plot_top_kinetics(growth, assign, wide, wavelength=740, n=TOP_N,
                      out=OUT_DIR / "top10_kinetics_740.png")
    plot_top_kinetics(growth, assign, wide, wavelength=590, n=TOP_N,
                      out=OUT_DIR / "top10_kinetics_590.png")
    plot_all_conditions_kinetics(growth, assign, wide, wavelength=740,
                                  out=OUT_DIR / "all_conditions_kinetics_740.png")
    plot_all_conditions_kinetics(growth, assign, wide, wavelength=590,
                                  out=OUT_DIR / "all_conditions_kinetics_590.png")
    plot_growth_rate(wide, wavelength=740, out=OUT_DIR / "growth_rate_740.png")
    plot_growth_rate(wide, wavelength=590, out=OUT_DIR / "growth_rate_590.png")
    plot_growth_rate_top_bottom(wide, wavelength=740,
                                 out=OUT_DIR / "growth_rate_top_bottom_740.png")
    plot_growth_rate_top_bottom(wide, wavelength=590,
                                 out=OUT_DIR / "growth_rate_top_bottom_590.png")

    # per-plate kinetics — one panel per (plate, wavelength), control highlighted
    for plate in sorted(growth["plate"].unique()):
        for wl in (740, 590):
            plot_per_plate_kinetics(growth, assign, plate, wl,
                                    OUT_DIR / f"per_plate_{plate}_{wl}.png")

    # all 4 plates' control overlaid → plate-to-plate variation of ctrl_media
    plot_control_across_plates(growth, assign,
                                OUT_DIR / "control_across_plates.png")

    # ---- report-style outputs (top 5 vs control) ----
    for wl in (740, 590):
        top5 = plot_top_n_vs_control_report(
            growth, assign, wide, n=5, wavelength=wl,
            out_png=OUT_DIR / f"report_top5_vs_control_{wl}.png",
            out_pdf=OUT_DIR / f"report_top5_vs_control_{wl}.pdf",
        )
        imp = improvement_vs_control_table(
            rep_df, top5, wavelength=wl,
            out_csv=OUT_DIR / f"improvement_vs_control_{wl}.csv",
        )
        print(f"\n% improvement of top 5 vs PLATE-MATCHED ctrl_media @ {wl} nm:")
        print(imp.to_string(index=False))

    plot_plate_concordance(rep_df, OUT_DIR / "plate_concordance.png")
    plot_metric_heatmap(wide, OUT_DIR / "heatmap_metrics.png")
    build_condition_viewer(growth, assign, wide,
                           out_path=OUT_DIR / "viewer_conditions.html")

    print(f"\nfigures + tables in {OUT_DIR}")


if __name__ == "__main__":
    main()
