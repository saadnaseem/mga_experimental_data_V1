"""
Biolog PM08 SN1 1c / 2c — quick-look growth-kinetics plotter.

Inputs (relative to this script):
    SN1_PlateInfo.csv
    SN1_RawReads/SN1_PM08_*_SN1_1c_*_RawReads.csv
    SN1_RawReads/SN1_PM08_*_SN1_2c_*_RawReads.csv

Outputs (./outputs/):
    growth_long.csv               tidy long: plate, well, row, col, wavelength, t_h, abs
    grid_<plate>.png              8x12 small-multiples overview, both wavelengths
    viewer.html                   interactive: pick wells via checkboxes / dropdowns

Run:
    python plot_growth.py
Then open outputs/viewer.html in a browser.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots


HERE = Path(__file__).resolve().parent
RAW_DIR = HERE / "SN1_RawReads"
OUT_DIR = HERE / "outputs"
OUT_DIR.mkdir(exist_ok=True)

WELL_RE = re.compile(r"^[A-H](0[1-9]|1[0-2])$")
ROWS = list("ABCDEFGH")
COLS = [f"{c:02d}" for c in range(1, 13)]
WELLS = [f"{r}{c}" for r in ROWS for c in COLS]

WL_COLORS = {590: "tab:red", 740: "tab:blue"}  # 590 = dye reduction, 740 = turbidity


# ---------- load + reshape ----------

def load_plate_info() -> pd.DataFrame:
    return pd.read_csv(HERE / "SN1_PlateInfo.csv")


def load_raw(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    well_cols = [c for c in df.columns if WELL_RE.match(c)]
    df["Read At"] = pd.to_datetime(df["Read At"])
    plate_id = df["PlateId"].iloc[0]

    long = df.melt(
        id_vars=["PlateId", "Wavelength", "Read At"],
        value_vars=well_cols,
        var_name="well",
        value_name="abs",
    )
    t0 = long["Read At"].min()
    long["t_h"] = (long["Read At"] - t0).dt.total_seconds() / 3600.0
    long["plate_uuid"] = plate_id
    long["row"] = long["well"].str[0]
    long["col"] = long["well"].str[1:].astype(int)
    long = long.rename(columns={"Wavelength": "wavelength"})
    return long[["plate_uuid", "well", "row", "col", "wavelength", "t_h", "abs"]]


def build_long_table() -> pd.DataFrame:
    info = load_plate_info()
    uuid_to_name = dict(zip(info["PlateId"], info["Sample"]))

    frames = []
    for f in sorted(RAW_DIR.glob("*_RawReads.csv")):
        long = load_raw(f)
        long["plate"] = long["plate_uuid"].map(uuid_to_name).fillna(long["plate_uuid"])
        frames.append(long)
    long = pd.concat(frames, ignore_index=True)
    return long.sort_values(["plate", "well", "wavelength", "t_h"]).reset_index(drop=True)


# ---------- static 8x12 grid ----------

def plot_grid(long: pd.DataFrame, plate: str, out_path: Path) -> None:
    sub = long[long["plate"] == plate]
    y_max = sub["abs"].max() * 1.05
    y_min = sub["abs"].min() * 0.95

    fig, axes = plt.subplots(8, 12, figsize=(20, 12), sharex=True, sharey=True)
    fig.suptitle(f"Plate {plate} — abs590 (red) / abs740 (blue) vs time (h)", fontsize=14)

    for i, row in enumerate(ROWS):
        for j, col in enumerate(COLS):
            ax = axes[i, j]
            well = f"{row}{col}"
            for wl, color in WL_COLORS.items():
                w = sub[(sub["well"] == well) & (sub["wavelength"] == wl)]
                if not w.empty:
                    ax.plot(w["t_h"].to_numpy(), w["abs"].to_numpy(),
                            color=color, lw=0.9)
            ax.set_title(well, fontsize=7, pad=1)
            ax.set_ylim(y_min, y_max)
            ax.tick_params(labelsize=5)
            if i != 7:
                ax.set_xticklabels([])
            if j != 0:
                ax.set_yticklabels([])

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------- interactive viewer ----------

def build_viewer(long: pd.DataFrame, out_path: Path) -> None:
    plates = sorted(long["plate"].unique())
    fig = make_subplots(
        rows=1, cols=len(plates),
        subplot_titles=[f"{p} — abs740 solid / abs590 dashed" for p in plates],
        shared_yaxes=True,
        horizontal_spacing=0.05,
    )

    # one trace per (plate, well, wavelength); legendgroup=well so 590+740 toggle together
    trace_meta: list[tuple[int, str, int]] = []  # (plate_subplot_idx, well, wavelength)
    for p_idx, plate in enumerate(plates, start=1):
        sub = long[long["plate"] == plate]
        for well in WELLS:
            for wl in (740, 590):
                w = sub[(sub["well"] == well) & (sub["wavelength"] == wl)].sort_values("t_h")
                if w.empty:
                    continue
                visible = (well == "A01")  # start with A01 only
                fig.add_trace(
                    go.Scatter(
                        x=w["t_h"].to_numpy(),
                        y=w["abs"].to_numpy(),
                        name=f"{plate} {well} {wl}",
                        legendgroup=well,
                        legendgrouptitle_text=well if (p_idx == 1 and wl == 740) else None,
                        line=dict(width=1.2, dash=("dash" if wl == 590 else "solid")),
                        mode="lines",
                        visible=visible,
                        hovertemplate=f"{plate} {well} {wl}nm<br>t=%{{x:.2f}} h<br>abs=%{{y:.3f}}<extra></extra>",
                    ),
                    row=1, col=p_idx,
                )
                trace_meta.append((p_idx, well, wl))

    n = len(trace_meta)

    def vis_mask(predicate) -> list[bool]:
        return [predicate(meta) for meta in trace_meta]

    well_buttons = [
        dict(label="All on",   method="update", args=[{"visible": [True] * n}]),
        dict(label="All off",  method="update", args=[{"visible": [False] * n}]),
        dict(label="A01 only", method="update",
             args=[{"visible": vis_mask(lambda m: m[1] == "A01")}]),
    ]
    for r in ROWS:
        well_buttons.append(dict(
            label=f"Row {r}", method="update",
            args=[{"visible": vis_mask(lambda m, rr=r: m[1].startswith(rr))}],
        ))
    for c in COLS:
        well_buttons.append(dict(
            label=f"Col {c}", method="update",
            args=[{"visible": vis_mask(lambda m, cc=c: m[1].endswith(cc))}],
        ))

    # wavelength dropdown — uses targeted restyle so it composes with the well filter
    idx_740 = [i for i, m in enumerate(trace_meta) if m[2] == 740]
    idx_590 = [i for i, m in enumerate(trace_meta) if m[2] == 590]
    wl_buttons = [
        dict(label="Show 740", method="restyle", args=[{"visible": True},  idx_740]),
        dict(label="Hide 740", method="restyle", args=[{"visible": False}, idx_740]),
        dict(label="Show 590", method="restyle", args=[{"visible": True},  idx_590]),
        dict(label="Hide 590", method="restyle", args=[{"visible": False}, idx_590]),
    ]

    fig.update_layout(
        height=650,
        title="SN1 Biolog PM08 — interactive viewer (legend = wells, dropdowns = quick picks)",
        legend=dict(groupclick="togglegroup", itemclick="toggle", itemdoubleclick="toggleothers"),
        updatemenus=[
            dict(type="dropdown", buttons=well_buttons,
                 x=0.00, y=1.13, xanchor="left", yanchor="top",
                 showactive=False, name="wells"),
            dict(type="dropdown", buttons=wl_buttons,
                 x=0.18, y=1.13, xanchor="left", yanchor="top",
                 showactive=False, name="wavelength"),
        ],
        annotations=list(fig.layout.annotations) + [
            dict(text="wells ▾",      x=0.00, y=1.18, xref="paper", yref="paper",
                 showarrow=False, xanchor="left", font=dict(size=11)),
            dict(text="wavelength ▾", x=0.18, y=1.18, xref="paper", yref="paper",
                 showarrow=False, xanchor="left", font=dict(size=11)),
        ],
        margin=dict(l=60, r=20, t=140, b=50),
    )
    for p_idx in range(1, len(plates) + 1):
        fig.update_xaxes(title_text="time (h)", row=1, col=p_idx)
    fig.update_yaxes(title_text="absorbance", row=1, col=1)

    # include_plotlyjs=True embeds plotly.js so the file works offline / behind
    # firewalls. Switch to "cdn" for a smaller file if internet is always available.
    fig.write_html(out_path, include_plotlyjs=True, full_html=True)


# ---------- main ----------

def main() -> None:
    long = build_long_table()
    out_csv = OUT_DIR / "growth_long.csv"
    long.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}  ({len(long):,} rows)")

    for plate in sorted(long["plate"].unique()):
        out = OUT_DIR / f"grid_{plate}.png"
        plot_grid(long, plate, out)
        print(f"wrote {out}")

    viewer = OUT_DIR / "viewer.html"
    build_viewer(long, viewer)
    print(f"wrote {viewer}  (open in browser)")


if __name__ == "__main__":
    main()
