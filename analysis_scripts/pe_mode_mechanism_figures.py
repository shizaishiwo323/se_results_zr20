#!/usr/bin/env python3
"""Generate mechanism figures for Pe-controlled dissolution and SE response."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import seismoelectric_offset_liu2018_spectral as se


RT_BASE = Path(r"C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI")
OUTDIR = ROOT / "paper_figure_sequence_analysis" / "pe_mode_mechanism"

RT_DIRS = {
    "Pe=0.1": RT_BASE
    / "dissolution_results-Da_0.0369_Pe_0.1000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random",
    "Pe=1": RT_BASE
    / "dissolution_results-Da_0.0369_Pe_1.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random",
    "Pe=10": RT_BASE
    / "dissolution_results-Da_0.0369_Pe_10.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random",
}

SE_DIRS = {
    "Pe=0.1": ROOT / "se_results_offset0p1",
    "Pe=1": ROOT / "se_results_offsetPep1",
    "Pe=10": ROOT / "se_results_offset10",
}

COLORS = {
    "Pe=0.1": "#267c7c",
    "Pe=1": "#c57f12",
    "Pe=10": "#9b3a8f",
}


def load_se() -> Dict[str, pd.DataFrame]:
    out = {}
    for pe, path in SE_DIRS.items():
        df = pd.read_csv(path / "seismoelectric_timeseries_results.csv")
        df = df[df["valid_poroelastic"].astype(bool)].copy()
        df["k_k0"] = df["Permeability_mD"] / df["Permeability_mD"].iloc[0]
        df["L_over_sigma"] = df["L_abs"] / df["sigma_abs"]
        out[pe] = df
    return out


def nearest_row(df: pd.DataFrame, phi: float) -> pd.Series:
    return df.loc[(df["Porosity_raw"] - phi).abs().idxmin()]


def coefficient_from_row(
    row: pd.Series,
    *,
    phi: float | None = None,
    k0: float | None = None,
    tau: float | None = None,
    cH: float | None = None,
) -> Dict[str, float]:
    cfg = se.SEConfig()
    coeff = se.se_coefficients(
        float(row["Porosity_used"] if phi is None else phi),
        float(row["k0_m2"] if k0 is None else k0),
        float(row["Tortuosity"] if tau is None else tau),
        float(row["OutletHConc_raw"] if cH is None else cH),
        float(row["omega0_rad_s"]),
        float(row["theta_deg"]),
        cfg,
    )
    return {
        "RE_abs": abs(coeff["R_E"]),
        "TTM_abs": abs(coeff["T_TM"]),
        "L_abs": abs(coeff["L"]),
        "sigma_abs": abs(coeff["sigma"]),
        "L_over_sigma": abs(coeff["L"]) / abs(coeff["sigma"]),
        "zeta": float(coeff["zeta"]),
        "pH": float(coeff["pH"]),
    }


def style_axis(ax, xlabel=None, ylabel=None, title=None, logy=False):
    ax.grid(True, which="major", alpha=0.22, linewidth=0.7)
    ax.grid(True, which="minor", alpha=0.10, linewidth=0.45)
    if logy:
        ax.set_yscale("log")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontweight="bold")


def savefig(fig: plt.Figure, name: str):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / name, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / name.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)


def draw_mechanism_cartoon():
    fig, ax = plt.subplots(figsize=(13.4, 7.8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(x, y, w, h, text, fc, ec="#2f2f2f", size=10.5, weight="normal"):
        rect = plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=1.25)
        ax.add_patch(rect)
        ax.text(
            x + w / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            fontsize=size,
            fontweight=weight,
            wrap=True,
        )
        return rect

    def arrow(x1, y1, x2, y2, color="#333333", lw=1.6, text=None, ty=0.0):
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops=dict(arrowstyle="-|>", lw=lw, color=color, shrinkA=2, shrinkB=2),
        )
        if text:
            ax.text((x1 + x2) / 2, (y1 + y2) / 2 + ty, text, ha="center", va="center", fontsize=9, color=color)

    ax.text(
        0.02,
        0.96,
        "Dissolution-regime control of seismoelectric interface EM response",
        fontsize=15,
        fontweight="bold",
        ha="left",
    )
    ax.text(
        0.02,
        0.915,
        "Main mechanism: channeling enhances hydraulic connectivity, but acid breakthrough collapses electrokinetic conversion.",
        fontsize=10.5,
        ha="left",
        color="#3b3b3b",
    )

    box(0.035, 0.72, 0.22, 0.12, "Low Pe\nfront dominated\nslow H+ breakthrough", "#dff3ef", size=9.2, weight="bold")
    box(0.035, 0.50, 0.22, 0.12, "Intermediate Pe\nchanneling begins\npartial acid breakthrough", "#fff0cf", size=9.2, weight="bold")
    box(0.035, 0.28, 0.22, 0.12, "High Pe\nfast channeling\nacid breakthrough", "#f1e2f2", size=9.2, weight="bold")

    box(0.33, 0.64, 0.19, 0.13, "Hydraulic branch\nk/k0 increases\ntortuosity decreases", "#e8edf6", size=9.2, weight="bold")
    box(0.33, 0.36, 0.19, 0.13, "Electrochemical branch\nH+ increases\nzeta -> 0 / sign change\nsigma increases", "#f7e5dd", size=8.4, weight="bold")

    box(0.62, 0.58, 0.18, 0.12, "Potential amplification\nstronger relative\nfluid-solid motion", "#edf6e9", size=9.0)
    box(0.62, 0.34, 0.18, 0.14, "Dominant suppression\nL/sigma collapses\nelectrical screening\nsource polarity changes", "#f8dede", size=8.4, weight="bold")

    box(0.865, 0.55, 0.12, 0.12, "Large SE\nresponse\nif H+ delayed", "#dff3ef", size=8.8, weight="bold")
    box(0.865, 0.32, 0.12, 0.15, "Weak / reversed\ninterface EM\nresponse\nafter breakthrough", "#f1e2f2", size=8.2, weight="bold")

    arrow(0.255, 0.78, 0.33, 0.70)
    arrow(0.255, 0.56, 0.33, 0.70)
    arrow(0.255, 0.34, 0.33, 0.70)
    arrow(0.255, 0.34, 0.33, 0.43, color="#9b3a2c", text="early H+ arrival", ty=-0.035)
    arrow(0.52, 0.70, 0.62, 0.64)
    arrow(0.52, 0.43, 0.62, 0.41, color="#9b3a2c")
    arrow(0.80, 0.64, 0.865, 0.61)
    arrow(0.80, 0.41, 0.865, 0.40, color="#9b3a2c")

    ax.plot([0.56, 0.56], [0.23, 0.80], color="#4f4f4f", lw=1, ls="--")
    ax.text(0.56, 0.82, "competition", fontsize=10, ha="center", color="#4f4f4f")

    ax.text(0.04, 0.13, "Paper-level claim:", fontsize=11.0, fontweight="bold")
    ax.text(
        0.24,
        0.13,
        "Dissolution can make the medium hydraulically more connected but seismoelectrically less active.",
        fontsize=10.6,
        color="#262626",
    )
    ax.text(
        0.24,
        0.07,
        "The controlling diagnostic is not porosity alone, but the coupled electrokinetic efficiency |L(omega)|/|sigma(omega)|.",
        fontsize=10.2,
        color="#262626",
    )
    savefig(fig, "fig1_mechanism_cartoon.png")


def plot_parameter_bridge(data: Dict[str, pd.DataFrame]):
    fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.6), constrained_layout=True)
    panels = [
        ("k_k0", "Permeability ratio k/k0", True),
        ("cH_molL", "H+ concentration (mol/L)", True),
        ("zeta", "Zeta potential (V)", False),
        ("sigma_abs", "|sigma(omega)| (S/m)", True),
        ("L_over_sigma", "|L(omega)| / |sigma(omega)|", True),
        ("Amax_waveform_spectral", "Waveform peak amplitude", True),
    ]
    for ax, (col, title, logy) in zip(axes.flat, panels):
        for pe, df in data.items():
            ax.plot(
                df["Porosity_raw"],
                df[col],
                marker="o",
                markersize=3.5,
                linewidth=1.8,
                label=pe,
                color=COLORS[pe],
            )
        if col == "zeta":
            ax.axhline(0, color="#333333", lw=1.0, ls="--")
        style_axis(ax, xlabel="Porosity", ylabel=title, title=title, logy=logy)
    axes[0, 0].legend(frameon=False, loc="best")
    fig.suptitle(
        "Same final state, different paths: acid breakthrough controls the electrokinetic bridge",
        fontsize=14,
        fontweight="bold",
    )
    savefig(fig, "fig2_same_porosity_parameter_bridge.png")


def plot_attribution(data: Dict[str, pd.DataFrame]):
    targets = [0.60, 0.80, 0.90]
    labels = ["Porosity", "Permeability", "Tortuosity", "Chemistry", "Structure", "Structure+chem"]
    keys = ["phi", "perm", "tau", "chem", "struct", "all"]
    rows = []
    for phi_target in targets:
        base = nearest_row(data["Pe=0.1"], phi_target)
        target = nearest_row(data["Pe=10"], phi_target)
        base_coeff = coefficient_from_row(base)
        variants = {
            "phi": coefficient_from_row(base, phi=float(target["Porosity_used"])),
            "perm": coefficient_from_row(base, k0=float(target["k0_m2"])),
            "tau": coefficient_from_row(base, tau=float(target["Tortuosity"])),
            "chem": coefficient_from_row(base, cH=float(target["OutletHConc_raw"])),
            "struct": coefficient_from_row(
                base,
                phi=float(target["Porosity_used"]),
                k0=float(target["k0_m2"]),
                tau=float(target["Tortuosity"]),
            ),
            "all": coefficient_from_row(
                base,
                phi=float(target["Porosity_used"]),
                k0=float(target["k0_m2"]),
                tau=float(target["Tortuosity"]),
                cH=float(target["OutletHConc_raw"]),
            ),
        }
        for key, value in variants.items():
            rows.append(
                {
                    "phi_target": phi_target,
                    "component": key,
                    "RE_ratio": value["RE_abs"] / base_coeff["RE_abs"],
                    "L_over_sigma_ratio": value["L_over_sigma"] / base_coeff["L_over_sigma"],
                    "zeta": value["zeta"],
                    "sigma_ratio": value["sigma_abs"] / base_coeff["sigma_abs"],
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "same_phi_pe10_vs_pe0p1_attribution.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6), sharey=True, constrained_layout=True)
    palette = ["#9ecae1", "#6baed6", "#3182bd", "#e6550d", "#756bb1", "#54278f"]
    for ax, phi_target in zip(axes, targets):
        part = out[out["phi_target"] == phi_target].set_index("component").loc[keys]
        x = np.arange(len(keys))
        bars = ax.bar(x, part["RE_ratio"], color=palette, edgecolor="#262626", linewidth=0.6)
        ax.set_yscale("log")
        ax.axhline(1, color="#333333", lw=1, ls="--")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_title(f"phi ~= {phi_target:.2f}", loc="left", fontweight="bold")
        ax.grid(True, axis="y", which="major", alpha=0.25)
        ax.grid(True, axis="y", which="minor", alpha=0.10)
        for b, value in zip(bars, part["RE_ratio"]):
            if value < 1e-2:
                ax.text(
                    b.get_x() + b.get_width() / 2,
                    value * 1.4,
                    f"{value:.1e}",
                    ha="center",
                    va="bottom",
                    fontsize=8.5,
                    rotation=90,
                )
    axes[0].set_ylabel("R_E ratio relative to Pe=0.1 baseline")
    fig.suptitle(
        "One-at-a-time attribution: Pe=10 chemistry, not permeability alone, collapses conversion",
        fontsize=14,
        fontweight="bold",
    )
    savefig(fig, "fig3_same_phi_parameter_attribution.png")


def plot_spatial_polarity():
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.7), constrained_layout=True)
    max_rows = []
    for pe, path in SE_DIRS.items():
        df = pd.read_csv(path / "waveform_spatial_peak_diagnostics.csv")
        norm = df["peak_signed"] / df["peak_abs"].abs().max()
        axes[0].plot(df["z_mm"], norm, color=COLORS[pe], lw=1.9, label=pe)
        for side in ["R_E", "T_E"]:
            sub = df[df["side"] == side]
            row = sub.loc[sub["peak_abs"].idxmax()]
            max_rows.append(
                {
                    "Pe": pe,
                    "side": side,
                    "peak_abs": row["peak_abs"],
                    "peak_signed": row["peak_signed"],
                    "z_mm": row["z_mm"],
                }
            )
    axes[0].axhline(0, color="#333333", lw=1)
    axes[0].axvline(0, color="#333333", lw=0.8, ls="--")
    style_axis(
        axes[0],
        xlabel="Receiver z (mm)",
        ylabel="Normalized signed peak",
        title="A. Signed spatial pattern",
    )
    axes[0].legend(frameon=False, loc="best")

    peaks = pd.DataFrame(max_rows)
    width = 0.26
    x = np.arange(3)
    pe_order = ["Pe=0.1", "Pe=1", "Pe=10"]
    for j, side in enumerate(["R_E", "T_E"]):
        values = [abs(peaks[(peaks["Pe"] == pe) & (peaks["side"] == side)]["peak_signed"].iloc[0]) for pe in pe_order]
        signs = [np.sign(peaks[(peaks["Pe"] == pe) & (peaks["side"] == side)]["peak_signed"].iloc[0]) for pe in pe_order]
        bars = axes[1].bar(x + (j - 0.5) * width, values, width=width, label=side, color=["#4c78a8", "#f58518"][j])
        for bar, sign in zip(bars, signs):
            axes[1].text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() * 1.12,
                "+" if sign > 0 else "-",
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
            )
    axes[1].set_yscale("log")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(pe_order)
    style_axis(axes[1], ylabel="Absolute signed peak", title="B. Peak magnitude and polarity")
    axes[1].legend(frameon=False)
    fig.suptitle("Finite-offset polarity reversal records the zeta-potential sign change", fontsize=14, fontweight="bold")
    peaks.to_csv(OUTDIR / "spatial_peak_polarity_summary.csv", index=False)
    savefig(fig, "fig4_spatial_polarity_reversal.png")


def build_snapshot_triptych():
    snapshots = {
        "Pe=0.1\nfront-dominated\nphi~0.79, t=7119.2 s": RT_DIRS["Pe=0.1"] / "timestep_0028.png",
        "Pe=1\nchanneling\nphi~0.80, t=400.6 s": RT_DIRS["Pe=1"] / "timestep_0024.png",
        "Pe=10\nfast acid breakthrough\nphi~0.81, t=56.4 s": RT_DIRS["Pe=10"] / "timestep_0017.png",
    }
    imgs = [Image.open(p).convert("RGB") for p in snapshots.values()]
    target_w = 760
    resized = []
    for img in imgs:
        ratio = target_w / img.width
        resized.append(img.resize((target_w, int(img.height * ratio)), Image.Resampling.LANCZOS))
    header_h = 92
    pad = 24
    width = target_w * 3 + pad * 4
    height = max(img.height for img in resized) + header_h + pad * 2
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font_title = ImageFont.truetype("arial.ttf", 28)
        font_sub = ImageFont.truetype("arial.ttf", 22)
    except OSError:
        font_title = ImageFont.load_default()
        font_sub = ImageFont.load_default()
    draw.text((pad, 18), "Dissolution modes at comparable porosity", fill=(30, 30, 30), font=font_title)
    for i, (label, img) in enumerate(zip(snapshots.keys(), resized)):
        x = pad + i * (target_w + pad)
        y = header_h
        canvas.paste(img, (x, y))
        draw.text((x + 10, y + 8), label, fill=(20, 20, 20), font=font_sub)
        draw.rectangle((x, y, x + target_w, y + img.height), outline=(70, 70, 70), width=2)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    canvas.save(OUTDIR / "fig5_dissolution_mode_snapshots_phi080.png")


def write_snapshot_summary(data: Dict[str, pd.DataFrame]):
    rows = []
    for pe, df in data.items():
        row = nearest_row(df, 0.75)
        rows.append(
            {
                "Pe": pe,
                "Time_s": row["Time_s"],
                "Porosity_raw": row["Porosity_raw"],
                "k_k0": row["k_k0"],
                "cH_molL": row["cH_molL"],
                "zeta": row["zeta"],
                "sigma_abs": row["sigma_abs"],
                "L_over_sigma": row["L_over_sigma"],
                "RE_abs": row["RE_abs"],
                "TTM_abs": row["TTM_abs"],
                "Amax_waveform_spectral": row["Amax_waveform_spectral"],
            }
        )
    pd.DataFrame(rows).to_csv(OUTDIR / "same_phi075_se_summary.csv", index=False)


def main():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    OUTDIR.mkdir(parents=True, exist_ok=True)
    data = load_se()
    write_snapshot_summary(data)
    draw_mechanism_cartoon()
    plot_parameter_bridge(data)
    plot_attribution(data)
    plot_spatial_polarity()
    build_snapshot_triptych()
    print("Wrote Pe-mode mechanism figures.")


if __name__ == "__main__":
    main()
