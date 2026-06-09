#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared analysis routines for the recommended manuscript figure sequence.

The routines in this file keep the paper-facing analyses tied to the current
Schakel/Pride/Liu forward model in ``seismoelectric_offset_liu2018_spectral``.
Each figure script is intentionally thin; this module owns data loading,
common diagnostics, and plotting.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import liu2018_fig2b_dipole_comparison as liu_fig2b
import parameter_sensitivity_analysis as psa
import seismoelectric_offset_liu2018_spectral as se


DEFAULT_OUTDIR = PROJECT_ROOT / "paper_figure_sequence_analysis"
IMAGEGEN_SOURCE = Path(
    r"C:\Users\imgw\.codex\generated_images\019e989d-7c5c-7e11-b6e3-788ce25ff052"
    r"\ig_067344ae93211b16016a22faa788cc8191b8a7d30ff452a4f8.png"
)


@dataclass(frozen=True)
class DissolutionCase:
    name: str
    pe: float
    path: Path
    color: str

    @property
    def input_path(self) -> Path:
        return self.path / "global_evolution.xlsx"

    @property
    def tortuosity_path(self) -> Path:
        return self.path / "tortuosity_segments.xlsx"


DEFAULT_CASES: Tuple[DissolutionCase, ...] = (
    DissolutionCase(
        "Pe0p1",
        0.1,
        Path(
            r"C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI"
            r"\dissolution_results-Da_0.0369_Pe_0.1000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random"
        ),
        "#2f6f9f",
    ),
    DissolutionCase(
        "Pe1",
        1.0,
        Path(
            r"C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI"
            r"\dissolution_results-Da_0.0369_Pe_1.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random"
        ),
        "#3f8f5f",
    ),
    DissolutionCase(
        "Pe10",
        10.0,
        Path(
            r"C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI"
            r"\dissolution_results-Da_0.0369_Pe_10.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random"
        ),
        "#b05a32",
    ),
)


def default_cfg() -> se.SEConfig:
    cfg = se.SEConfig()
    # The displayed waveform window is about 33 us wide.  A coarse uniform
    # frequency grid creates inverse-transform replicas with period 1/df; 128
    # samples push those aliases outside the T0 audit window for Figure 5.
    cfg.spectral_n_omega = 128
    cfg.spectral_n_k = 81
    cfg.waveform_nt = 720
    return cfg


def ensure_outdir(outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    for sub in ["tables", "waveform_cache", "figures", "discussion"]:
        (outdir / sub).mkdir(parents=True, exist_ok=True)
    return outdir


def load_case_table(case: DissolutionCase, cfg: se.SEConfig) -> pd.DataFrame:
    df = se.load_reactive_transport_table(case.input_path)
    ts = se.compute_time_series(df, cfg)
    out = df.merge(
        ts[
            [
                "Time_s",
                "Time_min",
                "Porosity_raw",
                "Porosity_used",
                "valid_poroelastic",
                "C_molL",
                "pH",
                "omega_t",
                "Lambda_m",
                "L_abs",
                "sigma_abs",
                "RE_real",
                "RE_imag",
                "RE_abs",
                "TTM_real",
                "TTM_imag",
                "TTM_abs",
                "matrix_cond",
            ]
        ],
        on="Time_s",
        how="left",
    )
    out["case"] = case.name
    out["Pe"] = case.pe
    out["k0_m2"] = out["Permeability_mD"].astype(float) * 9.869233e-16
    out["InjectedPV_norm"] = out["InjectedPV"] / out["InjectedPV"].iloc[-1]
    out["Time_norm"] = out["Time_s"] / out["Time_s"].iloc[-1]
    out["PoreVolumeToSurface_m"] = out.apply(psa.pore_volume_to_surface_ratio_m, axis=1)
    return out


def load_all_cases(cases: Sequence[DissolutionCase], cfg: se.SEConfig) -> Dict[str, pd.DataFrame]:
    return {case.name: load_case_table(case, cfg) for case in cases}


def valid_rows(df: pd.DataFrame, cfg: se.SEConfig) -> pd.DataFrame:
    return df[(df["Porosity"] > cfg.phi_min) & (df["Porosity"] < cfg.phi_max_valid)].copy()


def representative_indices(df: pd.DataFrame, cfg: se.SEConfig, n: int = 4) -> List[int]:
    valid = valid_rows(df, cfg)
    if valid.empty:
        return []
    targets = np.linspace(0.0, 1.0, n)
    idx: List[int] = []
    x_raw = valid["Time_s"].to_numpy(float)
    span = float(np.nanmax(x_raw) - np.nanmin(x_raw))
    if np.isfinite(span) and span > 0:
        x = (x_raw - np.nanmin(x_raw)) / span
    else:
        x = np.linspace(0.0, 1.0, len(valid))
    for target in targets:
        local_idx = int(np.nanargmin(np.abs(x - target)))
        idx.append(int(valid.index[local_idx]))
    return sorted(set(idx))


def row_dynamic_frequency_scan(
    row: pd.Series,
    cfg: se.SEConfig,
    freqs_hz: np.ndarray,
) -> pd.DataFrame:
    phi = float(np.clip(float(row["Porosity"]), cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
    tau = max(float(row["Tortuosity"]), 1.0 + 1e-6)
    cH = float(row["OutletHConc"])
    C_override = se.optional_float(row, "ElectrolyteConcentration_molL")
    sigma_f_override = se.optional_float(row, "FluidConductivity_S_m")
    rows = []
    for f in freqs_hz:
        dyn = se.dynamic_coefficients(
            phi,
            k0_m2,
            tau,
            cH,
            2.0 * math.pi * float(f),
            cfg,
            C_override_molL=C_override,
            sigma_f_override=sigma_f_override,
        )
        rows.append(
            {
                "Time_s": float(row["Time_s"]),
                "Porosity": float(row["Porosity"]),
                "frequency_Hz": float(f),
                "k_dyn_abs": abs(dyn["k_dyn"]),
                "L_abs": abs(dyn["L"]),
                "sigma_abs": abs(dyn["sigma"]),
                "omega_t": float(dyn["omega_t"]),
                "sigma_f": float(dyn["sigma_f"]),
                "C_molL": float(dyn["C_molL"]),
                "pH": float(dyn["pH"]),
            }
        )
    return pd.DataFrame(rows)


def coefficient_record_for_row(
    row: pd.Series,
    cfg: se.SEConfig,
    freq_hz: float,
    theta_deg: float,
) -> Dict[str, float]:
    phi = float(np.clip(float(row["Porosity"]), cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
    tau = max(float(row["Tortuosity"]), 1.0 + 1e-6)
    cH = float(row["OutletHConc"])
    coeff = se.se_coefficients(
        phi,
        k0_m2,
        tau,
        cH,
        2.0 * math.pi * float(freq_hz),
        theta_deg,
        cfg,
        C_override_molL=se.optional_float(row, "ElectrolyteConcentration_molL"),
        sigma_f_override=se.optional_float(row, "FluidConductivity_S_m"),
    )
    return {
        "frequency_Hz": float(freq_hz),
        "theta_deg": float(theta_deg),
        "Time_s": float(row["Time_s"]),
        "Porosity": float(row["Porosity"]),
        "RE_abs": abs(coeff["R_E"]),
        "TTM_abs": abs(coeff["T_TM"]),
        "RE_phase_rad": float(np.angle(coeff["R_E"])),
        "TTM_phase_rad": float(np.angle(coeff["T_TM"])),
        "L_abs": abs(coeff["L"]),
        "sigma_abs": abs(coeff["sigma"]),
        "matrix_cond": float(coeff["matrix_cond"]),
    }


def synthesize_or_load_waveform(
    case: DissolutionCase,
    row: pd.Series,
    cfg: se.SEConfig,
    outdir: Path,
    n_omega: int | None = None,
    n_k: int = 81,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Path]:
    n_omega = int(n_omega or cfg.spectral_n_omega)
    cache_dir = outdir / "waveform_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    name = f"{case.name}_step{int(row['TimeStep']):04d}_nw{n_omega}_nk{n_k}.npz"
    cache_path = cache_dir / name
    if cache_path.exists():
        arr = np.load(cache_path)
        return arr["z_m"], arr["t_s"], arr["U"], cache_path
    z, t, U = se.synthesize_waveforms_spectral(row, cfg, n_omega=n_omega, n_k=n_k)
    np.savez_compressed(
        cache_path,
        z_m=z,
        t_s=t,
        U=U,
        Time_s=float(row["Time_s"]),
        TimeStep=int(row["TimeStep"]),
        Pe=float(case.pe),
    )
    return z, t, U, cache_path


def post_t0_waveform_peak_metric(row: pd.Series, cfg: se.SEConfig, n_omega: int, n_k: int) -> float:
    z, t, U = se.synthesize_waveforms_spectral(row, cfg, n_omega=n_omega, n_k=n_k)
    del z
    T0 = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)
    mask = t >= T0
    val = float(np.nanmax(np.abs(U[:, mask]))) if np.any(mask) else np.nan
    return val if np.isfinite(val) else np.nan


def copy_imagegen_asset(outdir: Path) -> Path | None:
    if not IMAGEGEN_SOURCE.exists():
        return None
    target = outdir / "figures" / "figure1_imagegen_concept.png"
    shutil.copy2(IMAGEGEN_SOURCE, target)
    return target


def make_figure1(outdir: Path) -> Path:
    ensure_outdir(outdir)
    copy_imagegen_asset(outdir)
    fig, ax = plt.subplots(figsize=(12.5, 6.6), constrained_layout=True)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    def box(xy, wh, label, face, edge="#374151", fontsize=10):
        x, y = xy
        w, h = wh
        patch = plt.Rectangle((x, y), w, h, facecolor=face, edgecolor=edge, linewidth=1.2)
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fontsize)

    ax.text(0.06, 0.92, "Slow dissolution time", fontsize=15, weight="bold", color="#1f2937")
    ax.text(0.66, 0.92, "Fast waveform time", fontsize=15, weight="bold", color="#1f2937")
    ax.annotate(
        "Time_s: seconds to hours",
        xy=(0.08, 0.86),
        xytext=(0.36, 0.86),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="#2f6f9f"),
        fontsize=10,
        color="#2f6f9f",
        ha="center",
    )
    ax.annotate(
        "microseconds after acoustic excitation",
        xy=(0.72, 0.86),
        xytext=(0.94, 0.86),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="#b05a32"),
        fontsize=10,
        color="#b05a32",
        ha="center",
    )

    box((0.05, 0.57), (0.25, 0.19), "Reactive transport outputs\nphi, k0, alpha_inf\nH+, sigma_f", "#eaf3f9")
    box((0.05, 0.29), (0.25, 0.19), "Pore geometry and connectivity\nporosity, permeability,\ntortuosity, surface area", "#eaf3f9")
    box((0.39, 0.58), (0.23, 0.17), "Dynamic electrokinetic bridge\nk(omega), L(omega), sigma(omega)", "#eef8ef")
    box((0.39, 0.31), (0.23, 0.17), "Interface boundary problem\nSchakel coefficients\nR_E, T_TM, phase", "#eef8ef")
    box((0.72, 0.58), (0.22, 0.17), "Liu finite-offset waveform\nfrequency-wavenumber integral", "#fff3e8")
    box((0.72, 0.31), (0.22, 0.17), "Interface EM response\namplitude, polarity,\nreceiver-position pattern", "#fff3e8")

    for start, end in [
        ((0.30, 0.665), (0.39, 0.665)),
        ((0.30, 0.385), (0.39, 0.385)),
        ((0.505, 0.58), (0.505, 0.48)),
        ((0.62, 0.665), (0.72, 0.665)),
        ((0.62, 0.395), (0.72, 0.395)),
        ((0.83, 0.58), (0.83, 0.48)),
    ]:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.4, color="#4b5563"))

    ax.plot([0.05, 0.94], [0.18, 0.18], color="#9ca3af", lw=1.0)
    ax.text(
        0.50,
        0.12,
        "The dissolution clock changes material properties; each snapshot is then used for a separate microsecond-scale waveform simulation.",
        ha="center",
        fontsize=10.5,
        color="#374151",
    )
    outpath = outdir / "figures" / "figure1_conceptual_workflow_timescale_separation.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def make_figure2(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    tables = load_all_cases(cases, cfg)
    all_rows = pd.concat(tables.values(), ignore_index=True)
    all_rows.to_csv(outdir / "tables" / "figure2_reactive_transport_evolution.csv", index=False)
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 9.0), constrained_layout=True)
    panels = [
        ("Time_norm", "Porosity", "porosity", "linear"),
        ("Time_norm", "Permeability_mD", "permeability (mD)", "log"),
        ("Time_norm", "Tortuosity", "tortuosity", "linear"),
        ("Time_norm", "pH", "pH from OutletHConc", "linear"),
        ("Time_norm", "C_molL", "electrolyte concentration (mol/L)", "log"),
        ("Time_norm", "PoreVolumeToSurface_m", "pore volume / surface area (m)", "log"),
    ]
    for ax, (_, ycol, ylabel, yscale) in zip(axes.ravel(), panels):
        for case in cases:
            df = tables[case.name]
            ax.plot(df["Time_norm"], df[ycol], color=case.color, lw=1.6, marker=".", ms=3, label=f"Pe={case.pe:g}")
        ax.set_xlabel("normalized dissolution time")
        ax.set_ylabel(ylabel)
        if yscale == "log":
            ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=9)
    fig.suptitle("Figure 2. Reactive-transport and hydrological/electrical evolution")
    outpath = outdir / "figures" / "figure2_reactive_transport_evolution.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def make_figure3(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    freqs = np.logspace(math.log10(0.2 * cfg.f0), math.log10(3.0 * cfg.f0), 80)
    rows = []
    for case in cases:
        df = load_case_table(case, cfg)
        for idx in representative_indices(df, cfg, n=4):
            rec = row_dynamic_frequency_scan(df.loc[idx], cfg, freqs)
            rec["case"] = case.name
            rec["Pe"] = case.pe
            rec["Time_norm"] = float(df.loc[idx, "Time_norm"])
            rows.append(rec)
    scan = pd.concat(rows, ignore_index=True)
    scan.to_csv(outdir / "tables" / "figure3_dynamic_electrokinetic_bridge.csv", index=False)

    fig, axes = plt.subplots(3, len(cases), figsize=(13.5, 8.5), sharex=True, constrained_layout=True)
    metrics = [("k_dyn_abs", "|k(omega)| (m2)"), ("L_abs", "|L(omega)|"), ("sigma_abs", "|sigma(omega)| (S/m)")]
    for col, case in enumerate(cases):
        part_case = scan[scan["case"] == case.name]
        for row_i, (metric, ylabel) in enumerate(metrics):
            ax = axes[row_i, col]
            for _, part in part_case.groupby("Time_s"):
                label = f"t/T={part['Time_norm'].iloc[0]:.2f}"
                ax.plot(part["frequency_Hz"] / 1e3, part[metric], lw=1.2, label=label)
            ax.axvline(cfg.f0 / 1e3, color="0.25", ls=":", lw=1.0)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.22)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row_i == 0:
                ax.set_title(f"Pe={case.pe:g}")
            if row_i == 2:
                ax.set_xlabel("frequency (kHz)")
    axes[0, 0].legend(fontsize=7, loc="best")
    fig.suptitle("Figure 3. Frequency-dependent dynamic electrokinetic bridge")
    outpath = outdir / "figures" / "figure3_dynamic_electrokinetic_bridge.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def make_figure4(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    time_frames = []
    heat_frames = []
    freqs = np.linspace(0.4 * cfg.f0, 2.0 * cfg.f0, 24)
    angles = np.linspace(5.0, 75.0, 24)
    for case in cases:
        df = load_case_table(case, cfg)
        time_part = df[["Time_s", "Time_norm", "RE_abs", "TTM_abs", "RE_real", "RE_imag", "TTM_real", "TTM_imag"]].copy()
        time_part["case"] = case.name
        time_part["Pe"] = case.pe
        time_part["RE_phase_rad"] = np.angle(time_part["RE_real"] + 1j * time_part["RE_imag"])
        time_part["TTM_phase_rad"] = np.angle(time_part["TTM_real"] + 1j * time_part["TTM_imag"])
        time_frames.append(time_part)
        for idx in representative_indices(df, cfg, n=3):
            for f in freqs:
                for angle in angles:
                    rec = coefficient_record_for_row(df.loc[idx], cfg, f, angle)
                    rec["case"] = case.name
                    rec["Pe"] = case.pe
                    rec["stage_Time_norm"] = float(df.loc[idx, "Time_norm"])
                    heat_frames.append(rec)
    time_scan = pd.concat(time_frames, ignore_index=True)
    heat = pd.DataFrame(heat_frames)
    time_scan.to_csv(outdir / "tables" / "figure4_conversion_time_evolution.csv", index=False)
    heat.to_csv(outdir / "tables" / "figure4_frequency_angle_heatmap.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), constrained_layout=True)
    for case in cases:
        part = time_scan[time_scan["case"] == case.name]
        axes[0, 0].plot(part["Time_norm"], part["RE_abs"] / part["RE_abs"].iloc[0], color=case.color, label=f"Pe={case.pe:g}")
        axes[0, 1].plot(part["Time_norm"], part["TTM_abs"] / part["TTM_abs"].iloc[0], color=case.color)
        axes[1, 0].plot(part["Time_norm"], np.unwrap(part["RE_phase_rad"]), color=case.color)
        axes[1, 1].plot(part["Time_norm"], np.unwrap(part["TTM_phase_rad"]), color=case.color)
    for ax, ylabel in zip(
        axes.ravel(),
        ["|R_E| normalized", "|T_TM| normalized", "R_E phase (rad)", "T_TM phase (rad)"],
    ):
        ax.set_xlabel("normalized dissolution time")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
    axes[0, 0].set_yscale("log")
    axes[0, 1].set_yscale("log")
    axes[0, 0].legend(loc="best", fontsize=9)
    fig.suptitle("Figure 4a. Interface conversion coefficient evolution")
    fig.savefig(outdir / "figures" / "figure4a_conversion_coefficients_time.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(len(cases), 3, figsize=(12.0, 9.0), constrained_layout=True)
    for row_i, case in enumerate(cases):
        part_case = heat[heat["case"] == case.name]
        for col_i, stage in enumerate(sorted(part_case["stage_Time_norm"].unique())):
            ax = axes[row_i, col_i]
            part = part_case[np.isclose(part_case["stage_Time_norm"], stage)]
            grid = part.pivot_table(index="theta_deg", columns="frequency_Hz", values="RE_abs")
            im = ax.imshow(
                np.log10(grid.to_numpy()),
                origin="lower",
                aspect="auto",
                extent=[freqs.min() / 1e3, freqs.max() / 1e3, angles.min(), angles.max()],
                cmap="viridis",
            )
            ax.set_title(f"Pe={case.pe:g}, t/T={stage:.2f}")
            ax.set_xlabel("frequency (kHz)")
            if col_i == 0:
                ax.set_ylabel("incidence angle (deg)")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="log10 |R_E|")
    fig.suptitle("Figure 4b. Schakel-style frequency-angle conversion heatmaps")
    outpath = outdir / "figures" / "figure4b_conversion_frequency_angle_heatmaps.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def make_figure5(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    summary_rows = []
    for stale in (outdir / "tables").glob("figure5_*_peak_diagnostics.csv"):
        stale.unlink()
    fig, axes = plt.subplots(len(cases), 3, figsize=(13.0, 9.0), constrained_layout=True)
    for row_i, case in enumerate(cases):
        df = load_case_table(case, cfg)
        idxs = representative_indices(df, cfg, n=3)
        for col_i, idx in enumerate(idxs):
            row = df.loc[idx]
            z, t, U, cache_path = synthesize_or_load_waveform(case, row, cfg, outdir)
            diag = se.waveform_spatial_peak_diagnostics(z, t, U, cfg.offset_D)
            diag.to_csv(outdir / "tables" / f"figure5_{case.name}_step{int(row['TimeStep']):04d}_peak_diagnostics.csv", index=False)
            t_us = t * 1e6
            z_mm = z * 1e3
            U_plot = U.copy()
            vmax = np.nanmax(np.abs(U_plot))
            if np.isfinite(vmax) and vmax > 0:
                U_plot = U_plot / vmax
            ax = axes[row_i, col_i]
            im = ax.imshow(
                U_plot,
                origin="lower",
                aspect="auto",
                extent=[t_us.min(), t_us.max(), z_mm.min(), z_mm.max()],
                cmap="RdBu_r",
                vmin=-1,
                vmax=1,
            )
            T0_us = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl) * 1e6
            pre = float(np.nanmax(np.abs(U[:, t < T0_us * 1e-6])))
            post = float(np.nanmax(np.abs(U[:, t >= T0_us * 1e-6])))
            ax.axvline(T0_us, color="k", ls=":", lw=1.0)
            ax.axhline(0, color="k", lw=0.7)
            ax.set_title(f"Pe={case.pe:g}, t/T={row['Time_norm']:.2f}")
            ax.set_xlabel("waveform time (us)")
            if col_i == 0:
                ax.set_ylabel("receiver z (mm)")
            summary_rows.append(
                {
                    "case": case.name,
                    "Pe": case.pe,
                    "TimeStep": int(row["TimeStep"]),
                    "Time_s": float(row["Time_s"]),
                    "Time_norm": float(row["Time_norm"]),
                    "cache_path": str(cache_path),
                    "Amax": float(np.nanmax(np.abs(U))),
                    "pre_T0_max_abs": pre,
                    "post_T0_max_abs": post,
                    "pre_to_post_T0_ratio": pre / post if post > 0 else np.nan,
                }
            )
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.025, pad=0.01, label="normalized potential")
    fig.suptitle("Figure 5. Finite-offset waveform panels")
    outpath = outdir / "figures" / "figure5_finite_offset_waveform_panels.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    pd.DataFrame(summary_rows).to_csv(outdir / "tables" / "figure5_waveform_panel_summary.csv", index=False)
    return outpath


def make_figure6(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    rows = []
    fig, axes = plt.subplots(1, len(cases), figsize=(13.0, 4.2), constrained_layout=True)
    for ax, case in zip(axes, cases):
        df = load_case_table(case, cfg)
        idx = representative_indices(df, cfg, n=3)[-1]
        row = df.loc[idx]
        z, t, U, _ = synthesize_or_load_waveform(case, row, cfg, outdir)
        table = liu_fig2b.build_comparison_table(z, U, cfg, prefix="modeled")
        table["case"] = case.name
        table["Pe"] = case.pe
        table["Time_s"] = float(row["Time_s"])
        table["Time_norm"] = float(row["Time_norm"])
        rows.append(table)
        for side, color in [("R_E", "#2f6f9f"), ("T_E", "#b05a32")]:
            mask = np.where(z < 0, "R_E", np.where(z > 0, "T_E", "interface")) == side
            part = table[mask]
            ax.plot(
                np.abs(part["z_mm"]),
                part["modeled_abs_side_norm"],
                marker=".",
                lw=1.2,
                color=color,
                label=f"{side} spectral peak",
            )
        dip = table[table["z_m"] > 0]
        ax.plot(
            np.abs(dip["z_mm"]),
            dip["dipole_abs_side_norm"],
            color="0.2",
            ls="--",
            lw=1.2,
            label="Liu dipole norm.",
        )
        ax.axvline(abs(cfg.offset_D) / math.sqrt(2.0) * 1e3, color="0.35", ls=":", lw=1.0)
        ax.set_title(f"Pe={case.pe:g}, late valid stage")
        ax.set_xlabel("distance from interface (mm)")
        ax.set_ylabel("side-normalized peak")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, loc="best")
    all_rows = pd.concat(rows, ignore_index=True)
    all_rows.to_csv(outdir / "tables" / "figure6_spatial_peak_dipole_interpretation.csv", index=False)
    fig.suptitle("Figure 6. Spatial peak, polarity, and Liu dipole interpretation")
    outpath = outdir / "figures" / "figure6_spatial_peak_dipole_interpretation.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def make_figure7(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    rows = []
    contrib_rows = []
    for case in cases:
        df = se.load_reactive_transport_table(case.input_path)
        valid = valid_rows(df, cfg)
        base = valid.iloc[0]
        target = valid.iloc[-1]
        baseline = psa.row_coefficient_metrics(base, cfg)["RE_abs"]
        for _, target_row in valid.iterrows():
            full = psa.row_coefficient_metrics(target_row, cfg)
            full.update({"case": case.name, "Pe": case.pe, "Time_s": target_row["Time_s"], "component": "full_path"})
            rows.append(full)
            for group in psa.PARAMETER_GROUPS:
                hybrid = psa.make_one_at_a_time_row(base, target_row, group)
                rec = psa.row_coefficient_metrics(hybrid, cfg)
                rec.update({"case": case.name, "Pe": case.pe, "Time_s": target_row["Time_s"], "component": group})
                rows.append(rec)
        metrics = {
            "baseline": post_t0_waveform_peak_metric(base, cfg, 64, 61),
            "full_target": post_t0_waveform_peak_metric(target, cfg, 64, 61),
        }
        for group in psa.PARAMETER_GROUPS:
            metrics[group] = post_t0_waveform_peak_metric(psa.make_one_at_a_time_row(base, target, group), cfg, 64, 61)
        contrib = psa.build_log_contribution_table(metrics, list(psa.PARAMETER_GROUPS))
        contrib["case"] = case.name
        contrib["Pe"] = case.pe
        contrib["baseline_RE_abs"] = baseline
        contrib_rows.append(contrib)
    oat = pd.DataFrame(rows)
    oat["RE_abs_norm_to_first"] = oat.groupby(["case", "component"])["RE_abs"].transform(lambda x: x / x.iloc[0])
    contrib_all = pd.concat(contrib_rows, ignore_index=True)
    oat.to_csv(outdir / "tables" / "figure7_one_at_a_time_coefficients.csv", index=False)
    contrib_all.to_csv(outdir / "tables" / "figure7_waveform_peak_contribution.csv", index=False)

    fig, axes = plt.subplots(2, len(cases), figsize=(13.0, 7.2), constrained_layout=True)
    components = ["full_path", "porosity", "permeability", "tortuosity", "fluid_chemistry"]
    colors = {
        "full_path": "black",
        "porosity": "#4c78a8",
        "permeability": "#59a14f",
        "tortuosity": "#f28e2b",
        "fluid_chemistry": "#e15759",
    }
    for col, case in enumerate(cases):
        part_case = oat[oat["case"] == case.name]
        for comp in components:
            part = part_case[part_case["component"] == comp]
            axes[0, col].plot(part["Time_s"], part["RE_abs_norm_to_first"], color=colors[comp], lw=1.2, label=comp)
        axes[0, col].set_title(f"Pe={case.pe:g}")
        axes[0, col].set_yscale("log")
        axes[0, col].set_xlabel("dissolution time (s)")
        axes[0, col].set_ylabel("|R_E| normalized")
        axes[0, col].grid(True, which="both", alpha=0.25)
        part_contrib = contrib_all[(contrib_all["case"] == case.name) & ~contrib_all["component"].eq("full_observed_change")]
        bar_colors = ["#59a14f" if v > 0 else "#e15759" if v < 0 else "0.6" for v in part_contrib["delta_log10_metric"]]
        axes[1, col].bar(part_contrib["component"], part_contrib["delta_log10_metric"], color=bar_colors)
        axes[1, col].axhline(0, color="k", lw=0.8)
        axes[1, col].tick_params(axis="x", rotation=30)
        axes[1, col].set_ylabel("delta log10 waveform peak")
        axes[1, col].grid(True, axis="y", alpha=0.25)
    axes[0, 0].legend(fontsize=7, loc="best")
    fig.suptitle("Figure 7. Mechanism decomposition and one-at-a-time sensitivity")
    outpath = outdir / "figures" / "figure7_mechanism_decomposition.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def make_figure8(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    rows = []
    f_low = 0.6 * cfg.f0
    f_high = 1.6 * cfg.f0
    for case in cases:
        df = load_case_table(case, cfg)
        for _, row in valid_rows(df, cfg).iterrows():
            low = coefficient_record_for_row(row, cfg, f_low, cfg.coeff_theta_deg)
            high = coefficient_record_for_row(row, cfg, f_high, cfg.coeff_theta_deg)
            rows.append(
                {
                    "case": case.name,
                    "Pe": case.pe,
                    "Time_s": float(row["Time_s"]),
                    "Time_norm": float(row["Time_norm"]),
                    "Porosity": float(row["Porosity"]),
                    "RE_abs_f_low": low["RE_abs"],
                    "RE_abs_f_high": high["RE_abs"],
                    "RE_spectral_ratio_high_over_low": high["RE_abs"] / low["RE_abs"] if low["RE_abs"] > 0 else np.nan,
                    "TTM_spectral_ratio_high_over_low": high["TTM_abs"] / low["TTM_abs"] if low["TTM_abs"] > 0 else np.nan,
                    "L_sigma_index": row["L_abs"] / row["sigma_abs"] if row["sigma_abs"] > 0 else np.nan,
                    "RE_abs": float(row["RE_abs"]),
                    "TTM_abs": float(row["TTM_abs"]),
                }
            )
    metrics = pd.DataFrame(rows)
    for col in ["RE_abs", "TTM_abs", "RE_spectral_ratio_high_over_low", "L_sigma_index"]:
        metrics[col + "_norm"] = metrics.groupby("case")[col].transform(lambda x: x / x.iloc[0])
    metrics.to_csv(outdir / "tables" / "figure8_normalized_monitoring_metrics.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), constrained_layout=True)
    panels = [
        ("RE_abs_norm", "normalized |R_E|"),
        ("TTM_abs_norm", "normalized |T_TM|"),
        ("RE_spectral_ratio_high_over_low_norm", "normalized |R_E(f_high)|/|R_E(f_low)|"),
        ("L_sigma_index_norm", "normalized |L|/|sigma| index"),
    ]
    for ax, (col, ylabel) in zip(axes.ravel(), panels):
        for case in cases:
            part = metrics[metrics["case"] == case.name]
            ax.plot(part["Time_norm"], part[col], color=case.color, lw=1.5, marker=".", ms=3, label=f"Pe={case.pe:g}")
        ax.set_xlabel("normalized dissolution time")
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=9)
    fig.suptitle("Figure 8. Detectability-oriented normalized monitoring metrics")
    outpath = outdir / "figures" / "figure8_normalized_monitoring_metrics.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return outpath


def write_results_summary(cases: Sequence[DissolutionCase], cfg: se.SEConfig, outdir: Path) -> Path:
    summary = []
    for case in cases:
        df = load_case_table(case, cfg)
        valid = valid_rows(df, cfg)
        first = valid.iloc[0]
        last = valid.iloc[-1]
        summary.append(
            {
                "case": case.name,
                "Pe": case.pe,
                "n_rows": len(df),
                "valid_rows": len(valid),
                "final_Time_s": float(df["Time_s"].iloc[-1]),
                "valid_final_Time_s": float(last["Time_s"]),
                "porosity_initial": float(first["Porosity"]),
                "porosity_valid_final": float(last["Porosity"]),
                "permeability_factor_valid": float(last["Permeability_mD"] / first["Permeability_mD"]),
                "tortuosity_change_valid": float(last["Tortuosity"] - first["Tortuosity"]),
                "C_molL_factor_valid": float(last["C_molL"] / first["C_molL"]) if first["C_molL"] > 0 else np.nan,
                "RE_abs_factor_valid": float(last["RE_abs"] / first["RE_abs"]) if first["RE_abs"] > 0 else np.nan,
                "TTM_abs_factor_valid": float(last["TTM_abs"] / first["TTM_abs"]) if first["TTM_abs"] > 0 else np.nan,
                "L_abs_factor_valid": float(last["L_abs"] / first["L_abs"]) if first["L_abs"] > 0 else np.nan,
                "sigma_abs_factor_valid": float(last["sigma_abs"] / first["sigma_abs"]) if first["sigma_abs"] > 0 else np.nan,
            }
        )
    out = pd.DataFrame(summary)
    path = outdir / "tables" / "cross_case_key_findings.csv"
    out.to_csv(path, index=False)
    return path


def run_all(outdir: Path = DEFAULT_OUTDIR, cases: Sequence[DissolutionCase] = DEFAULT_CASES) -> List[Path]:
    cfg = default_cfg()
    ensure_outdir(outdir)
    paths = [
        make_figure1(outdir),
        make_figure2(cases, cfg, outdir),
        make_figure3(cases, cfg, outdir),
        make_figure4(cases, cfg, outdir),
        make_figure5(cases, cfg, outdir),
        make_figure6(cases, cfg, outdir),
        make_figure7(cases, cfg, outdir),
        make_figure8(cases, cfg, outdir),
    ]
    write_results_summary(cases, cfg, outdir)
    return paths


def parse_outdir() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    return parser.parse_args()
