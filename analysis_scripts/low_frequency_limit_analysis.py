#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Low-frequency / zero-frequency diagnostic for the seismoelectric model.

The Schakel/Pride frequency-domain wave solver contains terms proportional to
1/omega and 1/omega**2.  Therefore the exact omega=0 point is not evaluated by
the wave boundary-value solver.  Instead, this script reports the analytic
Darcy-flow limit of the dynamic electrokinetic properties at f=0 and compares it
with near-zero positive-frequency interface solutions.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SE_PATH = ROOT / "seismoelectric_offset_liu2018_spectral.py"
spec = importlib.util.spec_from_file_location("seismoelectric_offset_liu2018_spectral", SE_PATH)
se = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = se
spec.loader.exec_module(se)


DEFAULT_FREQUENCIES_HZ = [
    0.0,
    1.0e-6,
    1.0e-5,
    1.0e-4,
    1.0e-3,
    1.0e-2,
    1.0e-1,
    1.0,
    10.0,
    100.0,
    1000.0,
]


def row_transport_params(row: pd.Series, cfg: se.SEConfig) -> Dict[str, float | None]:
    phi = float(np.clip(float(row["Porosity"]), cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
    alpha_inf = max(float(row["Tortuosity"]), 1.0 + 1.0e-6)
    return {
        "phi": phi,
        "k0_m2": k0_m2,
        "alpha_inf": alpha_inf,
        "cH": float(row["OutletHConc"]),
        "C_override_molL": se.optional_float(row, "ElectrolyteConcentration_molL"),
        "sigma_f_override": se.optional_float(row, "FluidConductivity_S_m"),
    }


def darcy_dynamic_limit(
    phi: float,
    k0_m2: float,
    alpha_inf: float,
    cH: float,
    cfg: se.SEConfig,
    C_override_molL: float | None = None,
    sigma_f_override: float | None = None,
) -> Dict[str, complex | float]:
    """Analytic omega -> 0 limit of k(omega), L(omega), and sigma(omega)."""
    phi = float(np.clip(phi, cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(k0_m2), cfg.k0_min)
    alpha_inf = max(float(alpha_inf), 1.0 + 1.0e-6)
    ec = se.electrochemistry_from_h(cH, cfg, C_override_molL=C_override_molL)

    z_vals = np.array([cfg.z1, cfg.z2], dtype=float)
    b_vals = np.array([cfg.b1, cfg.b2], dtype=float)
    n_each = ec["C_molL"] * 1000.0 * cfg.N_A
    n_vals = np.array([n_each, n_each], dtype=float)

    omega_t = phi * cfg.eta / (alpha_inf * k0_m2 * cfg.rho_f)
    lambda_m = math.sqrt(max(8.0 * alpha_inf * k0_m2 / (phi * cfg.M_similarity), cfg.eps_complex))

    denom_debye = np.sum((cfg.e_charge * z_vals) ** 2 * n_vals) / (
        cfg.eps0 * cfg.eps_f * cfg.k_B * cfg.temperature
    )
    debye_d = math.sqrt(1.0 / max(float(np.real(denom_debye)), cfg.eps_complex))
    one_minus = 1.0 - 2.0 * debye_d / lambda_m

    k_dyn = complex(k0_m2, 0.0)
    L = -((phi / alpha_inf) * (cfg.eps0 * cfg.eps_f * ec["zeta"] / cfg.eta) * one_minus)

    sigma_f = float(np.sum((cfg.e_charge * z_vals) ** 2 * b_vals * n_vals))
    if sigma_f_override is not None and np.isfinite(sigma_f_override) and sigma_f_override > 0:
        sigma_f = float(sigma_f_override)

    if abs(ec["zeta"]) < 1.0e-12 or sigma_f <= 0.0:
        c_em = 0.0
        c_os = 0.0
        p_os = 0.0
    else:
        exp_terms = np.exp(-(cfg.e_charge * z_vals * ec["zeta"]) / (2.0 * cfg.k_B * cfg.temperature)) - 1.0
        c_em = float(2.0 * debye_d * np.sum((cfg.e_charge * z_vals) ** 2 * b_vals * n_vals * exp_terms))
        p_os = float(
            (8.0 * cfg.k_B * cfg.temperature * debye_d**2)
            / (cfg.eps0 * cfg.eps_f * ec["zeta"] ** 2)
            * np.sum(n_vals * exp_terms)
        )
        if abs(p_os) < 1.0e-30:
            c_os = 0.0
        else:
            c_os = ((cfg.eps0 * cfg.eps_f) ** 2 * ec["zeta"] ** 2 / (2.0 * debye_d * cfg.eta)) * p_os

    sigma = (phi * sigma_f / alpha_inf) * (
        1.0 + 2.0 * (c_em + c_os) / max(sigma_f * lambda_m, cfg.eps_complex)
    )
    darcy_response_index = complex(L, 0.0) / complex(sigma, 0.0) if sigma != 0.0 else np.nan + 0j

    return {
        **ec,
        "omega_t": omega_t,
        "Lambda": lambda_m,
        "debye_d": debye_d,
        "k_dyn": k_dyn,
        "L": complex(L, 0.0),
        "sigma_f": sigma_f,
        "C_em": c_em,
        "C_os": complex(c_os, 0.0),
        "sigma": complex(sigma, 0.0),
        "darcy_response_index": darcy_response_index,
    }


def _base_record(row: pd.Series, cfg: se.SEConfig, params: Dict[str, float | None]) -> Dict[str, float]:
    return {
        "Time_s": float(row["Time_s"]),
        "Time_min": float(row["Time_s"]) / 60.0,
        "Porosity_raw": float(row["Porosity"]),
        "Porosity_used": float(params["phi"]),
        "Permeability_mD": float(row["Permeability_mD"]),
        "k0_m2": float(params["k0_m2"]),
        "Tortuosity": float(params["alpha_inf"]),
        "OutletHConc_raw": float(row["OutletHConc"]),
    }


def _dynamic_record_fields(dyn: Dict[str, complex | float]) -> Dict[str, float]:
    idx = dyn["darcy_response_index"] if "darcy_response_index" in dyn else dyn["L"] / dyn["sigma"]
    return {
        "cH_molL": float(dyn["cH_molL"]),
        "pH": float(dyn["pH"]),
        "C_molL": float(dyn["C_molL"]),
        "zeta_V": float(dyn["zeta"]),
        "omega_t_rad_s": float(dyn["omega_t"]),
        "Lambda_m": float(dyn["Lambda"]),
        "debye_d_m": float(dyn["debye_d"]),
        "k_dyn_real": float(np.real(dyn["k_dyn"])),
        "k_dyn_imag": float(np.imag(dyn["k_dyn"])),
        "k_dyn_abs": float(abs(dyn["k_dyn"])),
        "L_real": float(np.real(dyn["L"])),
        "L_imag": float(np.imag(dyn["L"])),
        "L_abs": float(abs(dyn["L"])),
        "sigma_f_S_m": float(dyn["sigma_f"]),
        "sigma_real": float(np.real(dyn["sigma"])),
        "sigma_imag": float(np.imag(dyn["sigma"])),
        "sigma_abs": float(abs(dyn["sigma"])),
        "darcy_response_index_real": float(np.real(idx)),
        "darcy_response_index_imag": float(np.imag(idx)),
        "darcy_response_index_abs": float(abs(idx)),
    }


def low_frequency_records_for_row(
    row: pd.Series,
    cfg: se.SEConfig,
    frequencies_hz: Iterable[float],
    theta_deg: float,
) -> pd.DataFrame:
    params = row_transport_params(row, cfg)
    records: List[Dict[str, float | str]] = []
    base = _base_record(row, cfg, params)
    theta_rad = math.radians(theta_deg)
    vf = math.sqrt(cfg.K_fl / cfg.rho_fl)

    for freq_hz in frequencies_hz:
        freq_hz = float(freq_hz)
        if freq_hz == 0.0:
            dyn0 = darcy_dynamic_limit(**params, cfg=cfg)
            records.append(
                {
                    **base,
                    **_dynamic_record_fields(dyn0),
                    "frequency_Hz": 0.0,
                    "omega_rad_s": 0.0,
                    "omega_over_omega_t": 0.0,
                    "theta_deg": float(theta_deg),
                    "kx_m_inv": 0.0,
                    "regime": "Darcy limit",
                    "interface_solver_status": "not evaluated at singular omega=0",
                    "RE_real": np.nan,
                    "RE_imag": np.nan,
                    "RE_abs": np.nan,
                    "TTM_real": np.nan,
                    "TTM_imag": np.nan,
                    "TTM_abs": np.nan,
                    "Ru_abs": np.nan,
                    "Tu_abs": np.nan,
                    "matrix_cond": np.nan,
                }
            )
            continue

        omega = 2.0 * math.pi * freq_hz
        dyn = se.dynamic_coefficients(
            params["phi"],
            params["k0_m2"],
            params["alpha_inf"],
            params["cH"],
            omega,
            cfg,
            C_override_molL=params["C_override_molL"],
            sigma_f_override=params["sigma_f_override"],
        )
        rec: Dict[str, float | str] = {
            **base,
            **_dynamic_record_fields(dyn),
            "frequency_Hz": freq_hz,
            "omega_rad_s": omega,
            "omega_over_omega_t": omega / float(dyn["omega_t"]),
            "theta_deg": float(theta_deg),
            "kx_m_inv": omega / vf * math.sin(theta_rad),
            "regime": "near-zero frequency",
        }
        try:
            coeff = se.se_coefficients(
                params["phi"],
                params["k0_m2"],
                params["alpha_inf"],
                params["cH"],
                omega,
                theta_deg,
                cfg,
                C_override_molL=params["C_override_molL"],
                sigma_f_override=params["sigma_f_override"],
            )
            ru, tu = se.liu_electrical_potential_coefficients(coeff, rec["kx_m_inv"])
            rec.update(
                {
                    "interface_solver_status": "ok",
                    "RE_real": float(np.real(coeff["R_E"])),
                    "RE_imag": float(np.imag(coeff["R_E"])),
                    "RE_abs": float(abs(coeff["R_E"])),
                    "TTM_real": float(np.real(coeff["T_TM"])),
                    "TTM_imag": float(np.imag(coeff["T_TM"])),
                    "TTM_abs": float(abs(coeff["T_TM"])),
                    "Ru_abs": float(abs(ru)),
                    "Tu_abs": float(abs(tu)),
                    "matrix_cond": float(coeff["matrix_cond"]),
                }
            )
        except Exception as exc:
            rec.update(
                {
                    "interface_solver_status": f"failed: {type(exc).__name__}",
                    "RE_real": np.nan,
                    "RE_imag": np.nan,
                    "RE_abs": np.nan,
                    "TTM_real": np.nan,
                    "TTM_imag": np.nan,
                    "TTM_abs": np.nan,
                    "Ru_abs": np.nan,
                    "Tu_abs": np.nan,
                    "matrix_cond": np.nan,
                }
            )
        records.append(rec)

    return pd.DataFrame(records)


def build_low_frequency_table(
    df: pd.DataFrame,
    cfg: se.SEConfig,
    frequencies_hz: Iterable[float],
    theta_deg: float,
) -> pd.DataFrame:
    parts = [low_frequency_records_for_row(row, cfg, frequencies_hz, theta_deg) for _, row in df.iterrows()]
    out = pd.concat(parts, ignore_index=True)
    first_idx = out["frequency_Hz"] == 0.0
    ref = out.loc[first_idx, "darcy_response_index_abs"].replace(0.0, np.nan).iloc[0]
    out["darcy_response_index_abs_norm_to_initial"] = out["darcy_response_index_abs"] / ref
    return out


def representative_times(df: pd.DataFrame) -> List[float]:
    if len(df) <= 3:
        return [float(v) for v in df["Time_s"]]
    idx = sorted({0, len(df) // 2, len(df) - 1})
    return [float(df.iloc[i]["Time_s"]) for i in idx]


def plot_darcy_evolution(darcy: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True)
    x = darcy["Time_min"].to_numpy(float)
    axes[0].plot(x, darcy["darcy_response_index_abs_norm_to_initial"], marker="o", lw=1.5)
    axes[0].set_ylabel(r"$|L(0)/\sigma(0)|$ norm.")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(x, darcy["L_real"], marker="o", lw=1.5, label=r"$L(0)$")
    axes[1].plot(x, darcy["sigma_real"], marker="s", lw=1.5, label=r"$\sigma(0)$")
    axes[1].set_yscale("symlog", linthresh=1.0e-12)
    axes[1].set_xlabel("Dissolution time (min)")
    axes[1].set_ylabel("Darcy-limit value")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.suptitle("Zero-frequency Darcy electrokinetic response")
    fig.tight_layout()
    fig.savefig(outdir / "darcy_response_index_vs_dissolution_time.png", dpi=300)
    plt.close(fig)


def plot_frequency_convergence(table: pd.DataFrame, outdir: Path) -> None:
    times = representative_times(table.loc[table["frequency_Hz"] == 0.0])
    finite = table[table["frequency_Hz"] > 0.0].copy()
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True)
    metrics = [
        ("k_dyn_abs", r"$|k(\omega)|$"),
        ("L_abs", r"$|L(\omega)|$"),
        ("sigma_abs", r"$|\sigma(\omega)|$"),
        ("darcy_response_index_abs", r"$|L(\omega)/\sigma(\omega)|$"),
    ]
    for ax, (col, label) in zip(axes.ravel(), metrics):
        y_all = []
        for time_s in times:
            part = finite[np.isclose(finite["Time_s"], time_s)]
            zero = table[(np.isclose(table["Time_s"], time_s)) & (table["frequency_Hz"] == 0.0)][col].iloc[0]
            y = part[col].to_numpy(float) / zero if zero != 0.0 else part[col].to_numpy(float)
            y_all.extend([float(v) for v in y if np.isfinite(v)])
            ax.semilogx(part["frequency_Hz"], y, marker="o", lw=1.2, label=f"{time_s:g} s")
        ax.set_xscale("log")
        if y_all and min(y_all) >= 0.0:
            ax.set_ylim(0.0, max(1.05, max(y_all) * 1.05))
        ax.set_ylabel(f"{label} / 0 Hz")
        ax.grid(True, alpha=0.3)
    axes[-1, 0].set_xlabel("Frequency (Hz)")
    axes[-1, 1].set_xlabel("Frequency (Hz)")
    axes[0, 0].legend(title="Dissolution time")
    fig.suptitle("Near-zero dynamic-property convergence to Darcy limit")
    fig.tight_layout()
    fig.savefig(outdir / "low_frequency_dynamic_convergence.png", dpi=300)
    plt.close(fig)


def plot_interface_diagnostics(table: pd.DataFrame, outdir: Path) -> None:
    times = representative_times(table.loc[table["frequency_Hz"] == 0.0])
    finite = table[table["frequency_Hz"] > 0.0].copy()
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True)
    for ax, col, label in [
        (axes[0, 0], "RE_abs", r"$|R_E|$"),
        (axes[0, 1], "TTM_abs", r"$|T_{TM}|$"),
        (axes[1, 0], "Ru_abs", r"$|R_u|$"),
        (axes[1, 1], "matrix_cond", "matrix condition number"),
    ]:
        for time_s in times:
            part = finite[np.isclose(finite["Time_s"], time_s)]
            ax.plot(part["frequency_Hz"], part[col], marker="o", lw=1.2, label=f"{time_s:g} s")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
    axes[1, 0].set_xlabel("Frequency (Hz)")
    axes[1, 1].set_xlabel("Frequency (Hz)")
    axes[0, 0].legend(title="Dissolution time")
    fig.suptitle("Near-zero interface solver diagnostics")
    fig.tight_layout()
    fig.savefig(outdir / "near_zero_interface_coefficients.png", dpi=300)
    plt.close(fig)


def write_summary(table: pd.DataFrame, outdir: Path, frequencies_hz: Iterable[float], theta_deg: float) -> None:
    darcy = table[table["frequency_Hz"] == 0.0].copy()
    finite = table[table["frequency_Hz"] > 0.0].copy()
    ratio0 = float(darcy["darcy_response_index_abs_norm_to_initial"].iloc[-1])
    min_cond = float(np.nanmin(finite["matrix_cond"]))
    max_cond = float(np.nanmax(finite["matrix_cond"]))
    text = f"""# Low-Frequency Limit Diagnostic

This directory evaluates the model behavior at the zero-frequency / Darcy-flow
limit and compares it with near-zero positive-frequency Schakel interface
solutions.

## Main interpretation

- The exact `frequency_Hz = 0` rows are analytic Darcy-limit diagnostics. They
  report `k(0) = k0`, `L(0)`, `sigma(0)`, and the steady electrokinetic response
  index `L(0)/sigma(0)`.
- The full Schakel/Pride wave boundary-value problem is not evaluated at
  `omega = 0`, because the equations contain conductive terms proportional to
  `1/omega` and poro-electromagnetic terms proportional to `1/omega^2`.
- Near-zero positive frequencies are retained as convergence diagnostics for
  the interface solver at theta = {theta_deg:g} degrees.
- Across the dissolution sequence, the final normalized Darcy response index is
  {ratio0:.6g} relative to the initial value.
- The near-zero interface matrix condition number ranges from {min_cond:.3e} to
  {max_cond:.3e}, so the interface coefficients near zero frequency should be
  interpreted as numerical-limit diagnostics rather than a time-domain waveform.

## Files

- `low_frequency_dynamic_and_interface.csv`: all frequencies and dissolution
  times, including the `frequency_Hz = 0` Darcy-limit rows.
- `darcy_limit_vs_dissolution_time.csv`: only the exact zero-frequency Darcy
  rows.
- `low_frequency_run_summary.csv`: run settings and output overview.
- `darcy_response_index_vs_dissolution_time.png`: evolution of the steady
  Darcy electrokinetic response index.
- `low_frequency_dynamic_convergence.png`: near-zero convergence of dynamic
  properties to their Darcy limits.
- `near_zero_interface_coefficients.png`: interface coefficient and matrix
  conditioning diagnostics for positive low frequencies.

The sampled frequencies are: {', '.join(f'{float(f):g}' for f in frequencies_hz)} Hz.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def run_analysis(input_path: Path, outdir: Path, frequencies_hz: Iterable[float], theta_deg: float) -> pd.DataFrame:
    cfg = se.SEConfig()
    df = se.load_reactive_transport_table(input_path)
    outdir.mkdir(parents=True, exist_ok=True)

    table = build_low_frequency_table(df, cfg, frequencies_hz, theta_deg)
    darcy = table[table["frequency_Hz"] == 0.0].copy()

    table.to_csv(outdir / "low_frequency_dynamic_and_interface.csv", index=False)
    darcy.to_csv(outdir / "darcy_limit_vs_dissolution_time.csv", index=False)
    pd.DataFrame(
        [
            {
                "input": str(input_path),
                "n_dissolution_steps": len(df),
                "n_frequency_samples": len(list(frequencies_hz)),
                "theta_deg": float(theta_deg),
                "frequency_min_positive_Hz": min(float(f) for f in frequencies_hz if float(f) > 0.0),
                "frequency_max_Hz": max(float(f) for f in frequencies_hz),
                "output_directory": str(outdir),
            }
        ]
    ).to_csv(outdir / "low_frequency_run_summary.csv", index=False)

    plot_darcy_evolution(darcy, outdir)
    plot_frequency_convergence(table, outdir)
    plot_interface_diagnostics(table, outdir)
    write_summary(table, outdir, frequencies_hz, theta_deg)
    return table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate the low-frequency / Darcy limit of the SE model.")
    parser.add_argument("--input", type=Path, default=ROOT / "global_evolution.xlsx")
    parser.add_argument("--outdir", type=Path, default=ROOT / "low_frequency_limit_results")
    parser.add_argument("--theta-deg", type=float, default=45.0)
    parser.add_argument(
        "--frequencies-hz",
        type=float,
        nargs="+",
        default=DEFAULT_FREQUENCIES_HZ,
        help="Frequencies to sample. Include 0 for the analytic Darcy limit.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    freqs = sorted(set(float(f) for f in args.frequencies_hz))
    if 0.0 not in freqs:
        freqs = [0.0] + freqs
    run_analysis(args.input, args.outdir, freqs, args.theta_deg)


if __name__ == "__main__":
    main()
