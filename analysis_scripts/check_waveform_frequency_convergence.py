#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Frequency-sampling convergence check for the finite-offset waveform panels."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from recommended_figures_common import DEFAULT_CASES, default_cfg, ensure_outdir, load_case_table, parse_outdir
import seismoelectric_offset_liu2018_spectral as se


def main() -> None:
    args = parse_outdir()
    outdir = ensure_outdir(Path(args.outdir))
    cfg = default_cfg()
    case = DEFAULT_CASES[0]
    row = load_case_table(case, cfg).iloc[0]
    rows = []
    for n_omega in [16, 32, 64, 128, 256]:
        z, t, U = se.synthesize_waveforms_spectral(row, cfg, n_omega=n_omega, n_k=cfg.spectral_n_k)
        del z
        T0 = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)
        pre = float(abs(U[:, t < T0]).max())
        post = float(abs(U[:, t >= T0]).max())
        df_hz = (cfg.spectral_f_max_factor - cfg.spectral_f_min_factor) * cfg.f0 / (n_omega - 1)
        rows.append(
            {
                "case": case.name,
                "TimeStep": int(row["TimeStep"]),
                "n_omega": n_omega,
                "n_k": cfg.spectral_n_k,
                "frequency_spacing_Hz": df_hz,
                "time_alias_period_us": 1.0 / df_hz * 1e6,
                "post_T0_peak": post,
                "post_T0_peak_norm_to_n256": None,
                "pre_to_post_T0_ratio": pre / post,
            }
        )
    df = pd.DataFrame(rows)
    ref = float(df.loc[df["n_omega"] == 256, "post_T0_peak"].iloc[0])
    df["post_T0_peak_norm_to_n256"] = df["post_T0_peak"] / ref
    table_path = outdir / "tables" / "waveform_frequency_convergence.csv"
    df.to_csv(table_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), constrained_layout=True)
    axes[0].plot(df["n_omega"], df["post_T0_peak_norm_to_n256"], marker="o")
    axes[0].axhline(1.0, color="0.3", ls=":", lw=1.0)
    axes[0].set_xlabel("n_omega")
    axes[0].set_ylabel("post-T0 peak / n=256")
    axes[0].grid(True, alpha=0.25)
    axes[1].plot(df["n_omega"], df["pre_to_post_T0_ratio"], marker="o", color="#8a6f3f")
    axes[1].axhline(0.1, color="0.3", ls=":", lw=1.0)
    axes[1].set_xlabel("n_omega")
    axes[1].set_ylabel("pre/post T0 ratio")
    axes[1].set_yscale("log")
    axes[1].grid(True, which="both", alpha=0.25)
    outpath = outdir / "figures" / "waveform_frequency_convergence.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    print(table_path)
    print(outpath)


if __name__ == "__main__":
    main()

