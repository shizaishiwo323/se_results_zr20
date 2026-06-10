#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnose finite-window pre-T0 energy in the waveform panels.

This script does not alter waveform outputs.  It summarizes the ratio between
the largest absolute amplitude before the acoustic interface-arrival time T0
and the largest absolute amplitude after T0 for the cached Figure 5 panels.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from recommended_figures_common import ensure_outdir, parse_outdir


def main() -> None:
    args = parse_outdir()
    outdir = ensure_outdir(Path(args.outdir))
    summary_path = outdir / "tables" / "figure5_waveform_panel_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"{summary_path} does not exist. Run figure5_waveform_panels.py first."
        )
    df = pd.read_csv(summary_path)
    df["pre_to_post_T0_ratio"] = df["pre_T0_max_abs"] / df["post_T0_max_abs"]
    diag_path = outdir / "tables" / "waveform_t0_leakage_diagnostics.csv"
    df.to_csv(diag_path, index=False)

    labels = [f"{row.case}\nstep {int(row.TimeStep)}" for row in df.itertuples()]
    fig, ax = plt.subplots(figsize=(10.5, 4.2), constrained_layout=True)
    ax.bar(labels, df["pre_to_post_T0_ratio"], color="#8a6f3f")
    ax.axhline(0.1, color="0.2", ls=":", lw=1.1, label=r"$10\%$ reference")
    ax.axhline(1.0, color="#b05a32", ls="--", lw=1.1, label=r"pre-$T_0$ = post-$T_0$")
    ax.set_ylabel(r"$\max |u|_{\mathrm{pre}\!-\!T_0}/\max |u|_{\mathrm{post}\!-\!T_0}$")
    ax.set_title(r"Waveform $T_0$ leakage diagnostic for finite-band spectral panels")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    outpath = outdir / "figures" / "waveform_t0_leakage_diagnostics.png"
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    print(diag_path)
    print(outpath)


if __name__ == "__main__":
    main()
