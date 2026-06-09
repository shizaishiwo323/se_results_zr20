#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path

from recommended_figures_common import DEFAULT_CASES, default_cfg, ensure_outdir, make_figure8, parse_outdir


def main() -> None:
    args = parse_outdir()
    outdir = ensure_outdir(Path(args.outdir))
    print(make_figure8(DEFAULT_CASES, default_cfg(), outdir))


if __name__ == "__main__":
    main()

