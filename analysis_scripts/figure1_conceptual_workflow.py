#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path

from recommended_figures_common import make_figure1, parse_outdir


def main() -> None:
    args = parse_outdir()
    print(make_figure1(Path(args.outdir)))


if __name__ == "__main__":
    main()

