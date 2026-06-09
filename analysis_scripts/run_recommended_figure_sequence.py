#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path

from recommended_figures_common import DEFAULT_CASES, parse_outdir, run_all


def main() -> None:
    args = parse_outdir()
    paths = run_all(Path(args.outdir), DEFAULT_CASES)
    print("Generated recommended figure sequence:")
    for path in paths:
        print(f" - {path}")


if __name__ == "__main__":
    main()

