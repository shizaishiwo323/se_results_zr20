import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "analysis_scripts" / "recommended_figures_common.py"
spec = importlib.util.spec_from_file_location("recommended_figures_common", MODULE_PATH)
rfg = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = rfg
spec.loader.exec_module(rfg)


class RecommendedFigureSequenceTest(unittest.TestCase):
    def test_default_cases_point_to_existing_reactive_transport_inputs(self):
        for case in rfg.DEFAULT_CASES:
            self.assertTrue(case.input_path.exists(), case.input_path)
            self.assertTrue(case.tortuosity_path.exists(), case.tortuosity_path)

    def test_representative_indices_are_spread_across_valid_poroelastic_window(self):
        cfg = rfg.default_cfg()
        case = [c for c in rfg.DEFAULT_CASES if c.name == "Pe10"][0]
        df = rfg.load_case_table(case, cfg)

        idxs = rfg.representative_indices(df, cfg, n=3)

        self.assertEqual(len(idxs), 3)
        times = df.loc[idxs, "Time_s"].to_numpy(float)
        self.assertTrue(np.all(np.diff(times) > 0))
        self.assertLess(df.loc[idxs[-1], "Porosity"], cfg.phi_max_valid)

    def test_dynamic_frequency_scan_preserves_stage_identity(self):
        cfg = rfg.default_cfg()
        case = rfg.DEFAULT_CASES[0]
        df = rfg.load_case_table(case, cfg)
        row = df.loc[rfg.representative_indices(df, cfg, n=2)[0]]
        freqs = np.array([0.5 * cfg.f0, cfg.f0])

        scan = rfg.row_dynamic_frequency_scan(row, cfg, freqs)

        self.assertEqual(list(scan["frequency_Hz"]), list(freqs))
        self.assertTrue((scan["Time_s"] == row["Time_s"]).all())
        self.assertTrue((scan[["k_dyn_abs", "L_abs", "sigma_abs"]] > 0).all().all())


if __name__ == "__main__":
    unittest.main()

