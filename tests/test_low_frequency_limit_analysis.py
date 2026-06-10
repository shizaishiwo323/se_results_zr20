import importlib.util
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "analysis_scripts" / "low_frequency_limit_analysis.py"
spec = importlib.util.spec_from_file_location("low_frequency_limit_analysis", MODULE_PATH)
lf = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = lf
spec.loader.exec_module(lf)


class LowFrequencyLimitAnalysisTest(unittest.TestCase):
    def sample_row(self):
        return pd.Series(
            {
                "Time_s": 0.0,
                "Porosity": 0.25,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

    def test_darcy_limit_matches_near_zero_dynamic_terms(self):
        cfg = lf.se.SEConfig()
        row = self.sample_row()
        params = lf.row_transport_params(row, cfg)

        limit = lf.darcy_dynamic_limit(**params, cfg=cfg)
        near = lf.se.dynamic_coefficients(
            params["phi"],
            params["k0_m2"],
            params["alpha_inf"],
            params["cH"],
            2.0 * math.pi * 1.0e-9,
            cfg,
        )

        self.assertTrue(np.isclose(limit["k_dyn"].real, params["k0_m2"]))
        self.assertTrue(np.isclose(limit["L"].real, near["L"].real, rtol=1.0e-8))
        self.assertTrue(np.isclose(limit["sigma"].real, near["sigma"].real, rtol=1.0e-8))
        self.assertTrue(np.isfinite(limit["darcy_response_index"]))

    def test_zero_frequency_record_is_darcy_diagnostic_not_wave_solver(self):
        cfg = lf.se.SEConfig()
        records = lf.low_frequency_records_for_row(
            self.sample_row(),
            cfg,
            frequencies_hz=[0.0, 1.0e-3],
            theta_deg=45.0,
        )

        zero = records.loc[records["frequency_Hz"] == 0.0].iloc[0]
        near = records.loc[records["frequency_Hz"] == 1.0e-3].iloc[0]

        self.assertEqual(zero["regime"], "Darcy limit")
        self.assertEqual(zero["interface_solver_status"], "not evaluated at singular omega=0")
        self.assertTrue(np.isnan(zero["RE_abs"]))
        self.assertTrue(np.isfinite(zero["darcy_response_index_abs"]))
        self.assertEqual(near["regime"], "near-zero frequency")
        self.assertTrue(np.isfinite(near["RE_abs"]))


if __name__ == "__main__":
    unittest.main()
