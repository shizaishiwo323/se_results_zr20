import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "schakel2010_strict_sensitivity.py"
spec = importlib.util.spec_from_file_location("schakel2010_strict_sensitivity", MODULE_PATH)
strict = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = strict
spec.loader.exec_module(strict)


class Schakel2010StrictFluxTest(unittest.TestCase):
    def test_reproduction_script_is_independent_from_liu2018_forward_model(self):
        source = MODULE_PATH.read_text(encoding="utf-8")

        self.assertNotIn("seismoelectric_offset_liu2018_spectral", source)
        self.assertNotIn("parameter_sensitivity_analysis", source)

    def test_table_ii_energy_flux_coefficients_at_theta_30_for_10hz_and_500khz(self):
        table = strict.table_ii_validation_rows()
        by_frequency = {row["frequency_hz"]: row for row in table}

        low = by_frequency[10.0]
        high = by_frequency[500_000.0]

        self.assertTrue(np.isclose(low["RE_EE"], 7.4497e-11, rtol=3.0e-3))
        self.assertTrue(np.isclose(low["TE_TM_TM"], -2.6500e-10, rtol=3.0e-3))
        self.assertTrue(np.isclose(low["RE_Pr_Pr"], -4.0035e-1, rtol=3.0e-3))
        self.assertTrue(np.isclose(low["energy_balance_residual"], 0.0, atol=5.0e-5))

        self.assertTrue(np.isclose(high["RE_EE"], 1.1871e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(high["TE_TM_TM"], -4.9895e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(high["RE_Pr_Pr"], -3.2615e-1, rtol=3.0e-3))
        self.assertTrue(np.isclose(high["energy_balance_residual"], 0.0, atol=5.0e-5))
        self.assertEqual(high["TE_Pf_Pf_source"], "paper_table")
        self.assertTrue(np.isnan(high["TE_Pf_Pf_relative_error"]))
        self.assertEqual(high["TE_TM_TM_source"], "computed")

    def test_core_flux_api_does_not_expose_unverified_interference_terms(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=45.0)

        self.assertIn("TE_TM_TM", result)
        self.assertNotIn("TE_Pf_Pf", result)
        self.assertNotIn("TE", result)
        self.assertNotIn("energy_balance", result)

    def test_fig2_angle_scan_uses_schakel_potential_coefficients(self):
        df = strict.figure2_angle_scan(n=81)

        self.assertIn("sin_theta", df.columns)
        self.assertIn("R_E_abs", df.columns)
        self.assertIn("T_TM_abs", df.columns)
        self.assertIn("R_E_phase_rad", df.columns)
        self.assertIn("T_TM_phase_rad", df.columns)
        self.assertGreaterEqual(df["T_TM_phase_rad"].min(), 0.0)
        self.assertLessEqual(df["T_TM_phase_rad"].max(), 4.0 * np.pi + 1.0e-9)
        self.assertTrue(np.isclose(df["sin_theta"].iloc[0], 0.0))
        self.assertTrue(np.isclose(df["sin_theta"].iloc[-1], 2.0))
        self.assertLess(df.loc[df["sin_theta"].sub(0.0).abs().idxmin(), "R_E_abs"], 1.0e-6)
        self.assertLess(df.loc[df["sin_theta"].sub(1.0).abs().idxmin(), "R_E_abs"], 1.0e-6)
        self.assertGreater(df.loc[df["sin_theta"].sub(0.5).abs().idxmin(), "R_E_abs"], 1.0)
        self.assertGreater(df["T_TM_abs"].max(), 5.0e-10)
        self.assertLess(df["T_TM_abs"].max(), 1.0e-7)

    def test_table_iii_energy_flux_coefficients_at_omega_1e6_theta_45(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=45.0)

        self.assertTrue(np.isclose(result["RE_EE"], 5.4312e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["TE_TM_TM"], -1.6581e-6, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["RE_Pr_Pr"], -5.6982e-2, rtol=3.0e-3))

    def test_table_iii_energy_flux_coefficients_at_omega_1e6_theta_30(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=30.0)

        self.assertTrue(np.isclose(result["RE_EE"], 1.3503e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["TE_TM_TM"], -5.0880e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["RE_Pr_Pr"], -3.2877e-1, rtol=3.0e-3))

    def test_strict_coefficients_are_not_squared_potential_proxies(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=45.0)

        self.assertGreater(abs(result["R_E"]) ** 2, 1.0e16)
        self.assertLess(result["RE_EE"], 1.0e-5)
        self.assertLess(abs(result["TE_TM_TM"]), 1.0e-4)


if __name__ == "__main__":
    unittest.main()
