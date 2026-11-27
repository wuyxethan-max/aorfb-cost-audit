"""Unit and benchmark regression tests for aorfb_audit.

Compatible with both `pytest` and Python's built-in `unittest` runner.
"""

from math import isclose
import unittest

import numpy as np

from aorfb_audit import (
    FARADAY_C_PER_MOL,
    J_PER_KWH,
    YANG_BASE_UD,
    YANG_PUBLISHED_TARGETS,
    YANG_STAGE_INPUTS,
    HalfCell,
    benchmark_reconciliation,
    capital_cost_grid,
    eq1_capital,
    eq5_electrolyte,
    eq6_maintenance,
    evaluate_molecule_cost,
    load_emmel_dataset,
    load_yang_dataset,
    max_fade_for_cost,
    npv_electrolyte,
    polymorph_sensitivity_summary,
    protocol_sensitivity_summary,
    validate_datasets,
    yang_half_cell,
)


class TestAORFBCostAudit(unittest.TestCase):
    """Verify physical identities, dataset integrity, and Phase 3 audit benchmarks."""

    def test_physical_identities_and_edge_cases(self):
        self.assertTrue(isclose(FARADAY_C_PER_MOL / J_PER_KWH, 0.026801480, rel_tol=1e-7))
        self.assertTrue(isclose(eq1_capital(20, 800, 50), 270.0))

        simple_half = HalfCell(100, 1, 1, 1, 1, 0, 0, 0, solvent_usd_kg=0, soc=1)
        one_dollar_per_10_mol = eq5_electrolyte((simple_half,), voltage=1.0, tank_usd_m3=0)
        self.assertTrue(
            isclose(
                one_dollar_per_10_mol["active"],
                0.1 / ((FARADAY_C_PER_MOL / J_PER_KWH) * 0.94 * 0.97),
            )
        )
        self.assertTrue(
            isclose(
                sum(eq5_electrolyte((simple_half,), 2.0, 0).values()),
                sum(one_dollar_per_10_mol.values()) / 2.0,
            )
        )

        self.assertEqual(npv_electrolyte(60, 0), (0.0, 0))
        self.assertEqual(npv_electrolyte(60, 0.0015 / 100 * 365), (0.0, 0))
        self.assertEqual(npv_electrolyte(10, 0.02, discount=0), (20.0, 2))
        self.assertTrue(isclose(npv_electrolyte(10, 0, -100)[0], -100 / 1.02**20))
        self.assertTrue(isclose(eq6_maintenance(0, 0, 4 * 1.02**10)[0], 1.0))

    def test_dataset_integrity(self):
        emmel = load_emmel_dataset()
        yang = load_yang_dataset()
        n_emmel, n_yang = validate_datasets(emmel, yang)
        self.assertEqual(n_emmel, 38)
        self.assertEqual(n_yang, 19)

    def test_emmel_benchmark_reconciliation(self):
        emmel = load_emmel_dataset()
        summary, components = benchmark_reconciliation(emmel)
        self.assertTrue(
            np.allclose(summary["eq1_from_published"], summary["published"], rtol=0, atol=0.02)
        )
        self.assertAlmostEqual(summary.loc["VRFB", "independent"], 840.62, places=1)
        self.assertAlmostEqual(summary.loc["1,6-DPAP", "independent"], 459.63, places=1)
        self.assertAlmostEqual(summary.loc["TEMPTMA", "independent"], 3826.87, places=1)
        self.assertEqual(int(components.loc["1,6-DPAP", "n_exchanges"]), 0)
        self.assertEqual(int(components.loc["TEMPTMA", "n_exchanges"]), 98)

    def test_yang_phosphate_cost_reduction_and_benchmarks(self):
        o3 = yang_half_cell("O-3", "H3PO4", float(YANG_STAGE_INPUTS.loc["O-3", "solubility_M_phosphate"]))
        asym = yang_half_cell(
            "asym-O-5", "H3PO4", float(YANG_STAGE_INPUTS.loc["asym-O-5", "solubility_M_phosphate"])
        )
        o3_res = evaluate_molecule_cost(o3, 52.70, mode="sensitivity_variant", voltage=YANG_BASE_UD)
        asym_sens = evaluate_molecule_cost(asym, 0.45, mode="sensitivity_variant", voltage=YANG_BASE_UD)
        asym_ind = evaluate_molecule_cost(asym, 0.45, mode="independent", voltage=YANG_BASE_UD)

        reduction_pct = 100.0 * (1.0 - asym_sens["capital_usd_kWh"] / o3_res["capital_usd_kWh"])
        self.assertAlmostEqual(reduction_pct, 99.1, places=1)
        self.assertTrue(10700.0 < asym_sens["capital_usd_kWh"] < 10900.0)
        self.assertTrue(10700.0 < asym_ind["capital_usd_kWh"] < 10900.0)
        self.assertGreater(
            asym_ind["capital_usd_kWh"] - YANG_PUBLISHED_TARGETS["VRFB Future"], 10000.0
        )

    def test_inverse_fade_targets_and_grid_consistency(self):
        asym_sulfate = yang_half_cell("asym-O-5", "H2SO4", 0.44)
        for mode in ("sensitivity_variant", "independent"):
            lim_vrfb = max_fade_for_cost(
                asym_sulfate, YANG_PUBLISHED_TARGETS["VRFB Future"], YANG_BASE_UD, mode=mode
            )
            self.assertTrue(0.0025 <= lim_vrfb <= 0.0060)
            attained = evaluate_molecule_cost(
                asym_sulfate, lim_vrfb, mode=mode, voltage=YANG_BASE_UD
            )["capital_usd_kWh"]
            above = evaluate_molecule_cost(
                asym_sulfate, lim_vrfb + 1e-8, mode=mode, voltage=YANG_BASE_UD
            )["capital_usd_kWh"]
            self.assertLessEqual(attained, YANG_PUBLISHED_TARGETS["VRFB Future"])
            self.assertGreater(above, YANG_PUBLISHED_TARGETS["VRFB Future"])

        v_grid = np.linspace(1.0, 2.0, 21)
        f_grid = np.geomspace(0.001, 2.0, 30)
        grid = capital_cost_grid(asym_sulfate, v_grid, f_grid, mode="independent")
        self.assertTrue((np.diff(grid, axis=0) >= -1e-9).all())
        self.assertTrue((np.diff(grid, axis=1) < 0).all())
        point = evaluate_molecule_cost(
            asym_sulfate, float(f_grid[10]), mode="independent", voltage=float(v_grid[8])
        )["capital_usd_kWh"]
        self.assertAlmostEqual(float(grid[10, 8]), point, places=6)

    def test_protocol_and_polymorph_sensitivity_audit(self):
        asym_sulfate = yang_half_cell("asym-O-5", "H2SO4", 0.44)
        proto = protocol_sensitivity_summary(
            asym_sulfate, ps_fade_pct_day=0.46, gs_fade_pct_day=1.80, mode="independent"
        )
        self.assertAlmostEqual(proto["delta_gs_minus_ps_usd_kWh"], 33363.0, delta=150.0)

        o5_fresh = yang_half_cell("O-5", "H2SO4", 0.24)
        poly_ps = polymorph_sensitivity_summary(
            o5_fresh, aged_solubility_M=0.06, fade_pct_day=0.78, mode="independent"
        )
        self.assertAlmostEqual(poly_ps["delta_aged_minus_fresh_usd_kWh"], 50262.96, places=2)
        variant_poly_ps = polymorph_sensitivity_summary(
            o5_fresh, aged_solubility_M=0.06, fade_pct_day=0.78, mode="sensitivity_variant"
        )
        self.assertAlmostEqual(variant_poly_ps["delta_aged_minus_fresh_usd_kWh"], 50332.63, places=2)


if __name__ == "__main__":
    unittest.main()
