"""Core techno-economic and electrochemistry protocol-audit model for AORFBs.

Implements Eqs. 1, 2–3, 5 (S32–S36), and 6 (S37–S41) from Emmel et al.,
Nature Communications 14, 6672 (2023), together with hydroquinone posolyte
benchmarking, inverse fade-target bisection, and PS/GS protocol & polymorph
sensitivity audits using Yang et al., Molecules 26, 3823 (2021) and
Navarro García et al., Molecules 26, 1203 (2021).
"""

from dataclasses import dataclass, replace
from math import floor
from typing import Dict, Iterable, Sequence, Tuple

import numpy as np
import pandas as pd

# Physical and techno-economic model constants (Emmel et al., 2023)
FARADAY_C_PER_MOL = 96485.33212
J_PER_KWH = 3.6e6
DISCHARGE_H = 4
SYSTEM_KW = 1000
SYSTEM_EFFICIENCY = 0.94
COULOMBIC_EFFICIENCY = 0.97
DISCOUNT = 0.02
LIFE_YEARS = 20
STACK_YEARS = 10
REPLACEMENT_THRESHOLD = 0.80
INVERTER_USD_KW = 103
PUMP_USD_ITEM = 205
PUMPS_PER_STACK = 2

# Yang hydroquinone scenario constants
YANG_ORGANIC_PRICE_USD_KG = 3.48
YANG_STACK_PRICE_USD = 993.82
YANG_BASE_UD = 1.40
YANG_STACK_LOW = 430.17   # TEMPTMA WP1 Future Case stack count
YANG_STACK_HIGH = 438.64  # DHDMBS WP1 Future Case stack count
YANG_STACK_MID = (YANG_STACK_LOW + YANG_STACK_HIGH) / 2.0

YANG_PUBLISHED_TARGETS: Dict[str, float] = {
    "ammonium ferrocyanide": 651.7,
    "VRFB Future": 528.5,
    "VRFB Present": 676.7,
}

# Supporting electrolyte carrier presets
CARRIER_PRESETS: Dict[str, Dict[str, float]] = {
    "1 M H2SO4": dict(carrier_kg_L=1.06, salt_mol_L=1.0, salt_g_mol=98.079, salt_usd_kg=0.17),
    "1 M H3PO4": dict(carrier_kg_L=1.05, salt_mol_L=1.0, salt_g_mol=97.994, salt_usd_kg=0.50),
    "1 M KCl": dict(carrier_kg_L=1.0434, salt_mol_L=1.0, salt_g_mol=74.5483, salt_usd_kg=0.092),
    "0.3 M NaCl": dict(carrier_kg_L=1.01, salt_mol_L=0.3, salt_g_mol=58.443, salt_usd_kg=0.080),
}

# Short alias map matching notebook keys
YANG_CARRIERS: Dict[str, Dict[str, float]] = {
    "H3PO4": CARRIER_PRESETS["1 M H3PO4"],
    "H2SO4": CARRIER_PRESETS["1 M H2SO4"],
}

YANG_STAGE_INPUTS = pd.DataFrame(
    {
        "neutral_mass_g_mol": [251.33, 308.38, 336.43, 336.43],
        "e_per_mol": [2, 2, 2, 2],
        "ortho_sites_substituted": [4, 2, 4, 4],
        "E_half_V_vs_AgAgCl_sulfate": [0.54, 0.77, 0.64, 0.65],
        "solubility_M_phosphate": [0.48, 0.53, 0.40, 0.55],
        "solubility_M_sulfate": [0.48, np.nan, 0.24, 0.44],
        "fade_pct_day_phosphate_PS": [7.30, 52.70, 0.87, 0.45],
        "fade_pct_day_sulfate_PS": [5.70, np.nan, 0.78, 0.46],
        "fade_pct_day_sulfate_GS": [2.90, np.nan, 1.10, 1.80],
    },
    index=pd.Index(["O-1", "O-3", "O-5", "asym-O-5"], name="molecule"),
)

MODE_INDEX_AS_PRINTED: Dict[str, bool] = {
    "sensitivity_variant": True,
    "independent": False,
}


@dataclass(frozen=True)
class HalfCell:
    """Electrochemical and cost parameters for one flow-battery half-cell."""

    mass_g_mol: float
    electrons: int
    active_mol_L: float
    active_usd_kg: float
    carrier_kg_L: float
    salt_mol_L: float
    salt_g_mol: float
    salt_usd_kg: float
    solvent_usd_kg: float = 0.10
    soc: float = 0.80


def ideal_counter(side: HalfCell) -> HalfCell:
    """Construct Emmel's idealized counter-electrolyte (100 g/mol, 1 e-, 2 M)."""
    return HalfCell(
        mass_g_mol=100.0,
        electrons=1,
        active_mol_L=2.0,
        active_usd_kg=side.active_usd_kg,
        carrier_kg_L=side.carrier_kg_L,
        salt_mol_L=side.salt_mol_L,
        salt_g_mol=side.salt_g_mol,
        salt_usd_kg=side.salt_usd_kg,
        solvent_usd_kg=side.solvent_usd_kg,
        soc=side.soc,
    )


def eq1_capital(
    electrolyte_usd_kWh: float | np.ndarray,
    power_usd_kW: float | np.ndarray,
    maintenance_usd_kWh: float | np.ndarray,
    hours: float = DISCHARGE_H,
) -> float | np.ndarray:
    """Emmel Eq. 1: total system capital cost in USD/kWh."""
    return electrolyte_usd_kWh + power_usd_kW / hours + maintenance_usd_kWh


def eq5_electrolyte(
    half_cells: Sequence[HalfCell],
    voltage: float,
    tank_usd_m3: float = 411.0,
    max_fill: float = 0.80,
) -> Dict[str, float]:
    """Emmel Eq. 5 (S32–S36): electrolyte + tank cost breakdown in USD/kWh."""
    if voltage <= 0:
        raise ValueError("Discharge voltage must be positive.")
    totals = {component: 0.0 for component in ("active", "salt", "solvent", "tank")}
    kwh_per_mol_e = (
        FARADAY_C_PER_MOL * voltage * SYSTEM_EFFICIENCY * COULOMBIC_EFFICIENCY / J_PER_KWH
    )
    for side in half_cells:
        if side.electrons <= 0 or side.active_mol_L <= 0 or side.soc <= 0:
            raise ValueError("Electrons, active concentration, and SOC must be positive.")
        active_mol_per_mol_e = 1.0 / (side.electrons * side.soc)
        active_kg_per_mol_e = active_mol_per_mol_e * side.mass_g_mol / 1000.0
        liquid_L_per_mol_e = active_mol_per_mol_e / side.active_mol_L
        salt_kg_L = side.salt_mol_L * side.salt_g_mol / 1000.0
        solvent_kg_L = side.carrier_kg_L - salt_kg_L
        if solvent_kg_L < 0:
            raise ValueError(
                f"Carrier density ({side.carrier_kg_L} kg/L) is below salt mass density ({salt_kg_L:.4f} kg/L)."
            )
        totals["active"] += active_kg_per_mol_e * side.active_usd_kg
        totals["salt"] += liquid_L_per_mol_e * salt_kg_L * side.salt_usd_kg
        totals["solvent"] += liquid_L_per_mol_e * solvent_kg_L * side.solvent_usd_kg
        totals["tank"] += tank_usd_m3 / (side.active_mol_L * 1000.0 * max_fill)
    return {part: cost / kwh_per_mol_e for part, cost in totals.items()}


def power_cost(stacks: float, unit_stack_usd: float) -> Tuple[float, float]:
    """Emmel Eqs. 2–3 / Tables S1 & S3: power subsystem cost (USD/kW) and pump cost."""
    pump_usd_kW = PUMPS_PER_STACK * PUMP_USD_ITEM * stacks / SYSTEM_KW
    # Table S3: heat exchanger, inverter, piping, cabling, process control, pumps.
    bop_usd_kW = 11.0 + INVERTER_USD_KW + 25.0 + 2.0 + 153960.0 / SYSTEM_KW + pump_usd_kW
    stack_usd_kW = stacks * unit_stack_usd / SYSTEM_KW
    return stack_usd_kW + bop_usd_kW + 87.5, pump_usd_kW


def npv_electrolyte(
    electrolyte_usd_kWh: float,
    annual_f: float,
    residual_usd_kWh: float = 0.0,
    *,
    discount: float = DISCOUNT,
    life: int = LIFE_YEARS,
    threshold: float = REPLACEMENT_THRESHOLD,
    maintenance_fraction: float = 1.0,
    index_as_printed: bool = False,
) -> Tuple[float, int]:
    """Emmel Eq. 6 (S37–S38): NPV of electrolyte replacement events over project life."""
    if annual_f < 0 or not (0.0 <= threshold < 1.0):
        raise ValueError("annual_f must be >= 0 and threshold must be in [0, 1).")
    terminal = residual_usd_kWh / (1.0 + discount) ** life
    if annual_f == 0:
        return terminal, 0
    period_years = (1.0 - threshold) / annual_f
    n_changes = floor(life / period_years + 1e-10)
    index_start = 0 if index_as_printed else 1
    factor = sum(
        (1.0 + discount) ** (-i * period_years)
        for i in range(index_start, index_start + n_changes)
    )
    cost_per_change = residual_usd_kWh + electrolyte_usd_kWh * maintenance_fraction
    return terminal + factor * cost_per_change, n_changes


def eq6_maintenance(
    electrolyte_usd_kWh: float,
    annual_f: float,
    stack_replacement_usd_kW: float,
    **options,
) -> Tuple[float, float, float, int]:
    """Emmel Eq. 6 (S37–S41): total maintenance NPV = electrolyte NPV + 10-yr stack NPV."""
    electrolyte_npv, n_changes = npv_electrolyte(electrolyte_usd_kWh, annual_f, **options)
    discount = options.get("discount", DISCOUNT)
    stack_npv = stack_replacement_usd_kW / DISCHARGE_H / (1.0 + discount) ** STACK_YEARS
    return electrolyte_npv + stack_npv, electrolyte_npv, stack_npv, n_changes


def yang_half_cell(molecule: str, acid: str, solubility_M: float) -> HalfCell:
    """Build a HalfCell for one of the four Yang et al. (2021) hydroquinones."""
    inputs = YANG_STAGE_INPUTS.loc[molecule]
    carrier_key = acid if acid in YANG_CARRIERS else acid.replace("1 M ", "")
    return HalfCell(
        mass_g_mol=float(inputs["neutral_mass_g_mol"]),
        electrons=int(inputs["e_per_mol"]),
        active_mol_L=float(solubility_M),
        active_usd_kg=YANG_ORGANIC_PRICE_USD_KG,
        **YANG_CARRIERS[carrier_key],
    )


def fixed_cost_inputs(
    stacks: float,
    mode: str,
    stack_price_usd: float = YANG_STACK_PRICE_USD,
) -> Tuple[float, float]:
    """Return (power_usd_kW, stack_replacement_usd_kW) for a given audit mode."""
    if mode not in MODE_INDEX_AS_PRINTED:
        raise ValueError(f"Unknown mode '{mode}'. Choose from {list(MODE_INDEX_AS_PRINTED)}.")
    total_power, pump = power_cost(stacks, stack_price_usd)
    replacement_stacks = 1.0 if MODE_INDEX_AS_PRINTED[mode] else float(stacks)
    replacement_usd_kW = (
        replacement_stacks * stack_price_usd / SYSTEM_KW + INVERTER_USD_KW + pump
    )
    return total_power, replacement_usd_kW


def evaluate_molecule_cost(
    half_cell: HalfCell,
    fade_pct_day: float,
    mode: str = "independent",
    *,
    voltage: float = YANG_BASE_UD,
    stacks: float = YANG_STACK_MID,
    stack_price_usd: float = YANG_STACK_PRICE_USD,
) -> Dict[str, float]:
    """Evaluate full system cost breakdown (USD/kWh) for a candidate half-cell vs. ideal counter."""
    power_usd_kW, replacement_usd_kW = fixed_cost_inputs(
        stacks, mode, stack_price_usd=stack_price_usd
    )
    el_parts = eq5_electrolyte((half_cell, ideal_counter(half_cell)), voltage)
    electrolyte = sum(el_parts.values())
    annual_f = fade_pct_day / 100.0 * 365.0
    maintenance, exchange_npv, stack_npv, n_exchanges = eq6_maintenance(
        electrolyte,
        annual_f,
        replacement_usd_kW,
        index_as_printed=MODE_INDEX_AS_PRINTED[mode],
    )
    return dict(
        active_usd_kWh=el_parts["active"],
        salt_usd_kWh=el_parts["salt"],
        solvent_usd_kWh=el_parts["solvent"],
        tank_usd_kWh=el_parts["tank"],
        eq5_usd_kWh=electrolyte,
        power_usd_kWh=power_usd_kW / DISCHARGE_H,
        exchange_npv_usd_kWh=exchange_npv,
        stack_npv_usd_kWh=stack_npv,
        maintenance_usd_kWh=maintenance,
        n_exchanges=int(n_exchanges),
        annual_f=annual_f,
        capital_usd_kWh=float(eq1_capital(electrolyte, power_usd_kW, maintenance)),
    )


def capital_cost_grid(
    half_cell: HalfCell,
    voltages: Sequence[float] | np.ndarray,
    fade_rates_pct_day: Sequence[float] | np.ndarray,
    mode: str = "independent",
    *,
    stacks: float = YANG_STACK_MID,
    stack_price_usd: float = YANG_STACK_PRICE_USD,
) -> np.ndarray:
    """Vectorized 2D capital cost surface over (fade_rates x voltages) in USD/kWh."""
    power_usd_kW, replacement_usd_kW = fixed_cost_inputs(
        stacks, mode, stack_price_usd=stack_price_usd
    )
    el_at_1V = sum(eq5_electrolyte((half_cell, ideal_counter(half_cell)), 1.0).values())
    maintenance = [
        eq6_maintenance(
            1.0,
            float(f) / 100.0 * 365.0,
            replacement_usd_kW,
            index_as_printed=MODE_INDEX_AS_PRINTED[mode],
        )
        for f in fade_rates_pct_day
    ]
    exchange_multipliers = np.array([row[1] for row in maintenance], dtype=float)
    stack_npv = float(maintenance[0][2])
    electrolyte = el_at_1V / np.asarray(voltages, dtype=float)
    return np.asarray(
        eq1_capital(
            electrolyte[None, :],
            power_usd_kW,
            exchange_multipliers[:, None] * electrolyte[None, :] + stack_npv,
        ),
        dtype=float,
    )


def max_fade_for_cost(
    half_cell: HalfCell,
    target_usd_kWh: float,
    voltage: float,
    mode: str = "independent",
    *,
    stacks: float = YANG_STACK_MID,
    stack_price_usd: float = YANG_STACK_PRICE_USD,
    upper_bound_pct_day: float = 2.0,
) -> float:
    """Find the maximum daily fade (%/day) that stays at or below target_usd_kWh via bisection."""
    lower = 0.0
    upper = float(upper_bound_pct_day)
    zero_cost = evaluate_molecule_cost(
        half_cell, lower, mode, voltage=voltage, stacks=stacks, stack_price_usd=stack_price_usd
    )["capital_usd_kWh"]
    if zero_cost > target_usd_kWh:
        return float("nan")
    while (
        evaluate_molecule_cost(
            half_cell, upper, mode, voltage=voltage, stacks=stacks, stack_price_usd=stack_price_usd
        )["capital_usd_kWh"]
        <= target_usd_kWh
    ):
        upper *= 2.0
        if upper > 100.0:
            return float("nan")
    for _ in range(35):
        middle = (lower + upper) / 2.0
        cost = evaluate_molecule_cost(
            half_cell, middle, mode, voltage=voltage, stacks=stacks, stack_price_usd=stack_price_usd
        )["capital_usd_kWh"]
        if cost <= target_usd_kWh:
            lower = middle
        else:
            upper = middle
    return lower


def protocol_sensitivity_summary(
    half_cell: HalfCell,
    ps_fade_pct_day: float,
    gs_fade_pct_day: float,
    mode: str = "independent",
    *,
    voltage: float = YANG_BASE_UD,
    stacks: float = YANG_STACK_MID,
) -> Dict[str, float]:
    """Compare system cost under PS vs. GS fade rates at fixed composition."""
    ps_res = evaluate_molecule_cost(
        half_cell, ps_fade_pct_day, mode, voltage=voltage, stacks=stacks
    )
    gs_res = evaluate_molecule_cost(
        half_cell, gs_fade_pct_day, mode, voltage=voltage, stacks=stacks
    )
    return dict(
        ps_usd_kWh=ps_res["capital_usd_kWh"],
        gs_usd_kWh=gs_res["capital_usd_kWh"],
        delta_gs_minus_ps_usd_kWh=gs_res["capital_usd_kWh"] - ps_res["capital_usd_kWh"],
        ps_exchanges=ps_res["n_exchanges"],
        gs_exchanges=gs_res["n_exchanges"],
    )


def polymorph_sensitivity_summary(
    fresh_half_cell: HalfCell,
    aged_solubility_M: float,
    fade_pct_day: float,
    mode: str = "independent",
    *,
    voltage: float = YANG_BASE_UD,
    stacks: float = YANG_STACK_MID,
) -> Dict[str, float]:
    """Compare system cost between fresh and aged-polymorph solubility at fixed fade."""
    aged_half_cell = replace(fresh_half_cell, active_mol_L=float(aged_solubility_M))
    fresh_res = evaluate_molecule_cost(
        fresh_half_cell, fade_pct_day, mode, voltage=voltage, stacks=stacks
    )
    aged_res = evaluate_molecule_cost(
        aged_half_cell, fade_pct_day, mode, voltage=voltage, stacks=stacks
    )
    return dict(
        fresh_eq5_usd_kWh=fresh_res["eq5_usd_kWh"],
        aged_eq5_usd_kWh=aged_res["eq5_usd_kWh"],
        fresh_total_usd_kWh=fresh_res["capital_usd_kWh"],
        aged_total_usd_kWh=aged_res["capital_usd_kWh"],
        delta_aged_minus_fresh_usd_kWh=(
            aged_res["capital_usd_kWh"] - fresh_res["capital_usd_kWh"]
        ),
    )


def benchmark_reconciliation(emmel: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Reconcile published Emmel et al. (2023) VRFB, 1,6-DPAP, and TEMPTMA benchmarks."""
    by_name = emmel.set_index("molecule")
    v_density = 20.52 * (90.939 / 1000.0) * 1.6 / 2.41
    vanadium = HalfCell(90.939, 1, 1.6, 20.52, v_density, 4.0, 98.079, 0.170)
    dpap = HalfCell(354.37, 2, 1.005, 3.48, 1.0434, 1.0, 74.5483, 0.092)
    temp = HalfCell(249.803, 1, 3.2, 3.48, 1.01, 0.3, 58.443, 0.080)

    benchmarks = [
        (
            "VRFB",
            "Present",
            676.70,
            0.74,
            442.78,
            293.13,
            58.70,
            1299.47,
            1661.02,
            0.0,
            0.0,
            (vanadium, vanadium),
        ),
    ]
    for molecule, side in [("1,6-DPAP", dpap), ("TEMPTMA", temp)]:
        p = by_name.loc[molecule]
        benchmarks.append(
            (
                molecule,
                "Future",
                float(p["capital_usd_kWh"]),
                float(p["discharge_voltage_V"]),
                float(p["stacks"]),
                float(p["electrolyte_usd_kWh"]),
                float(p["maintenance_usd_kWh"]),
                float(p["power_usd_kW"]),
                993.82,
                float(p["f_annual_from_daily_linear"]),
                float(p["f_annual_table_S8"]),
                (side, ideal_counter(side)),
            )
        )

    records, parts = [], []
    for (
        name,
        scenario,
        target,
        voltage,
        stacks,
        published_el,
        published_maint,
        published_power,
        stack_price,
        daily_annualized,
        table_annual,
        half_cells,
    ) in benchmarks:
        electrolyte_parts = eq5_electrolyte(half_cells, voltage)
        electrolyte = sum(electrolyte_parts.values())
        power, pump = power_cost(stacks, stack_price)
        peripherals = INVERTER_USD_KW + pump
        stack_full = stacks * stack_price / SYSTEM_KW + peripherals
        maintenance, repl_npv, stack_npv, events = eq6_maintenance(
            electrolyte, daily_annualized, stack_full
        )
        independent = float(eq1_capital(electrolyte, power, maintenance))
        reconciliation = float(eq1_capital(published_el, published_power, published_maint))
        stack_single = stack_price / SYSTEM_KW + peripherals
        diagnostic_maint = eq6_maintenance(
            electrolyte, table_annual, stack_single, index_as_printed=True
        )[0]
        diagnostic = float(eq1_capital(electrolyte, power, diagnostic_maint))
        records.append(
            dict(
                system=name,
                scenario=scenario,
                published=target,
                eq1_from_published=reconciliation,
                independent=independent,
                independent_error=independent - target,
                sensitivity_variant=diagnostic,
                sensitivity_error=diagnostic - target,
            )
        )
        parts.append(
            dict(
                system=name,
                **electrolyte_parts,
                eq5=electrolyte,
                published_eq5=published_el,
                independent_power_usd_kW=power,
                published_power_usd_kW=published_power,
                exchange_npv=repl_npv,
                stack_npv=stack_npv,
                n_exchanges=events,
                eq6=maintenance,
                published_eq6=published_maint,
                sensitivity_eq6=diagnostic_maint,
            )
        )

    return pd.DataFrame(records).set_index("system"), pd.DataFrame(parts).set_index("system")
