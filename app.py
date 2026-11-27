"""Interactive Streamlit web calculator for AORFB cost & protocol sensitivity auditing."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from aorfb_audit import (
    CARRIER_PRESETS,
    YANG_BASE_UD,
    YANG_ORGANIC_PRICE_USD_KG,
    YANG_PUBLISHED_TARGETS,
    YANG_STACK_HIGH,
    YANG_STACK_LOW,
    YANG_STACK_MID,
    YANG_STAGE_INPUTS,
    HalfCell,
    benchmark_reconciliation,
    capital_cost_grid,
    evaluate_molecule_cost,
    load_emmel_dataset,
    load_yang_dataset,
    max_fade_for_cost,
    polymorph_sensitivity_summary,
    protocol_sensitivity_summary,
    yang_half_cell,
)

st.set_page_config(
    page_title="AORFB Cost & Protocol Audit",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Preset scenarios grounded in Yang et al. (2021) & Navarro Garcia et al. (2021)
MOLECULE_PRESETS = {
    "asym-O-5 (1 M H2SO4, PS)": dict(
        name="asym-O-5",
        mass=336.43,
        electrons=2,
        carrier="1 M H2SO4",
        solubility=0.44,
        voltage=1.40,
        fade=0.46,
        protocol="Potentiostatic (PS)",
        ps_fade=0.46,
        gs_fade=1.80,
        aged_sol=0.44,
    ),
    "asym-O-5 (1 M H3PO4, PS)": dict(
        name="asym-O-5",
        mass=336.43,
        electrons=2,
        carrier="1 M H3PO4",
        solubility=0.55,
        voltage=1.40,
        fade=0.45,
        protocol="Potentiostatic (PS)",
        ps_fade=0.45,
        gs_fade=1.80,
        aged_sol=0.55,
    ),
    "O-5 fresh (1 M H2SO4, PS)": dict(
        name="O-5 (fresh)",
        mass=336.43,
        electrons=2,
        carrier="1 M H2SO4",
        solubility=0.24,
        voltage=1.40,
        fade=0.78,
        protocol="Potentiostatic (PS)",
        ps_fade=0.78,
        gs_fade=1.10,
        aged_sol=0.06,
    ),
    "O-5 aged polymorph (0.06 M, PS)": dict(
        name="O-5 (aged polymorph)",
        mass=336.43,
        electrons=2,
        carrier="1 M H2SO4",
        solubility=0.06,
        voltage=1.40,
        fade=0.78,
        protocol="Potentiostatic (PS)",
        ps_fade=0.78,
        gs_fade=1.10,
        aged_sol=0.06,
    ),
    "O-1 (1 M H2SO4, PS)": dict(
        name="O-1",
        mass=251.33,
        electrons=2,
        carrier="1 M H2SO4",
        solubility=0.48,
        voltage=1.40,
        fade=5.70,
        protocol="Potentiostatic (PS)",
        ps_fade=5.70,
        gs_fade=2.90,
        aged_sol=0.48,
    ),
    "O-3 (1 M H3PO4, PS — stress case)": dict(
        name="O-3",
        mass=308.38,
        electrons=2,
        carrier="1 M H3PO4",
        solubility=0.53,
        voltage=1.40,
        fade=52.70,
        protocol="Potentiostatic (PS)",
        ps_fade=52.70,
        gs_fade=48.50,
        aged_sol=0.53,
    ),
    "Custom candidate molecule": dict(
        name="Custom candidate",
        mass=300.0,
        electrons=2,
        carrier="1 M H2SO4",
        solubility=1.00,
        voltage=1.40,
        fade=0.05,
        protocol="Unspecified",
        ps_fade=0.05,
        gs_fade=0.20,
        aged_sol=0.50,
    ),
}


@st.cache_data
def get_datasets():
    emmel = load_emmel_dataset()
    yang = load_yang_dataset()
    summary, components = benchmark_reconciliation(emmel)
    return emmel, yang, summary, components


emmel_df, yang_df, recon_summary, recon_components = get_datasets()

# ---------------------------------------------------------------------------
# Sidebar: Candidate Molecule & Protocol Configuration
# ---------------------------------------------------------------------------
st.sidebar.title("Candidate Parameters")

preset_label = st.sidebar.selectbox(
    "Load literature preset",
    options=list(MOLECULE_PRESETS.keys()),
    index=0,
)
preset = MOLECULE_PRESETS[preset_label]

mode_label = st.sidebar.radio(
    "Maintenance & discounting convention",
    options=[
        "Independent (full-stack, threshold-crossing NPV)",
        "Sensitivity variant (single-stack, Eq. S38 t-1 index)",
    ],
    index=0,
    help=(
        "Independent replaces all stacks at year 10 and discounts the first electrolyte "
        "renewal at threshold crossing. Sensitivity variant reproduces the printed Eq. S38 "
        "t-1 index and single-stack maintenance diagnostic."
    ),
)
mode = "independent" if mode_label.startswith("Independent") else "sensitivity_variant"

st.sidebar.markdown("---")
mol_name = st.sidebar.text_input("Molecule label", value=preset["name"])
mass_g_mol = st.sidebar.number_input(
    "Molar mass Mw (g/mol, neutral basis)",
    min_value=50.0,
    max_value=1500.0,
    value=float(preset["mass"]),
    step=1.0,
)
electrons = st.sidebar.selectbox(
    "Electrons transferred (n)",
    options=[1, 2],
    index=0 if preset["electrons"] == 1 else 1,
)
carrier_name = st.sidebar.selectbox(
    "Supporting electrolyte carrier",
    options=list(CARRIER_PRESETS.keys()),
    index=list(CARRIER_PRESETS.keys()).index(preset["carrier"]),
)
solubility_M = st.sidebar.number_input(
    "Active solubility c (mol/L)",
    min_value=0.01,
    max_value=6.0,
    value=float(preset["solubility"]),
    step=0.02,
    format="%.3f",
)
voltage_V = st.sidebar.slider(
    "Full-cell discharge voltage Ud (V)",
    min_value=0.80,
    max_value=2.00,
    value=float(preset["voltage"]),
    step=0.01,
)
fade_pct_day = st.sidebar.number_input(
    "Daily capacity fade (%/day)",
    min_value=0.0001,
    max_value=100.0,
    value=float(preset["fade"]),
    step=0.01,
    format="%.4f",
)
protocol_tag = st.sidebar.selectbox(
    "Reported cycling protocol",
    options=["Potentiostatic (PS)", "Galvanostatic (GS)", "Unspecified"],
    index=["Potentiostatic (PS)", "Galvanostatic (GS)", "Unspecified"].index(
        preset["protocol"]
    ),
)

with st.sidebar.expander("Advanced techno-economic inputs", expanded=False):
    active_price = st.number_input(
        "Organic active price (USD/kg)",
        min_value=0.50,
        max_value=50.0,
        value=float(YANG_ORGANIC_PRICE_USD_KG),
        step=0.10,
    )
    stacks = st.slider(
        "WP1 stack count (TEMPTMA–DHDMBS range)",
        min_value=float(YANG_STACK_LOW),
        max_value=float(YANG_STACK_HIGH),
        value=float(YANG_STACK_MID),
        step=0.5,
    )
    soc_window = st.slider("Accessible SOC window", 0.50, 0.95, 0.80, 0.05)

candidate = HalfCell(
    mass_g_mol=float(mass_g_mol),
    electrons=int(electrons),
    active_mol_L=float(solubility_M),
    active_usd_kg=float(active_price),
    soc=float(soc_window),
    **CARRIER_PRESETS[carrier_name],
)

res = evaluate_molecule_cost(
    candidate,
    fade_pct_day,
    mode=mode,
    voltage=voltage_V,
    stacks=stacks,
)

# ---------------------------------------------------------------------------
# Header & Scientific Boundary Warnings
# ---------------------------------------------------------------------------
st.title("AORFB Cost & Cycling-Protocol Audit Calculator")
st.caption(
    "Techno-economic benchmarking (Emmel et al., Nat. Commun. 2023) applied to aqueous "
    "organic redox flow battery active materials and cycling-protocol sensitivity "
    "(Yang et al., Molecules 2021; Navarro García et al., Molecules 2021)."
)

# Carbon HER OCV guide for asym-O-5 posolyte (~0.65 + 0.20 + 0.60 = 1.45 V)
CARBON_OCV_GUIDE = 1.45
warnings = []
if voltage_V > CARBON_OCV_GUIDE:
    warnings.append(
        f"**Exploratory voltage regime ({voltage_V:.2f} V > 1.45 V):** "
        "Exceeds the approximate 1.45 V open-circuit guide for acidic hydroquinone posolytes "
        "paired against a carbon hydrogen-evolution bound (−0.60 V vs. SHE at pH 0)."
    )
if fade_pct_day >= 1.0:
    warnings.append(
        f"**Rapid-fade extrapolation ({fade_pct_day:.2f}%/day):** "
        "Laboratory symmetric-cell fade extrapolated to saturation and 20-year operation "
        "should be interpreted as a stress-case cost amplification."
    )
if protocol_tag == "Unspecified":
    warnings.append(
        "**Protocol not specified:** PS vs. GS cycling, concentration, and SOC window "
        "differences can shift modeled capital cost by >$30,000/kWh for the same molecule."
    )
for w in warnings:
    st.warning(w)

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab1, tab2, tab3 = st.tabs(
    [
        "1. Single-Molecule Audit",
        "2. 42-Molecule Benchmark & Contours",
        "3. Literature Audit & Reproducibility",
    ]
)

# =========================== TAB 1 =========================================
with tab1:
    vrfb_target = YANG_PUBLISHED_TARGETS["VRFB Future"]
    ferro_target = YANG_PUBLISHED_TARGETS["ammonium ferrocyanide"]

    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Total Capital Cost", f"${res['capital_usd_kWh']:,.1f} / kWh")
    k2.metric(
        "20-yr Electrolyte Renewals",
        f"{res['n_exchanges']:,}",
        help=f"Linearized annual replacement fraction f = {res['annual_f']:.3f} yr^-1",
    )
    k3.metric(
        "Gap vs. VRFB Future ($528.5)",
        f"{res['capital_usd_kWh'] - vrfb_target:+,.1f} / kWh",
    )
    k4.metric(
        "Gap vs. Ferrocyanide ($651.7)",
        f"{res['capital_usd_kWh'] - ferro_target:+,.1f} / kWh",
    )

    c_left, c_right = st.columns([1.25, 1.0])

    with c_left:
        st.subheader("Capital Cost Breakdown (USD/kWh)")
        breakdown_df = pd.DataFrame(
            {
                "Component": [
                    "Active material (Eq. 5)",
                    "Supporting salt (Eq. 5)",
                    "Solvent (Eq. 5)",
                    "Tanks (Eq. 5)",
                    "Power subsystem (Eq. 2–3 / 4 h)",
                    "10-yr stack maintenance NPV (Eq. 6)",
                    "20-yr electrolyte renewal NPV (Eq. 6)",
                ],
                "USD/kWh": [
                    res["active_usd_kWh"],
                    res["salt_usd_kWh"],
                    res["solvent_usd_kWh"],
                    res["tank_usd_kWh"],
                    res["power_usd_kWh"],
                    res["stack_npv_usd_kWh"],
                    res["exchange_npv_usd_kWh"],
                ],
            }
        )
        breakdown_df["Share (%)"] = (
            100.0 * breakdown_df["USD/kWh"] / res["capital_usd_kWh"]
        )

        fig_bar = go.Figure(
            go.Bar(
                y=breakdown_df["Component"][::-1],
                x=breakdown_df["USD/kWh"][::-1],
                orientation="h",
                marker_color=[
                    "#c73329",
                    "#4a6fa5",
                    "#2b5876",
                    "#6c8ebf",
                    "#8faadc",
                    "#b4c7e7",
                    "#164273",
                ],
                text=[f"${v:,.1f}" for v in breakdown_df["USD/kWh"][::-1]],
                textposition="auto",
            )
        )
        fig_bar.update_layout(
            height=340,
            margin=dict(l=10, r=20, t=10, b=30),
            xaxis_title="Cost (USD/kWh)",
            yaxis_title="",
            template="plotly_white",
        )
        st.plotly_chart(fig_bar, use_container_width=True)

    with c_right:
        st.subheader("Inverse Stability Targets")
        max_fade_ferro = max_fade_for_cost(
            candidate, ferro_target, voltage_V, mode=mode, stacks=stacks
        )
        max_fade_vrfb = max_fade_for_cost(
            candidate, vrfb_target, voltage_V, mode=mode, stacks=stacks
        )

        target_rows = []
        for label, tgt, lim in [
            ("Ammonium ferrocyanide ($651.7/kWh)", ferro_target, max_fade_ferro),
            ("VRFB Future ($528.5/kWh)", vrfb_target, max_fade_vrfb),
        ]:
            if np.isfinite(lim) and lim > 0:
                factor = fade_pct_day / lim
                target_rows.append(
                    {
                        "Benchmark": label,
                        "Max Fade (%/day)": f"≤ {lim:.5f}",
                        "Improvement Needed": (
                            f"{factor:,.1f}x slower" if factor > 1.0 else "Already meets target"
                        ),
                    }
                )
            else:
                target_rows.append(
                    {
                        "Benchmark": label,
                        "Max Fade (%/day)": "Unreachable at 0% fade",
                        "Improvement Needed": "Requires higher Ud or lower fixed cost",
                    }
                )
        st.dataframe(pd.DataFrame(target_rows), use_container_width=True, hide_index=True)

        st.markdown("**Component Table**")
        st.dataframe(
            breakdown_df.style.format({"USD/kWh": "{:,.2f}", "Share (%)": "{:.2f}%"}),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("---")
    st.subheader("Protocol & Polymorph Stress Test")
    st.caption(
        "Quantify how switching cycling protocol (PS vs. GS) or losing solubility through "
        "polymorph aging shifts system capital cost at identical voltage and stack count."
    )
    sc1, sc2 = st.columns(2)

    with sc1:
        st.markdown("**A. Cycling Protocol Bias (PS vs. GS)**")
        p1, p2 = st.columns(2)
        ps_fade_in = p1.number_input(
            "PS fade (%/day)",
            min_value=0.0001,
            max_value=100.0,
            value=float(preset["ps_fade"]),
            step=0.01,
            format="%.4f",
            key="ps_stress",
        )
        gs_fade_in = p2.number_input(
            "GS fade (%/day)",
            min_value=0.0001,
            max_value=100.0,
            value=float(preset["gs_fade"]),
            step=0.01,
            format="%.4f",
            key="gs_stress",
        )
        proto_res = protocol_sensitivity_summary(
            candidate, ps_fade_in, gs_fade_in, mode=mode, voltage=voltage_V, stacks=stacks
        )
        m1, m2, m3 = st.columns(3)
        m1.metric("PS Total", f"${proto_res['ps_usd_kWh']:,.0f}")
        m2.metric("GS Total", f"${proto_res['gs_usd_kWh']:,.0f}")
        m3.metric(
            "GS − PS Shift",
            f"{proto_res['delta_gs_minus_ps_usd_kWh']:+,.0f} $/kWh",
        )

    with sc2:
        st.markdown("**B. Polymorph Solubility Aging (Fixed Fade)**")
        s1, s2 = st.columns(2)
        fresh_sol_in = s1.number_input(
            "Fresh solubility (M)",
            min_value=0.01,
            max_value=6.0,
            value=float(solubility_M),
            step=0.02,
            format="%.3f",
            key="fresh_sol",
        )
        aged_sol_in = s2.number_input(
            "Aged polymorph solubility (M)",
            min_value=0.01,
            max_value=6.0,
            value=float(min(preset["aged_sol"], fresh_sol_in)),
            step=0.01,
            format="%.3f",
            key="aged_sol",
        )
        fresh_side = replace(candidate, active_mol_L=fresh_sol_in)
        poly_res = polymorph_sensitivity_summary(
            fresh_side, aged_sol_in, fade_pct_day, mode=mode, voltage=voltage_V, stacks=stacks
        )
        pm1, pm2, pm3 = st.columns(3)
        pm1.metric("Fresh Total", f"${poly_res['fresh_total_usd_kWh']:,.0f}")
        pm2.metric("Aged Total", f"${poly_res['aged_total_usd_kWh']:,.0f}")
        pm3.metric(
            "Aged − Fresh Shift",
            f"{poly_res['delta_aged_minus_fresh_usd_kWh']:+,.0f} $/kWh",
        )

    # Export audit CSV
    audit_export = pd.DataFrame(
        [
            {
                "molecule": mol_name,
                "mode": mode,
                "protocol": protocol_tag,
                "mass_g_mol": mass_g_mol,
                "electrons": electrons,
                "carrier": carrier_name,
                "solubility_M": solubility_M,
                "voltage_V": voltage_V,
                "fade_pct_day": fade_pct_day,
                "eq5_usd_kWh": round(res["eq5_usd_kWh"], 2),
                "power_usd_kWh": round(res["power_usd_kWh"], 2),
                "exchange_npv_usd_kWh": round(res["exchange_npv_usd_kWh"], 2),
                "stack_npv_usd_kWh": round(res["stack_npv_usd_kWh"], 2),
                "n_exchanges": res["n_exchanges"],
                "capital_usd_kWh": round(res["capital_usd_kWh"], 2),
                "max_fade_vrfb_pct_day": round(max_fade_vrfb, 6)
                if np.isfinite(max_fade_vrfb)
                else None,
                "max_fade_ferrocyanide_pct_day": round(max_fade_ferro, 6)
                if np.isfinite(max_fade_ferro)
                else None,
                "protocol_delta_gs_minus_ps_usd_kWh": round(
                    proto_res["delta_gs_minus_ps_usd_kWh"], 2
                ),
                "polymorph_delta_aged_minus_fresh_usd_kWh": round(
                    poly_res["delta_aged_minus_fresh_usd_kWh"], 2
                ),
            }
        ]
    )
    st.download_button(
        label="Download Single-Molecule Audit Summary (CSV)",
        data=audit_export.to_csv(index=False).encode("utf-8"),
        file_name=f"{mol_name.lower().replace(' ', '_')}_aorfb_audit.csv",
        mime="text/csv",
    )

# =========================== TAB 2 =========================================
with tab2:
    st.subheader("42-Molecule Benchmark Map (Emmel 38 + Yang 4 Hydroquinones)")
    st.caption(
        "Compares published Future Case capital costs from Emmel et al. (2023) against "
        "modeled Yang et al. (2021) phosphate-PS hydroquinones (1.40 V) and your active candidate."
    )

    # Compute Yang 4 phosphate hydroquinones under current mode
    yang_points = []
    for name, row in YANG_STAGE_INPUTS.iterrows():
        hc = yang_half_cell(name, "H3PO4", float(row["solubility_M_phosphate"]))
        y_res = evaluate_molecule_cost(
            hc,
            float(row["fade_pct_day_phosphate_PS"]),
            mode=mode,
            voltage=YANG_BASE_UD,
            stacks=stacks,
        )
        yang_points.append(
            {
                "molecule": f"{name} (Yang H3PO4 PS)",
                "fade_pct_day": float(row["fade_pct_day_phosphate_PS"]),
                "capital_usd_kWh": y_res["capital_usd_kWh"],
                "solubility_M": float(row["solubility_M_phosphate"]),
                "voltage_V": YANG_BASE_UD,
            }
        )
    yang_pts_df = pd.DataFrame(yang_points)

    fig_bench = go.Figure()
    for role, color in [("negolyte", "#164273"), ("posolyte", "#c73329")]:
        sub = emmel_df[emmel_df["electrode_role"] == role]
        fig_bench.add_trace(
            go.Scatter(
                x=sub["fade_pct_day"],
                y=sub["capital_usd_kWh"],
                mode="markers",
                name=f"Emmel {role} (n={len(sub)})",
                marker=dict(size=9, color=color, opacity=0.82),
                text=sub["molecule"],
                customdata=np.stack(
                    [
                        sub["material_family"],
                        sub["solubility_mol_L"],
                        sub["discharge_voltage_V"],
                    ],
                    axis=-1,
                ),
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "Family: %{customdata[0]}<br>"
                    "Fade: %{x:.4g} %/day<br>"
                    "Capital cost: $%{y:,.1f}/kWh<br>"
                    "Solubility: %{customdata[1]:.2f} M | Ud: %{customdata[2]:.2f} V"
                    "<extra></extra>"
                ),
            )
        )

    fig_bench.add_trace(
        go.Scatter(
            x=yang_pts_df["fade_pct_day"],
            y=yang_pts_df["capital_usd_kWh"],
            mode="markers+text",
            name="Yang hydroquinones (H3PO4 PS, 1.40 V)",
            text=[m.split(" ")[0] for m in yang_pts_df["molecule"]],
            textposition="top center",
            marker=dict(size=11, symbol="diamond", color="#e07a22", line=dict(width=1, color="#1f2933")),
            hovertemplate=(
                "<b>%{text}</b> (Yang H3PO4 PS)<br>"
                "Fade: %{x:.3g} %/day<br>"
                "Modeled cost: $%{y:,.1f}/kWh<extra></extra>"
            ),
        )
    )

    fig_bench.add_trace(
        go.Scatter(
            x=[fade_pct_day],
            y=[res["capital_usd_kWh"]],
            mode="markers",
            name=f"Active candidate: {mol_name}",
            marker=dict(
                size=16,
                symbol="star",
                color="#00a878",
                line=dict(width=1.5, color="#0b1f1a"),
            ),
            hovertemplate=(
                f"<b>{mol_name}</b> ({protocol_tag})<br>"
                "Fade: %{x:.4g} %/day<br>"
                "Modeled cost: $%{y:,.1f}/kWh<extra></extra>"
            ),
        )
    )

    fig_bench.add_hline(
        y=YANG_PUBLISHED_TARGETS["VRFB Future"],
        line_dash="dot",
        line_color="#3a4653",
        annotation_text="VRFB Future ($528.5/kWh)",
        annotation_position="bottom right",
    )
    fig_bench.add_hline(
        y=YANG_PUBLISHED_TARGETS["ammonium ferrocyanide"],
        line_dash="dash",
        line_color="#6b7280",
        annotation_text="Ammonium ferrocyanide ($651.7/kWh)",
        annotation_position="top right",
    )
    fig_bench.update_layout(
        height=470,
        xaxis_type="log",
        yaxis_type="log",
        xaxis_title="Daily capacity fade (%/day, log scale)",
        yaxis_title="Capital cost (USD/kWh, log scale)",
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(l=20, r=20, t=40, b=30),
    )
    st.plotly_chart(fig_bench, use_container_width=True)

    st.markdown("---")
    show_contours = st.checkbox(
        f"Show 2D inverse cost-target contour map for {mol_name} (Voltage vs. Daily Fade)",
        value=True,
    )
    if show_contours:
        v_grid = np.linspace(1.0, 2.0, 81)
        f_grid = np.geomspace(0.001, 2.0, 120)
        z_grid = capital_cost_grid(candidate, v_grid, f_grid, mode=mode, stacks=stacks)
        log_z = np.log10(np.clip(z_grid, 200.0, 2e5))

        fig_cnt = go.Figure(
            data=go.Heatmap(
                x=v_grid,
                y=f_grid,
                z=log_z,
                colorscale="YlGnBu",
                colorbar=dict(
                    title="Cost ($/kWh)",
                    tickvals=[np.log10(v) for v in [500, 1000, 3000, 10000, 30000]],
                    ticktext=["500", "1,000", "3,000", "10,000", "30,000"],
                ),
                customdata=z_grid,
                hovertemplate=(
                    "Ud: %{x:.2f} V<br>Fade: %{y:.4f} %/day<br>"
                    "Capital cost: $%{customdata:,.0f}/kWh<extra></extra>"
                ),
            )
        )

        # Overlay red inverse target curves
        v_curve = np.linspace(1.0, 2.0, 41)
        for bench_name, dash_style in [
            ("VRFB Future", "solid"),
            ("ammonium ferrocyanide", "dash"),
        ]:
            tgt_val = YANG_PUBLISHED_TARGETS[bench_name]
            lim_curve = [
                max_fade_for_cost(candidate, tgt_val, float(v), mode=mode, stacks=stacks)
                for v in v_curve
            ]
            fig_cnt.add_trace(
                go.Scatter(
                    x=v_curve,
                    y=lim_curve,
                    mode="lines",
                    name=f"${tgt_val:g}/kWh ({bench_name})",
                    line=dict(color="#cf2026", width=2.5, dash=dash_style),
                )
            )

        fig_cnt.add_vline(
            x=CARBON_OCV_GUIDE,
            line_dash="dot",
            line_color="#1f2933",
            annotation_text="1.45 V carbon OCV guide",
            annotation_position="top right",
        )
        if 1.0 <= voltage_V <= 2.0 and 0.001 <= fade_pct_day <= 2.0:
            fig_cnt.add_trace(
                go.Scatter(
                    x=[voltage_V],
                    y=[fade_pct_day],
                    mode="markers",
                    name=f"{mol_name} operating point",
                    marker=dict(
                        size=12,
                        symbol="diamond",
                        color="white",
                        line=dict(width=2, color="#172d3f"),
                    ),
                )
            )
        fig_cnt.update_layout(
            height=440,
            yaxis_type="log",
            xaxis_title="Full-cell discharge voltage Ud (V)",
            yaxis_title="Daily capacity fade (%/day, log scale)",
            template="plotly_white",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
            margin=dict(l=20, r=20, t=40, b=30),
        )
        st.plotly_chart(fig_cnt, use_container_width=True)

    with st.expander("Explore full 38-molecule Emmel benchmark table", expanded=False):
        role_filter = st.multiselect(
            "Filter by electrode role",
            options=["negolyte", "posolyte"],
            default=["negolyte", "posolyte"],
        )
        display_cols = [
            "molecule",
            "material_family",
            "electrode_role",
            "molar_mass_g_mol",
            "n_electrons",
            "solubility_mol_L",
            "fade_pct_day",
            "discharge_voltage_V",
            "capital_usd_kWh",
        ]
        filtered = emmel_df[emmel_df["electrode_role"].isin(role_filter)][display_cols]
        st.dataframe(filtered, use_container_width=True, hide_index=True)

# =========================== TAB 3 =========================================
with tab3:
    st.subheader("Published Benchmark Reconciliation (Emmel et al., 2023)")
    st.caption(
        "Independent recalculation of Eqs. 1, 5, and 6 from documented inputs compared against "
        "published Present Case VRFB and Future Case 1,6-DPAP / TEMPTMA totals."
    )
    st.dataframe(
        recon_summary.reset_index().style.format(
            {
                "published": "{:,.2f}",
                "eq1_from_published": "{:,.2f}",
                "independent": "{:,.2f}",
                "independent_error": "{:+,.2f}",
                "sensitivity_variant": "{:,.2f}",
                "sensitivity_error": "{:+,.2f}",
            }
        ),
        use_container_width=True,
        hide_index=True,
    )

    st.subheader("Key Methodological Audit Findings & Scope Boundaries")
    st.markdown(
        """
- **Potential Reference Offset (+0.60 V):** Emmel Fig. 2 plots redox potentials roughly `+0.60 V` above conventional RHE values for `1,6-DPAP` and `TEMPTMA`, consistent with referencing potentials to the `-0.60 V` carbon hydrogen-evolution bound at pH 0 (Eq. S42). Both digitized and conventional RHE values are retained in `data/emmel_38_molecules_audit.csv`.
- **Eq. S40 Maintenance Diagnostic:** Charging single-stack replacement alongside the printed Eq. S38 `t - 1` discount index closely reproduces two published maintenance totals, whereas full-stack replacement at year 10 increases maintenance NPV. Both modes are exposed in this calculator because the original published implementation was not verified directly.
- **Pricing & Extrapolation Caveats:** Active organic materials are priced on a neutral-molecule mass basis (`$3.48/kg`) paired with an idealized counter-electrolyte (`100 g/mol`, `1 e-`, `2 M`). Fade rates measured in dilute symmetric laboratory cells (`0.05–0.10 M`) are extrapolated to saturation solubility and 20-year operation as a scenario audit, not a direct field lifetime prediction.
- **Protocol & Polymorph Sensitivity:** For matched-formulation `asym-O-5` (`0.10 M`, `1 M H2SO4`), switching potentiostatic (`0.46%/day`) to galvanostatic (`1.8%/day`) fade increases modeled system cost by `+$33,363/kWh`. For `O-5`, aging-induced polymorph solubility loss (`0.24 M -> 0.06 M`) adds `+$50,333/kWh` at fixed PS fade.
        """
    )

    st.subheader("Primary Sources & Attribution")
    st.markdown(
        """
1. **Yang et al. (2021)**, *Molecules* **26**(13), 3823. [https://doi.org/10.3390/molecules26133823](https://doi.org/10.3390/molecules26133823)
2. **Navarro García et al. (2021)**, *Molecules* **26**(5), 1203. [https://doi.org/10.3390/molecules26051203](https://doi.org/10.3390/molecules26051203)
3. **Emmel et al. (2023)**, *Nature Communications* **14**, 6672. [https://doi.org/10.1038/s41467-023-42450-9](https://doi.org/10.1038/s41467-023-42450-9)
        """
    )
