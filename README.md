# AORFB Cost & Cycling-Protocol Audit (`aorfb-cost-audit`)

> **Interactive Web App:** [Open the Streamlit calculator](https://aorfb-cost-audit.streamlit.app/)
>
> **Archival DOI:** `[Paste your Zenodo DOI badge after creating GitHub Release v1.0.0]`

An open-source electrochemistry modeling package, reproducible audit notebook, and browser-accessible **Streamlit calculator** that asks a practical question for aqueous organic redox flow battery (AORFB) design: **how much must positive-electrolyte capacity fade improve for substituted hydroquinones to approach inorganic cost benchmarks, and how severely do inconsistent laboratory protocols (PS vs. GS, concentration, polymorph aging) skew system-level rankings?**

This project independently audits the techno-economic model of [Emmel et al. (2023)](https://doi.org/10.1038/s41467-023-42450-9) across 38 organic active materials and applies it to four ortho-substituted hydroquinone posolytes studied by [Yang et al. (2021)](https://doi.org/10.3390/molecules26133823) and [Navarro García et al. (2021)](https://doi.org/10.3390/molecules26051203).

---

## Key Chemistry & Audit Findings

- **The phosphate comparison shows a 99.1% reduction in modeled cost:** At an assumed `1.40 V`, using phosphate PS fade rates and saturation solubilities gives about `$1.20M/kWh` for `O-3` and `$10,735–$10,798/kWh` for `asym-O-5`. In the separate sulfate PS scenario, `asym-O-5` would need **42–56×** slower fade to match the published ammonium ferrocyanide benchmark (`$651.7/kWh`) and **84–168×** slower fade to match Future Case VRFB (`$528.5/kWh`). The latter requires a maximum of about `0.00274–0.00548%/day` across the two maintenance conventions.
- **Cycling-rate sensitivity for the same formulation:** Yang et al. measured `0.46%/day` under PS and `1.80%/day` under GS for `0.10 M` `asym-O-5` in `1 M H2SO4`. Substituting the GS fade rate in the `0.44 M` saturation scenario increases modeled cost by **$33,363/kWh**.
- **Polymorph solubility loss amplifies modeled replacement cost:** Reducing `O-5` sulfate solubility from `0.24 M` to `0.06 M` adds **$50,263–$50,333/kWh** at the fixed PS fade rate and **$70,880–$70,950/kWh** at the fixed GS fade rate across the two maintenance conventions.
- **Published benchmark reconciliation:** Published components sum to `676.7`, `370.0`, and `3785.5 USD/kWh` for Present Case VRFB, Future Case `1,6-DPAP`, and Future Case `TEMPTMA`, respectively, at one decimal place. Independent totals are `$840.62`, `$459.63`, and `$3826.87/kWh`. The audit also documents the approximate `+0.60 V` reference offset for two points in Emmel Fig. 2 and tests the effect of a single-stack replacement assumption against the system-wide replacement specified in Eq. S40.

---

## Quickstart (Local App, Tests & Notebook Reproduction)

Requires **Python 3.11+**.

```bash
# 1. Install dependencies
python -m pip install -r requirements.txt

# 2. Run unit & benchmark regression tests
pytest -v

# 3. Launch the interactive Streamlit calculator locally
streamlit run app.py

# 4. (Optional) Headlessly re-execute the original audit notebook and rebuild figures/
python reproduce.py
```

---

## What the Interactive Web Tool (`app.py`) Does

1. **Single-Molecule Cost & Protocol Audit:** Load literature presets (`O-1`, `O-3`, `O-5`, `asym-O-5`) or enter custom $(M_w, n, U_d, c, \text{fade})$ parameters to inspect full Eq. 1/5/6 cost breakdowns, 20-year electrolyte renewal counts, inverse fade targets, and one-click PS/GS & polymorph stress tests. Export any audit run as CSV.
2. **42-Molecule Benchmark Explorer & Inverse Contours:** Visualize your candidate overlaid against all 38 Emmel molecules and 4 Yang hydroquinones on log-log axes, plus an interactive 2D voltage-vs.-fade contour map with carbon-HER OCV guide (`1.45 V`).
3. **Methodology, Audit Notes & Reproducibility:** Inspect the live benchmark reconciliation table, assumption boundaries, and full primary literature citations.

---

## Repository Structure

```text
aorfb-cost-audit/
├── app.py                              # 3-tab Streamlit interactive calculator
├── aorfb_audit/
│   ├── __init__.py                     # Public Python API
│   ├── model.py                        # Eqs. 1, 2-3, 5, 6, bisection targets & audits
│   └── data_loader.py                  # CSV loaders & physical-bounds validation
├── data/
│   ├── emmel_38_molecules_audit.csv    # Transcribed Emmel Tables S8/S10 + digitized Fig. S9/S10
│   └── yang_tables_1_to_3.csv          # Transcribed Yang Tables 1-3 cycling & capacity data
├── tests/
│   └── test_model.py                   # Unit & benchmark regression test suite
├── figures/                            # Publication-style audit figures rebuilt by reproduce.py
├── flow_battery_benchmark.ipynb        # Original self-contained audit notebook (Sections 0-5)
├── reproduce.py                        # Headless notebook runner
├── pyproject.toml                      # Package metadata
├── requirements.txt                    # Python dependencies (local + Streamlit Cloud)
├── .zenodo.json                        # Zenodo archival metadata for DOI minting
└── LICENSE                             # MIT License
```

---

## Scope & Modeling Caveats

- Scenario calculations assume a 4 h system, neutral-molecule active pricing (`$3.48/kg`), an idealized counter-electrolyte (`100 g/mol`, `1 e-`, `2 M`), and an illustrative WP1 stack count (`430.17–438.64`).
- The model extrapolates linear fade from dilute symmetric cells (`0.05–0.10 M`) to saturation and 20-year operation. Each renewal is charged on the full Eq. 5 cost basis, including tanks, following the printed Eq. S38 expression. These costs quantify sensitivity to scenario assumptions.
- The `1.45 V` open-circuit guide is an estimate for `asym-O-5` against the pH 0 carbon hydrogen-evolution bound. The app marks discharge voltages above this guide as exploratory.

---

## Primary Sources & Licensing

- **Yang et al.**, *Molecules* **26**, 3823 (2021). [DOI: 10.3390/molecules26133823](https://doi.org/10.3390/molecules26133823)
- **Navarro García et al.**, *Molecules* **26**, 1203 (2021). [DOI: 10.3390/molecules26051203](https://doi.org/10.3390/molecules26051203)
- **Emmel et al.**, *Nature Communications* **14**, 6672 (2023). [DOI: 10.1038/s41467-023-42450-9](https://doi.org/10.1038/s41467-023-42450-9)

Original code is released under the **MIT License**. Transcribed or adapted source data retain attribution to the cited studies under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
