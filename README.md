# AORFB Cost & Cycling-Protocol Audit (`aorfb-cost-audit`)

> **Interactive Web App:** `[Paste your Streamlit Community Cloud URL after deployment]`  
> **Archival DOI:** `[Paste your Zenodo DOI badge after creating GitHub Release v1.0.0]`

An open-source electrochemistry modeling package, reproducible audit notebook, and browser-accessible **Streamlit calculator** that asks a practical question for aqueous organic redox flow battery (AORFB) design: **how much must positive-electrolyte capacity fade improve for substituted hydroquinones to approach inorganic cost benchmarks, and how severely do inconsistent laboratory protocols (PS vs. GS, concentration, polymorph aging) skew system-level rankings?**

This project independently audits the techno-economic model of [Emmel et al. (2023)](https://doi.org/10.1038/s41467-023-42450-9) across 38 organic active materials and applies it to four ortho-substituted hydroquinone posolytes studied by [Yang et al. (2021)](https://doi.org/10.3390/molecules26133823) and [Navarro García et al. (2021)](https://doi.org/10.3390/molecules26051203).

---

## Key Chemistry & Audit Findings

- **Ortho-blocking cuts modeled cost by 99.1%, yet still trails inorganic benchmarks by >$10,000/kWh:** Under common phosphate PS conditions (`1.40 V`), `asym-O-5` lowers modeled capital cost from `$1.20M/kWh` (`O-3`) to `$10,735–$10,798/kWh`, but still requires **42–56x** slower fade to match ammonium ferrocyanide (`$651.7/kWh`) and **84–168x** slower fade (`≤ 0.00274–0.00548%/day`) to match Future Case VRFB (`$528.5/kWh`).
- **Cycling protocol bias alters cost by >$33,000/kWh on the same molecule:** For matched-formulation `asym-O-5` (`0.10 M` in `1 M H2SO4`), substituting galvanostatic (`1.8%/day`) for potentiostatic (`0.46%/day`) fade increases modeled capital cost by **+$33,363/kWh**.
- **Polymorph solubility loss amplifies replacement cost:** Keeping `O-5` fade fixed while reducing sulfate solubility from fresh (`0.24 M`) to aged polymorph (`0.06 M`) adds **+$50,333/kWh** (PS) or **+$70,950/kWh** (GS).
- **Published benchmark reconciliation:** Published component sums recover `676.7`, `370.0`, and `3785.5 USD/kWh` for Present Case VRFB, Future Case `1,6-DPAP`, and Future Case `TEMPTMA`. Independent recalculation documents both the `+0.60 V` carbon-HER reference offset in Emmel Fig. 2 and the single-stack vs. full-stack Eq. S40 maintenance diagnostic.

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
- Extrapolating dilute symmetric-cell laboratory fade (`0.05–0.10 M`) to saturation solubility and 20-year operation is a scenario sensitivity audit rather than a field lifetime forecast.
- Voltages above `1.45 V` for acidic hydroquinone posolytes exceed the approximate carbon hydrogen-evolution OCV guide and are marked as exploratory.

---

## Primary Sources & Licensing

- **Yang et al.**, *Molecules* **26**, 3823 (2021). [DOI: 10.3390/molecules26133823](https://doi.org/10.3390/molecules26133823)
- **Navarro García et al.**, *Molecules* **26**, 1203 (2021). [DOI: 10.3390/molecules26051203](https://doi.org/10.3390/molecules26051203)
- **Emmel et al.**, *Nature Communications* **14**, 6672 (2023). [DOI: 10.1038/s41467-023-42450-9](https://doi.org/10.1038/s41467-023-42450-9)

Original code is released under the **MIT License**. Transcribed or adapted source data retain attribution to the cited studies under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
