# Coherence-Shadow

**Galactic coherence correlates with dark matter fraction**

**Author:** Harley Robinson  
**Date:** March 2026  
**Data:** SPARC (Spitzer Photometry and Accurate Rotation Curves), N=175 galaxies

---

## Summary

We find that galactic coherence — measured two independent ways — correlates negatively with inferred dark matter fraction.

| Coherence Measure | β | p-value |
|-------------------|---------|---------|
| Kinematic (rotation curve smoothness) | -0.068 | 0.139 |
| **Morphological (light concentration)** | **-0.354** | **0.005** |

More ordered galaxies have less dark matter. The two measures are uncorrelated with each other (r = 0.13), yet both point the same direction.

---

## The Finding

**Concentration index** (R₅₀/R₉₀, ratio of radii containing 50% vs 90% of disk light) predicts dark matter fraction at p = 0.005.

This correlation is:
- Independent of stellar mass (r = -0.02)
- Independent of Hubble morphology (r = -0.04)
- Independent of kinematic coherence (r = 0.13)

High-concentration galaxies have **7% less inferred dark matter** than low-concentration galaxies (p = 0.007).

---

## Replication

### Requirements
- Python 3.8+
- numpy, pandas, scipy
- SPARC data (download from [SPARC database](http://astroweb.cwru.edu/SPARC/))

### Steps

```bash
# Clone this repository
git clone https://github.com/EntropyWizardchaos/coherence-shadow.git
cd coherence-shadow

# Install dependencies
pip install -r analysis/requirements.txt

# Download SPARC Rotmod_LTG data from http://astroweb.cwru.edu/SPARC/
# Extract to a local directory, e.g., ./sparc_data/Rotmod_LTG/

# Run full analysis (produces sparc_results.csv with all columns)
python analysis/sparc_coherence.py --sparc-dir /path/to/Rotmod_LTG/ --output data/sparc_results.csv

# Optional: provide SPARC summary table for Hubble type controls
python analysis/sparc_coherence.py --sparc-dir /path/to/Rotmod_LTG/ \
    --sparc-table /path/to/SPARC_Lelli2016c.mrt \
    --output data/sparc_results.csv

# Run standalone concentration test
python analysis/sparc_concentration.py --sparc-dir /path/to/Rotmod_LTG/

# Or run concentration test on pre-computed CSV
python analysis/sparc_concentration.py --csv data/sparc_results.csv
```

### What You Should Get

If the result replicates:
- `C_conc → f_dm`: β ≈ -0.35, p < 0.01
- High vs Low concentration Δf_dm ≈ -0.07, p < 0.01
- `C_kin → f_dm`: β ≈ -0.07, p ≈ 0.14 (correct sign, marginal)

If it doesn't replicate, please open an issue.

### Note on Included Data

`data/sparc_results.csv` contains pre-computed results for verification. To reproduce from scratch, you must download the SPARC Rotmod_LTG files and run the analysis scripts. The scripts compute all values from the raw rotation curve data — nothing is hard-coded.

---

## Key Results

### 1. Concentration Predicts Dark Matter Fraction

```
f_dm ~ C_conc
β = -0.354 ± 0.125
p = 0.005
```

### 2. Effect Is Not Explained by Confounders

| Control Variable | Correlation with C_conc |
|------------------|------------------------|
| Stellar mass     | r = -0.02 (p = 0.77)  |
| Hubble type      | r = -0.04 (p = 0.66)  |
| Kinematic coherence | r = 0.13 (p = 0.13) |

### 3. Two Independent Measures Agree

Kinematic and morphological coherence are weakly correlated, yet both predict lower dark matter fraction:

```
Kinematic:     β = -0.068 (p = 0.14)
Morphological: β = -0.354 (p = 0.005)
```

### 4. Gas Fraction Split

The effect is **stronger in low-gas galaxies** (β = -0.10) than high-gas galaxies (β = -0.004). This is opposite to what baryonic feedback models predict.

---

## Files

```
coherence-shadow/
├── README.md                       # This file
├── LICENSE                         # MIT License
├── methodology.md                  # How coherence and f_dm are computed
├── results.md                      # Full results with tables and interpretation
├── analysis/
│   ├── sparc_coherence.py          # Full pipeline: both coherence measures,
│   │                               #   radial bins, gas splits, regressions
│   ├── sparc_concentration.py      # Standalone morphological coherence test
│   └── requirements.txt            # Python dependencies
├── data/
│   └── sparc_results.csv           # Pre-computed per-galaxy values
└── theory/
    └── cmd_framework.md            # Theoretical context (optional reading)
```

---

## Theoretical Context

These tests were motivated by a framework predicting that coherence (organizational order) should correlate with reduced dark matter fraction. See `theory/cmd_framework.md` for background.

However, **the empirical result stands independent of theory**. The correlation either replicates or it doesn't.

---

## Citation

If you use this work, please cite:

```
Robinson, H. (2026). Coherence-Shadow: Galactic coherence correlates with dark matter fraction.
GitHub repository: https://github.com/EntropyWizardchaos/coherence-shadow
```

---

## Contact

Questions, replications, or critiques welcome via GitHub issues.

---

*"The test is the truth. Everything else is conversation."*
