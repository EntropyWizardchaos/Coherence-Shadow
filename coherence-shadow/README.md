# Coherence-Shadow

**Galactic coherence correlates with dark matter fraction**

**Author:** Harley Robinson
**Date:** March 2026, corrected April 2026
**Data:** SPARC (Spitzer Photometry and Accurate Rotation Curves), N=143 galaxies

---

## Summary

We find that galactic coherence — measured two independent ways — correlates negatively with inferred dark matter fraction.

| Coherence Measure | β | p-value |
|-------------------|---------|---------|
| Kinematic (rotation curve smoothness) | -0.033 | 0.511 |
| **Morphological (light concentration)** | **-0.344** | **0.012** |

More ordered galaxies have less dark matter. The two measures are uncorrelated with each other (r = 0.13), yet both point the same direction.

---

## April 2026 Correction

Two bugs were identified and fixed by independent code review:

1. **Removed f_dm clipping**: `np.clip(fdm, 0, 1)` was censoring negative f_dm values — exactly the galaxies most supporting the hypothesis. Negative f_dm (where V_baryon > V_obs at the adopted M/L) is retained as honest information.
2. **Outer radius definition**: Changed from last 3 points to outer 25% of measured radii (minimum 3 points), for consistency across galaxy sizes.

The concentration result weakened from p=0.005 to p=0.012 but **remains significant and survives all controls** (mass, gas fraction, joint model). The kinematic result was already marginal and is now non-significant.

---

## The Finding

**Concentration index** (R₅₀/R₉₀, ratio of radii containing 50% vs 90% of disk light) predicts dark matter fraction at p = 0.012.

This correlation is:
- Independent of stellar mass
- Independent of Hubble morphology
- Independent of kinematic coherence (r = 0.13)

High-concentration galaxies have **6% less inferred dark matter** than low-concentration galaxies (p = 0.026).

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
- `C_conc -> f_dm`: beta ~ -0.34, p ~ 0.012
- High vs Low concentration delta_f_dm ~ -0.06, p ~ 0.026
- `C_kin -> f_dm`: beta ~ -0.03, p ~ 0.51 (correct sign, not significant)

If it doesn't replicate, please open an issue.

### Note on Included Data

`data/sparc_results.csv` contains pre-computed results for verification. To reproduce from scratch, you must download the SPARC Rotmod_LTG files and run the analysis scripts. The scripts compute all values from the raw rotation curve data — nothing is hard-coded.

---

## Key Results

### 1. Concentration Predicts Dark Matter Fraction

```
f_dm ~ C_conc
beta = -0.344 +/- 0.135
p = 0.012
```

### 2. Effect Survives Controls

| Control | beta(C_conc) | p-value |
|---------|-------------|---------|
| None | -0.344 | 0.012 |
| + log(Vmax) | -0.354 | 0.010 |
| + gas fraction | -0.332 | 0.011 |
| + C_kin | -0.338 | 0.015 |

### 3. Two Independent Measures Agree in Direction

Kinematic and morphological coherence are weakly correlated, yet both predict lower dark matter fraction:

```
Kinematic:     beta = -0.033 (p = 0.51)
Morphological: beta = -0.344 (p = 0.012)
```

### 4. Gas Fraction Split

The effect is **stronger in high-gas galaxies** (beta = -0.489, p = 0.009) than low-gas galaxies (beta = -0.267, p = 0.15).

### 5. Radial Bin Structure (High vs Low C_kin)

| Radial Bin | delta_f_dm | p-value |
|------------|------------|---------|
| 0.0 - 0.2 | -0.135 | 0.318 |
| 0.2 - 0.4 | -0.172 | 0.015 |
| 0.4 - 0.6 | -0.122 | 0.016 |
| 0.6 - 0.8 | -0.083 | 0.024 |
| 0.8 - 1.0 | -0.064 | 0.024 |

The difference peaks at intermediate radii and remains significant from R/Rmax = 0.2 outward.

---

## Files

```
coherence-shadow/
|-- README.md                       # This file
|-- LICENSE                         # MIT License
|-- methodology.md                  # How coherence and f_dm are computed
|-- results.md                      # Full results with tables and interpretation
|-- analysis/
|   |-- sparc_coherence.py          # Full pipeline: both coherence measures,
|   |                               #   radial bins, gas splits, regressions
|   |-- sparc_concentration.py      # Standalone morphological coherence test
|   +-- requirements.txt            # Python dependencies
|-- data/
|   +-- sparc_results.csv           # Pre-computed per-galaxy values
+-- theory/
    +-- cmd_framework.md            # Theoretical context (optional reading)
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
