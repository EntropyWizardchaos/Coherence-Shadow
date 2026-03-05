# Methodology

## Data Source

**SPARC** (Spitzer Photometry and Accurate Rotation Curves)  
Lelli, McGaugh, & Schombert (2016)  
http://astroweb.cwru.edu/SPARC/

175 disk galaxies with:
- High-quality rotation curves (HI + Hα)
- Spitzer 3.6μm photometry (stellar mass tracer)
- Mass decomposition into gas, disk, bulge components

---

## Coherence Metrics

### 1. Kinematic Coherence (C_kin)

Measures rotation curve smoothness.

```python
def kinematic_coherence(R, V):
    """
    R: galactocentric radius (kpc)
    V: observed rotation velocity (km/s)
    Returns: C_kin in [0, 1], higher = smoother
    """
    dVdR = np.gradient(V, R)
    d2VdR2 = np.gradient(dVdR, R)
    
    roughness = np.mean(np.abs(d2VdR2))
    roughness_normalized = roughness * (R.max() / V.max())
    
    return np.exp(-roughness_normalized)
```

**Interpretation:** High C_kin means the rotation curve has low second-derivative magnitude — it rises and flattens smoothly without wiggles.

### 2. Morphological Coherence (C_conc)

Measures light concentration — how organized the stellar distribution is.

```python
def concentration_index(R, SB):
    """
    R: galactocentric radius (kpc)
    SB: disk surface brightness (L_sun/pc^2)
    Returns: C_conc = R_50 / R_90
    """
    # Sort by radius
    order = np.argsort(R)
    R_s, SB_s = R[order], SB[order]
    
    # Cumulative light
    dR = np.gradient(R_s)
    dL = SB_s * 2 * np.pi * R_s * np.abs(dR)
    cumL = np.cumsum(dL)
    cumL_norm = cumL / cumL[-1]
    
    # Find R_50 and R_90
    R50 = np.interp(0.5, cumL_norm, R_s)
    R90 = np.interp(0.9, cumL_norm, R_s)
    
    return R50 / R90
```

**Interpretation:** High C_conc means more light is concentrated toward the center — the disk has a well-organized, steeply declining profile.

---

## Dark Matter Fraction

Standard mass decomposition assuming Newtonian gravity:

```python
def dark_matter_fraction(V_obs, V_gas, V_disk, V_bul, Y_star=0.5):
    """
    V_obs: observed rotation velocity
    V_gas: gas contribution (includes 1.33× for He)
    V_disk: disk contribution (for M/L = 1)
    V_bul: bulge contribution (for M/L = 1)
    Y_star: stellar mass-to-light ratio at 3.6μm
    
    Returns: f_dm = 1 - V_bar²/V_obs²
    """
    V_bar_sq = V_gas**2 + Y_star * (V_disk**2 + V_bul**2)
    V_obs_sq = V_obs**2
    
    f_dm = 1 - V_bar_sq / V_obs_sq
    return np.clip(f_dm, 0, 1)
```

**Y_star = 0.5** is the standard value for 3.6μm stellar mass-to-light ratio (Lelli et al. 2016).

**Per-galaxy f_dm:** We use the mean f_dm at outer radii (last 3 data points), where dark matter dominates.

---

## Statistical Methods

### Simple Linear Regression

```python
slope, intercept, r, p, se = scipy.stats.linregress(C, f_dm)
```

### Multiple Regression

For controlled analyses:

```python
X = np.column_stack([np.ones(N), C, confounders...])
beta = np.linalg.lstsq(X, f_dm, rcond=None)[0]
# Standard errors from covariance matrix
```

### High vs Low Comparison

Split at median, compare means via t-test:

```python
high = galaxies[C >= median(C)]
low = galaxies[C < median(C)]
t, p = scipy.stats.ttest_ind(high['f_dm'], low['f_dm'])
```

---

## Sample Selection

- All 175 SPARC galaxies processed
- Excluded: galaxies with < 8 data points in rotation curve
- Final sample: 143 galaxies with valid coherence and f_dm measurements

---

## Reproducibility

All results in this repository are generated from the raw SPARC Rotmod_LTG rotation curve files using the scripts in `analysis/`. To reproduce:

1. Download SPARC data from http://astroweb.cwru.edu/SPARC/
2. Run `python analysis/sparc_coherence.py --sparc-dir /path/to/Rotmod_LTG/`
3. The script computes C_kin, C_conc, f_dm, f_gas, and radial bin values for each galaxy
4. All statistical tests (regressions, splits, controls) are computed and printed

The included `data/sparc_results.csv` contains pre-computed values for quick verification. The concentration test can also be run standalone: `python analysis/sparc_concentration.py --csv data/sparc_results.csv`

---

## Predictions Tested

Prior to analysis, the framework predicted:

1. **Sign:** β(C → f_dm) < 0 (more coherent → less dark matter)
2. **Magnitude:** |β| ~ 0.1
3. **Radial structure:** Effect peaks at intermediate radii
4. **Independence:** Should hold after controlling for mass, morphology

All predictions were specified before running the analysis on SPARC data.
