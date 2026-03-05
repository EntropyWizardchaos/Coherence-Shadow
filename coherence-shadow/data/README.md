# Data Notes

## sparc_results.csv

This CSV contains pre-computed per-galaxy values from the SPARC dataset.

### Columns
- `galaxy` — SPARC galaxy name
- `C_kin` — Kinematic coherence (rotation curve smoothness, 0-1)
- `C_conc` — Morphological coherence (concentration index R50/R90)
- `f_dm` — Dark matter fraction (outer radii, Y* = 0.5)
- `Vmax` — Maximum observed rotation velocity (km/s)
- `Rmax` — Maximum radius (kpc)
- `N_points` — Number of data points in rotation curve

### Key Result

The concentration test can be verified directly from this CSV:

```python
import pandas as pd
from scipy import stats

df = pd.read_csv("sparc_results.csv")
slope, intercept, r, p, se = stats.linregress(df["C_conc"], df["f_dm"])
print(f"C_conc → f_dm: β = {slope:.4f}, p = {p:.6f}")
# Expected: β ≈ -0.354, p ≈ 0.005
```

### To Regenerate from Raw SPARC Data

Download SPARC from http://astroweb.cwru.edu/SPARC/ and run:

```bash
python analysis/sparc_coherence.py --sparc-dir /path/to/Rotmod_LTG/ --output data/sparc_results.csv
```
