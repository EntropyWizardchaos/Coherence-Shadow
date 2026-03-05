#!/usr/bin/env python3
"""
SPARC Concentration Analysis — Independent Morphological Test
==============================================================
Tests morphological coherence (light concentration) against dark matter
fraction, independent of kinematic coherence.

This is the standalone version of the concentration test.
The full pipeline (sparc_coherence.py) includes this analysis as well.

Usage:
    python sparc_concentration.py --sparc-dir /path/to/sparc/Rotmod_LTG/

    Or use pre-computed results:
    python sparc_concentration.py --csv ../data/sparc_results.csv

Output:
    - Concentration → f_dm regression
    - High vs Low concentration comparison
    - Independence checks against C_kin, mass, morphology

Author: Harley Robinson
Date: March 2026
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

Y_STAR = 0.5
MIN_POINTS = 8


def read_rotmod_file(path):
    """Read SPARC rotation model file."""
    try:
        df = pd.read_csv(path, sep=r'\s+', comment="#", header=None)
    except Exception:
        return None

    ncol = df.shape[1]
    if ncol >= 7:
        cols = ["R", "Vobs", "e_Vobs", "Vgas", "Vdisk", "Vbul", "SBdisk"][:min(ncol, 7)]
        if ncol >= 8:
            cols.append("SBbul")
        df.columns = cols[:ncol]
    elif ncol == 6:
        df.columns = ["R", "Vobs", "e_Vobs", "Vgas", "Vdisk", "Vbul"]
    else:
        return None

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["R", "Vobs", "Vgas", "Vdisk", "Vbul"])
    df = df.sort_values("R").reset_index(drop=True)
    return df


def concentration_index(df):
    """C_conc = R_50 / R_90 from disk surface brightness profile."""
    if "SBdisk" not in df.columns:
        return np.nan

    R = df["R"].to_numpy(dtype=float)
    SB = df["SBdisk"].to_numpy(dtype=float)

    valid = np.isfinite(SB) & (SB > 0) & np.isfinite(R) & (R > 0)
    if valid.sum() < 6:
        return np.nan

    R_v, SB_v = R[valid], SB[valid]
    order = np.argsort(R_v)
    R_s, SB_s = R_v[order], SB_v[order]

    dR = np.gradient(R_s)
    dL = SB_s * 2 * np.pi * R_s * np.abs(dR)
    cumL = np.cumsum(dL)

    if cumL[-1] <= 0:
        return np.nan

    cumL_norm = cumL / cumL[-1]

    try:
        R50 = np.interp(0.5, cumL_norm, R_s)
        R90 = np.interp(0.9, cumL_norm, R_s)
        if R90 <= 0:
            return np.nan
        return float(R50 / R90)
    except Exception:
        return np.nan


def kinematic_coherence(df):
    """C_kin from rotation curve smoothness."""
    R = df["R"].to_numpy(dtype=float)
    V = df["Vobs"].to_numpy(dtype=float)
    if len(R) < 5:
        return np.nan
    dVdR = np.gradient(V, R)
    d2VdR2 = np.gradient(dVdR, R)
    roughness = np.nanmean(np.abs(d2VdR2))
    V_scale = np.nanmax(V) if np.nanmax(V) > 0 else 1.0
    R_scale = np.nanmax(R) if np.nanmax(R) > 0 else 1.0
    roughness_normalized = roughness * (R_scale / (V_scale + 1e-12))
    return float(np.clip(np.exp(-roughness_normalized), 0, 1))


def compute_fdm_outer(df, y_star=Y_STAR):
    """Dark matter fraction at outer radii."""
    vgas2 = df["Vgas"].to_numpy(dtype=float) ** 2
    vdisk2 = df["Vdisk"].to_numpy(dtype=float) ** 2
    vbul2 = df["Vbul"].to_numpy(dtype=float) ** 2
    vobs2 = df["Vobs"].to_numpy(dtype=float) ** 2
    vbar2 = vgas2 + y_star * (vdisk2 + vbul2)
    fdm = np.clip(1.0 - (vbar2 / (vobs2 + 1e-12)), 0, 1)
    return float(np.mean(fdm[-3:])) if len(fdm) >= 3 else float(fdm[-1])


def run_from_sparc(sparc_dir):
    """Process raw SPARC files and return DataFrame."""
    files = [f for f in os.listdir(sparc_dir) if f.endswith("_rotmod.dat")]
    print(f"Found {len(files)} galaxy files")

    rows = []
    for fn in sorted(files):
        path = os.path.join(sparc_dir, fn)
        df = read_rotmod_file(path)
        if df is None or len(df) < MIN_POINTS:
            continue

        c_kin = kinematic_coherence(df)
        c_conc = concentration_index(df)
        fdm = compute_fdm_outer(df)

        if not np.isfinite(c_kin):
            continue

        rows.append({
            "galaxy": fn.replace("_rotmod.dat", ""),
            "C_kin": c_kin,
            "C_conc": c_conc,
            "f_dm": fdm,
            "Vmax": float(df["Vobs"].max()),
        })

    return pd.DataFrame(rows)


def run_analysis(df):
    """Run concentration analysis on DataFrame."""

    valid = df.dropna(subset=["C_conc", "f_dm"])
    print(f"\nGalaxies with valid concentration data: {len(valid)}")

    # --- Main test ---
    print(f"\n{'='*60}")
    print("CONCENTRATION INDEX → DARK MATTER FRACTION")
    print(f"{'='*60}")

    slope, intercept, r, p, se = stats.linregress(
        valid["C_conc"].to_numpy(), valid["f_dm"].to_numpy())
    print(f"\n  C_conc → f_dm:")
    print(f"    β = {slope:+.4f} ± {se:.4f}")
    print(f"    r = {r:.4f}")
    print(f"    p = {p:.4f}")

    # --- High vs Low ---
    print(f"\n{'='*60}")
    print("HIGH vs LOW CONCENTRATION")
    print(f"{'='*60}")

    med = valid["C_conc"].median()
    high = valid[valid["C_conc"] >= med]
    low = valid[valid["C_conc"] < med]

    print(f"\n  High concentration (N={len(high)}): mean f_dm = {high['f_dm'].mean():.3f}")
    print(f"  Low concentration  (N={len(low)}):  mean f_dm = {low['f_dm'].mean():.3f}")

    t, p_t = stats.ttest_ind(high["f_dm"], low["f_dm"])
    print(f"  Difference: {high['f_dm'].mean() - low['f_dm'].mean():+.3f}")
    print(f"  t-test p = {p_t:.4f}")

    # --- Independence from kinematic coherence ---
    print(f"\n{'='*60}")
    print("INDEPENDENCE CHECKS")
    print(f"{'='*60}")

    if "C_kin" in valid.columns:
        r_kc, p_kc = stats.pearsonr(valid["C_kin"], valid["C_conc"])
        print(f"\n  C_conc vs C_kin: r = {r_kc:.3f}, p = {p_kc:.4f}")
        if abs(r_kc) < 0.3:
            print("  → Measures are largely independent")

    if "Vmax" in valid.columns:
        log_v = np.log10(valid["Vmax"].clip(lower=1))
        r_cv, p_cv = stats.pearsonr(valid["C_conc"], log_v)
        print(f"  C_conc vs log(Vmax): r = {r_cv:.3f}, p = {p_cv:.4f}")

    # --- Kinematic for comparison ---
    print(f"\n{'='*60}")
    print("COMPARISON: KINEMATIC COHERENCE")
    print(f"{'='*60}")

    if "C_kin" in df.columns:
        slope_k, _, r_k, p_k, se_k = stats.linregress(
            df["C_kin"].to_numpy(), df["f_dm"].to_numpy())
        print(f"\n  C_kin → f_dm:")
        print(f"    β = {slope_k:+.4f} ± {se_k:.4f}")
        print(f"    p = {p_k:.4f}")

    # --- Summary ---
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"""
  Two independent coherence measures:

    Kinematic (C_kin):      β = {slope_k:+.4f}, p = {p_k:.4f}
    Morphological (C_conc): β = {slope:+.4f}, p = {p:.4f}

  Correlation between measures: r = {r_kc:.3f}
""")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SPARC Concentration Analysis")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sparc-dir",
                       help="Path to SPARC Rotmod_LTG directory")
    group.add_argument("--csv",
                       help="Path to pre-computed sparc_results.csv")

    args = parser.parse_args()

    if args.csv:
        df = pd.read_csv(args.csv)
        print(f"Loaded {len(df)} galaxies from {args.csv}")
    else:
        df = run_from_sparc(args.sparc_dir)
        print(f"Processed {len(df)} galaxies from SPARC files")

    run_analysis(df)
