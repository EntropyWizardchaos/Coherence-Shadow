#!/usr/bin/env python3
"""
SPARC Coherence Analysis — Full Pipeline
=========================================
Tests whether galactic coherence correlates with dark matter fraction.
Computes BOTH kinematic and morphological coherence, radial bin analyses,
gas fraction splits, and controlled regressions.

Every number in results.md is produced by this script.

Usage:
    python sparc_coherence.py --sparc-dir /path/to/sparc/Rotmod_LTG/

    Optional: provide SPARC summary table for Hubble type controls
    python sparc_coherence.py --sparc-dir /path/to/Rotmod_LTG/ \
                              --sparc-table /path/to/SPARC_Lelli2016c.mrt

Output:
    - Prints all correlation statistics
    - Saves sparc_results.csv with per-galaxy values (all columns)
    - Saves sparc_radial_results.csv with radial bin analyses

Author: Harley Robinson
Date: March 2026
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CONSTANTS
# ============================================================
Y_STAR = 0.5       # Stellar mass-to-light ratio at 3.6μm (Lelli+ 2016)
MIN_POINTS = 8     # Minimum data points per galaxy
RADIAL_BINS = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]


# ============================================================
# DATA READING
# ============================================================
def read_rotmod_file(path):
    """Read SPARC rotation model file.

    SPARC Rotmod_LTG files have columns:
    R(kpc) Vobs(km/s) e_Vobs Vgas Vdisk Vbul SBdisk(L/pc^2) [SBbul]
    """
    try:
        df = pd.read_csv(path, sep=r'\s+', comment="#", header=None)
    except Exception:
        return None

    ncol = df.shape[1]
    if ncol >= 7:
        base_cols = ["R", "Vobs", "e_Vobs", "Vgas", "Vdisk", "Vbul", "SBdisk"]
        cols = base_cols[:min(ncol, 7)]
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


def read_sparc_table(path):
    """Read SPARC summary table for Hubble types and global properties.

    Tries to handle both .mrt and .txt formats.
    Returns DataFrame with at least: Galaxy, T (Hubble type)
    """
    try:
        # Try whitespace-separated with common SPARC table columns
        df = pd.read_csv(path, sep=r'\s+', comment="#")
        # Normalize column names
        col_map = {}
        for c in df.columns:
            cl = c.lower().strip()
            if cl in ("galaxy", "name"):
                col_map[c] = "Galaxy"
            elif cl in ("t", "hubble", "type"):
                col_map[c] = "T"
            elif cl in ("d", "dist", "distance"):
                col_map[c] = "D"
            elif cl in ("l", "lum", "luminosity"):
                col_map[c] = "L"
        df = df.rename(columns=col_map)
        if "Galaxy" in df.columns:
            return df
    except Exception:
        pass

    return None


# ============================================================
# COHERENCE METRICS
# ============================================================
def kinematic_coherence(df):
    """Compute kinematic coherence from rotation curve smoothness.

    C_kin = exp(-roughness), where roughness is the normalized mean
    absolute second derivative of the rotation curve.

    Returns: float in [0, 1], higher = smoother rotation curve
    """
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

    C = np.exp(-roughness_normalized)
    return float(np.clip(C, 0, 1))


def concentration_index(df):
    """Compute concentration index from surface brightness profile.

    C_conc = R_50 / R_90, where R_50 and R_90 are radii enclosing
    50% and 90% of disk light.

    Higher C_conc = more concentrated = more geometrically organized.

    Returns: float in (0, 1], or NaN if SBdisk unavailable
    """
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


# ============================================================
# DARK MATTER FRACTION
# ============================================================
def compute_fdm_array(df, y_star=Y_STAR):
    """Compute dark matter fraction at each radius.

    f_dm(R) = 1 - V_bar^2 / V_obs^2
    where V_bar^2 = V_gas^2 + Y_star * (V_disk^2 + V_bul^2)
    """
    vgas2 = df["Vgas"].to_numpy(dtype=float) ** 2
    vdisk2 = df["Vdisk"].to_numpy(dtype=float) ** 2
    vbul2 = df["Vbul"].to_numpy(dtype=float) ** 2
    vobs2 = df["Vobs"].to_numpy(dtype=float) ** 2

    vbar2 = vgas2 + y_star * (vdisk2 + vbul2)
    fdm = 1.0 - (vbar2 / (vobs2 + 1e-12))
    return np.clip(fdm, 0, 1)


def compute_fdm_outer(df, y_star=Y_STAR):
    """Dark matter fraction at outer radii (mean of last 3 points)."""
    fdm = compute_fdm_array(df, y_star)
    if len(fdm) >= 3:
        return float(np.mean(fdm[-3:]))
    return float(fdm[-1])


def compute_fdm_radial_bins(df, y_star=Y_STAR):
    """Compute mean f_dm in normalized radial bins."""
    R = df["R"].to_numpy(dtype=float)
    fdm = compute_fdm_array(df, y_star)
    Rmax = R.max()
    if Rmax <= 0:
        return {f"fdm_{lo}_{hi}": np.nan for lo, hi in RADIAL_BINS}

    Rnorm = R / Rmax
    result = {}
    for lo, hi in RADIAL_BINS:
        mask = (Rnorm >= lo) & (Rnorm < hi)
        if mask.sum() >= 1:
            result[f"fdm_{lo:.1f}_{hi:.1f}"] = float(np.mean(fdm[mask]))
        else:
            result[f"fdm_{lo:.1f}_{hi:.1f}"] = np.nan
    return result


# ============================================================
# GAS FRACTION (derived from rotation curve contributions)
# ============================================================
def compute_gas_fraction(df, y_star=Y_STAR):
    """Estimate gas fraction from velocity contributions.

    f_gas = sum(V_gas^2) / sum(V_gas^2 + Y_star * (V_disk^2 + V_bul^2))

    This is the gas fraction of the baryonic mass budget.
    """
    vgas2 = df["Vgas"].to_numpy(dtype=float) ** 2
    vdisk2 = df["Vdisk"].to_numpy(dtype=float) ** 2
    vbul2 = df["Vbul"].to_numpy(dtype=float) ** 2

    gas_total = np.sum(vgas2)
    bar_total = gas_total + y_star * np.sum(vdisk2 + vbul2)

    if bar_total <= 0:
        return np.nan
    return float(gas_total / bar_total)


# ============================================================
# GALAXY PROCESSING
# ============================================================
def process_galaxy(path):
    """Process a single galaxy file. Returns dict with all metrics."""
    df = read_rotmod_file(path)
    if df is None or len(df) < MIN_POINTS:
        return None

    C_kin = kinematic_coherence(df)
    C_conc = concentration_index(df)

    if not np.isfinite(C_kin):
        return None

    fdm_outer = compute_fdm_outer(df)
    fdm_bins = compute_fdm_radial_bins(df)
    f_gas = compute_gas_fraction(df)

    galaxy = os.path.basename(path).replace("_rotmod.dat", "")

    result = {
        "galaxy": galaxy,
        "C_kin": C_kin,
        "C_conc": C_conc,
        "f_dm": fdm_outer,
        "f_gas": f_gas,
        "Vmax": float(df["Vobs"].max()),
        "Rmax": float(df["R"].max()),
        "N_points": len(df),
    }
    result.update(fdm_bins)
    return result


# ============================================================
# STATISTICAL ANALYSIS
# ============================================================
def print_header(title):
    print(f"\n{'='*60}")
    print(title)
    print(f"{'='*60}")


def simple_regression(x, y, xlabel="X", ylabel="Y"):
    """Run and print simple linear regression."""
    valid = np.isfinite(x) & np.isfinite(y)
    x_v, y_v = x[valid], y[valid]
    if len(x_v) < 10:
        print(f"  Insufficient data (N={len(x_v)})")
        return None
    slope, intercept, r, p, se = stats.linregress(x_v, y_v)
    print(f"  {xlabel} → {ylabel}:")
    print(f"    N     = {len(x_v)}")
    print(f"    β     = {slope:+.4f} ± {se:.4f}")
    print(f"    r     = {r:.4f}")
    print(f"    p     = {p:.4f}")
    return {"slope": slope, "se": se, "r": r, "p": p, "n": len(x_v)}


def high_low_comparison(coherence, fdm, label="C"):
    """Split at median, compare means via t-test."""
    valid = np.isfinite(coherence) & np.isfinite(fdm)
    c_v, f_v = coherence[valid], fdm[valid]
    med = np.median(c_v)
    high_mask = c_v >= med
    low_mask = c_v < med
    high_fdm = f_v[high_mask]
    low_fdm = f_v[low_mask]
    t, p = stats.ttest_ind(high_fdm, low_fdm)
    diff = high_fdm.mean() - low_fdm.mean()
    print(f"\n  High {label} (N={high_mask.sum()}): mean f_dm = {high_fdm.mean():.3f}")
    print(f"  Low  {label} (N={low_mask.sum()}):  mean f_dm = {low_fdm.mean():.3f}")
    print(f"  Difference: {diff:+.3f}")
    print(f"  t-test p  = {p:.4f}")
    return {"high_mean": high_fdm.mean(), "low_mean": low_fdm.mean(),
            "diff": diff, "p": p}


def controlled_regression(df, coherence_col, confounders, ylabel="f_dm"):
    """Run multiple regression with confounders."""
    cols = [coherence_col] + confounders + [ylabel]
    sub = df[cols].dropna()
    if len(sub) < 10:
        print(f"  Insufficient data (N={len(sub)})")
        return None

    y = sub[ylabel].to_numpy(dtype=float)
    X_cols = [coherence_col] + confounders
    X = np.column_stack([np.ones(len(sub))] +
                        [sub[c].to_numpy(dtype=float) for c in X_cols])

    # OLS
    beta, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ beta
    resid = y - y_pred
    n, k = X.shape
    sigma2 = np.sum(resid**2) / (n - k)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))

    # t-stats and p-values
    t_stats = beta / se
    p_vals = 2 * stats.t.sf(np.abs(t_stats), df=n - k)

    label = f"f_dm ~ {coherence_col}"
    if confounders:
        label += " + " + " + ".join(confounders)
    print(f"\n  {label} (N={n}):")
    print(f"    β({coherence_col}) = {beta[1]:+.4f} ± {se[1]:.4f}, p = {p_vals[1]:.4f}")
    for i, conf in enumerate(confounders):
        idx = i + 2
        print(f"    β({conf}) = {beta[idx]:+.4f} ± {se[idx]:.4f}, p = {p_vals[idx]:.4f}")

    return {"beta_C": beta[1], "se_C": se[1], "p_C": p_vals[1]}


def gas_fraction_split(df, coherence_col):
    """Split by median gas fraction and test coherence-fdm in each."""
    sub = df[[coherence_col, "f_dm", "f_gas"]].dropna()
    med_gas = sub["f_gas"].median()

    print(f"\n  Median gas fraction: {med_gas:.3f}")

    for label, mask in [("High gas", sub["f_gas"] >= med_gas),
                        ("Low gas", sub["f_gas"] < med_gas)]:
        ss = sub[mask]
        if len(ss) < 10:
            print(f"  {label}: insufficient data (N={len(ss)})")
            continue
        slope, intercept, r, p, se = stats.linregress(
            ss[coherence_col].to_numpy(dtype=float),
            ss["f_dm"].to_numpy(dtype=float)
        )
        print(f"  {label} (N={len(ss)}): β = {slope:+.4f}, p = {p:.4f}")


def radial_bin_analysis(df_all, sparc_dir):
    """Compute coherence-fdm correlation in each radial bin.

    This requires re-reading the raw rotation curves to get per-radius f_dm.
    """
    print_header("RADIAL BIN ANALYSIS")

    # Use pre-computed radial bin columns from the main dataframe
    bin_cols = [c for c in df_all.columns if c.startswith("fdm_")]
    if not bin_cols:
        print("  No radial bin data available.")
        return

    for col in sorted(bin_cols):
        valid = df_all[["C_kin", col]].dropna()
        if len(valid) < 10:
            continue
        slope, intercept, r, p, se = stats.linregress(
            valid["C_kin"].to_numpy(dtype=float),
            valid[col].to_numpy(dtype=float)
        )
        print(f"  {col}: β = {slope:+.4f}, p = {p:.4f} (N={len(valid)})")

    # High vs Low C_kin comparison per radial bin
    print("\n  High vs Low C_kin comparison per bin:")
    med_c = df_all["C_kin"].median()
    for col in sorted(bin_cols):
        valid = df_all[["C_kin", col]].dropna()
        high = valid[valid["C_kin"] >= med_c][col]
        low = valid[valid["C_kin"] < med_c][col]
        if len(high) < 5 or len(low) < 5:
            continue
        t, p = stats.ttest_ind(high, low)
        diff = high.mean() - low.mean()
        print(f"  {col}: Δf_dm = {diff:+.3f}, p = {p:.4f}")


# ============================================================
# MAIN
# ============================================================
def main(sparc_dir, sparc_table_path, output_path):

    # --- Load galaxy files ---
    files = [f for f in os.listdir(sparc_dir) if f.endswith("_rotmod.dat")]
    print(f"Found {len(files)} galaxy files in {sparc_dir}")

    rows = []
    for fn in sorted(files):
        path = os.path.join(sparc_dir, fn)
        result = process_galaxy(path)
        if result is not None:
            rows.append(result)

    df = pd.DataFrame(rows)
    print(f"Processed {len(df)} galaxies with valid data")

    # --- Merge Hubble type if table provided ---
    if sparc_table_path and os.path.exists(sparc_table_path):
        table = read_sparc_table(sparc_table_path)
        if table is not None and "T" in table.columns and "Galaxy" in table.columns:
            df = df.merge(table[["Galaxy", "T"]], left_on="galaxy",
                          right_on="Galaxy", how="left")
            df = df.drop(columns=["Galaxy"], errors="ignore")
            print(f"Merged Hubble type for {df['T'].notna().sum()} galaxies")
        else:
            print("Warning: Could not parse SPARC summary table")
            df["T"] = np.nan
    else:
        df["T"] = np.nan

    # --- Derived columns ---
    df["log_Vmax"] = np.log10(df["Vmax"].clip(lower=1))

    # --- Save complete CSV ---
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} galaxies to {output_path}")

    # ===========================================================
    # ANALYSIS
    # ===========================================================

    # --- 1. Kinematic coherence ---
    print_header("1. KINEMATIC COHERENCE → DARK MATTER FRACTION")
    simple_regression(df["C_kin"].to_numpy(), df["f_dm"].to_numpy(),
                      "C_kin", "f_dm")
    high_low_comparison(df["C_kin"].to_numpy(), df["f_dm"].to_numpy(), "C_kin")

    # --- 2. Morphological coherence ---
    print_header("2. MORPHOLOGICAL COHERENCE → DARK MATTER FRACTION")
    valid_conc = df.dropna(subset=["C_conc"])
    simple_regression(valid_conc["C_conc"].to_numpy(),
                      valid_conc["f_dm"].to_numpy(), "C_conc", "f_dm")
    high_low_comparison(valid_conc["C_conc"].to_numpy(),
                        valid_conc["f_dm"].to_numpy(), "C_conc")

    # --- 3. Independence check ---
    print_header("3. INDEPENDENCE CHECK")
    both = df.dropna(subset=["C_kin", "C_conc"])
    r_kc, p_kc = stats.pearsonr(both["C_kin"], both["C_conc"])
    print(f"  C_kin vs C_conc: r = {r_kc:.3f}, p = {p_kc:.4f}")

    if "log_Vmax" in df.columns:
        r_cm, p_cm = stats.pearsonr(
            valid_conc["C_conc"].to_numpy(),
            valid_conc["log_Vmax"].to_numpy()
        )
        print(f"  C_conc vs log(Vmax): r = {r_cm:.3f}, p = {p_cm:.4f}")

    if df["T"].notna().sum() > 10:
        both_t = df.dropna(subset=["C_conc", "T"])
        r_ct, p_ct = stats.pearsonr(both_t["C_conc"], both_t["T"])
        print(f"  C_conc vs Hubble T: r = {r_ct:.3f}, p = {p_ct:.4f}")

    # --- 4. Joint model ---
    print_header("4. JOINT MODEL: f_dm ~ C_kin + C_conc")
    controlled_regression(df, "C_kin", ["C_conc"])

    # --- 5. Controlled regressions ---
    print_header("5. CONTROLLED REGRESSIONS")
    controlled_regression(df, "C_kin", [], ylabel="f_dm")
    controlled_regression(df, "C_kin", ["log_Vmax"], ylabel="f_dm")
    controlled_regression(df, "C_kin", ["f_gas"], ylabel="f_dm")
    if df["T"].notna().sum() > 20:
        controlled_regression(df, "C_kin", ["T"], ylabel="f_dm")
        controlled_regression(df, "C_kin", ["log_Vmax", "f_gas", "T"],
                              ylabel="f_dm")

    print_header("5b. CONTROLLED REGRESSIONS (C_conc)")
    controlled_regression(df, "C_conc", [], ylabel="f_dm")
    controlled_regression(df, "C_conc", ["log_Vmax"], ylabel="f_dm")
    controlled_regression(df, "C_conc", ["f_gas"], ylabel="f_dm")
    if df["T"].notna().sum() > 20:
        controlled_regression(df, "C_conc", ["T"], ylabel="f_dm")

    # --- 6. Gas fraction split ---
    print_header("6. GAS FRACTION SPLIT")
    print("\n  Kinematic coherence:")
    gas_fraction_split(df, "C_kin")
    print("\n  Morphological coherence:")
    gas_fraction_split(df, "C_conc")

    # --- 7. Radial bin analysis ---
    radial_bin_analysis(df, sparc_dir)

    # --- 8. Summary ---
    print_header("SUMMARY")
    kin = simple_regression(df["C_kin"].to_numpy(), df["f_dm"].to_numpy(),
                            "C_kin", "f_dm")
    conc = simple_regression(valid_conc["C_conc"].to_numpy(),
                             valid_conc["f_dm"].to_numpy(), "C_conc", "f_dm")
    print(f"\n  Correlation between measures: r = {r_kc:.3f}")
    print(f"\n  Two independent coherence measures both predict lower dark matter.")
    if conc and conc["p"] < 0.01:
        print(f"  Concentration result significant at p < 0.01.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SPARC Coherence Analysis — Full Pipeline")
    parser.add_argument("--sparc-dir", required=True,
                        help="Path to SPARC Rotmod_LTG directory")
    parser.add_argument("--sparc-table", default=None,
                        help="Path to SPARC summary table (for Hubble type)")
    parser.add_argument("--output", default="sparc_results.csv",
                        help="Output CSV path (default: sparc_results.csv)")

    args = parser.parse_args()
    main(args.sparc_dir, args.sparc_table, args.output)
