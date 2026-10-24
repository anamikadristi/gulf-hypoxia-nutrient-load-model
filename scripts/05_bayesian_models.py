"""Step 5: Bayesian regression of hypoxic area on May nitrate load (PyMC), and the
May nitrate load needed to bring the expected mid-summer area down to the
Hypoxia Task Force goal of 5,000 km2, with posterior uncertainty.

Models (Student-t likelihood, so single odd years such as 1988 do not dominate):
  B1 linear : area = a + b * (N - 130)                       undisturbed years
  B2 log    : area = a + b * ln(N / 130)                     undisturbed years
  B3 storm  : area = (a + b * (N - 130)) * (1 - s * D)       all years, D = storm/wind flag
N = May nitrate+nitrite load to the Gulf (thousand metric tons, kt).

Outputs: outputs/tables/bayes_parameters.csv, bayes_loo.csv, load_reduction.csv,
         prob_below_goal_curve.csv, outputs/tables/bayes_draws_B1_B2.csv (thinned draws)
Run:     python -I scripts/05_bayesian_models.py
"""
from pathlib import Path
import os

# Pure-Python PyTensor backend: slower but avoids C compiler setup problems.
os.environ.setdefault("PYTENSOR_FLAGS", "cxx=")

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "outputs" / "tables"
GOAL = 5000.0
NC = 130.0  # centering value for May NO3 load, kt
SEED = 20261023


def sample(model):
    with model:
        idata = pm.sample(2000, tune=2000, chains=4, target_accept=0.95, random_seed=SEED,
                          progressbar=False)
        pm.compute_log_likelihood(idata, progressbar=False)
    return idata


def priors():
    a = pm.Normal("a", 14000, 6000)          # area at N = 130 kt
    b = pm.Normal("b", 0, 300)               # km2 per kt (B1) or per ln unit scaled below
    sigma = pm.HalfNormal("sigma", 6000)
    nu = pm.Gamma("nu", 2, 0.1)
    return a, b, sigma, nu


def main():
    d = pd.read_csv(ROOT / "data" / "model_dataset.csv").dropna(subset=["area_km2", "no3_kt"])
    calm = d[d.disturbed == 0]
    N, A = calm.no3_kt.values, calm.area_km2.values

    with pm.Model() as b1:
        a, b, sigma, nu = priors()
        pm.StudentT("area", nu=nu, mu=a + b * (N - NC), sigma=sigma, observed=A)
    with pm.Model() as b2:
        a = pm.Normal("a", 14000, 6000)
        b = pm.Normal("b", 0, 30000)         # km2 per unit ln(load)
        sigma = pm.HalfNormal("sigma", 6000)
        nu = pm.Gamma("nu", 2, 0.1)
        pm.StudentT("area", nu=nu, mu=a + b * np.log(N / NC), sigma=sigma, observed=A)
    Na, Aa, D = d.no3_kt.values, d.area_km2.values, d.disturbed.values
    with pm.Model() as b3:
        a, b, sigma, nu = priors()
        s = pm.Beta("s", 2, 2)               # fractional reduction of area in disturbed years
        pm.StudentT("area", nu=nu, mu=(a + b * (Na - NC)) * (1 - s * D), sigma=sigma, observed=Aa)

    res = {}
    for name, m in [("B1 linear (undisturbed)", b1), ("B2 log load (undisturbed)", b2),
                    ("B3 linear + storm factor (all years)", b3)]:
        print("sampling", name)
        res[name] = sample(m)

    # parameter summaries and diagnostics
    rows = []
    for name, idata in res.items():
        sm = az.summary(idata, ci_prob=0.90)
        sm.insert(0, "model", name)
        rows.append(sm.reset_index().rename(columns={"index": "param"}))
    pars = pd.concat(rows, ignore_index=True)
    pars.to_csv(TAB / "bayes_parameters.csv", index=False)
    print(pars[["model", "param", "mean", "sd", "eti90_lb", "eti90_ub", "r_hat", "ess_bulk"]].to_string(index=False))

    # PSIS-LOO (approximate leave-one-year-out) for the two undisturbed-year models
    loo_rows = []
    for name in ["B1 linear (undisturbed)", "B2 log load (undisturbed)"]:
        lo = az.loo(res[name], pointwise=True)
        k = np.asarray(lo.pareto_k).ravel()
        loo_rows.append({"model": name, "elpd_loo": float(lo.elpd), "se": float(lo.se), "p_loo": float(lo.p),
                         "max_pareto_k": float(np.max(k)), "good_k_threshold": float(lo.good_k),
                         "n_k_over_threshold": int((k > lo.good_k).sum())})
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(TAB / "bayes_loo.csv", index=False)
    print(loo.round(2).to_string(index=False))

    # Load needed for expected (storm-free) area = 5,000 km2
    base_9601 = d[(d.year >= 1985) & (d.year <= 1996)].no3_kt.mean()
    allq = pd.read_csv(ROOT / "data" / "usgs_may_loads_wrtdsk.csv")
    base_8096 = allq[(allq.year >= 1980) & (allq.year <= 1996)].may_no3_t.mean() / 1000
    recent = allq[(allq.year >= 2016) & (allq.year <= 2025)].may_no3_t.mean() / 1000
    out, draws = [], {}
    for name, kind in [("B1 linear (undisturbed)", "lin"), ("B2 log load (undisturbed)", "log")]:
        post = res[name].posterior
        a = post["a"].values.ravel()
        b = post["b"].values.ravel()
        if kind == "lin":
            nstar = NC + (GOAL - a) / b
        else:
            nstar = NC * np.exp((GOAL - a) / b)
        nstar = np.where(b > 0, nstar, np.nan)
        nstar = np.clip(nstar, 0, None)  # a load below zero is not physical
        draws[name] = (a, b, post["sigma"].values.ravel(), post["nu"].values.ravel(), nstar)
        q = np.nanpercentile(nstar, [5, 50, 95])
        for ref_name, ref in [("1980-1996 mean May load", base_8096), ("2016-2025 mean May load", recent)]:
            red = 100 * (1 - nstar / ref)
            rq = np.nanpercentile(red, [5, 50, 95])
            out.append({"model": name, "reference": ref_name, "reference_kt": ref,
                        "required_load_kt_median": q[1], "required_load_kt_5": q[0], "required_load_kt_95": q[2],
                        "reduction_pct_median": rq[1], "reduction_pct_5": rq[0], "reduction_pct_95": rq[2],
                        "share_draws_b_le_0": float(np.mean(b <= 0)),
                        "share_draws_needing_load_at_or_below_0": float(np.nanmean(nstar <= 0))})
    lr = pd.DataFrame(out)
    lr.to_csv(TAB / "load_reduction.csv", index=False)
    print(f"\nreference May NO3 loads: 1980-1996 {base_8096:.1f} kt, 2016-2025 {recent:.1f} kt "
          f"(1985-1996 {base_9601:.1f} kt)")
    print(lr.round(1).to_string(index=False))

    # Probability that a storm-free year's area is below 5,000 km2, as a function of load
    rng = np.random.default_rng(SEED)
    grid = np.arange(5, 261, 5.0)
    curve = {"no3_kt": grid}
    for name, (a, b, sig, nu, _) in draws.items():
        mu = a[:, None] + (b[:, None] * (grid - NC) if "linear" in name else b[:, None] * np.log(grid / NC))
        eps = sig[:, None] * rng.standard_t(nu[:, None], size=mu.shape)
        ypp = mu + eps
        tag = "B1" if name.startswith("B1") else "B2"
        curve[f"{tag}_p_below_goal"] = np.mean(ypp < GOAL, axis=0)
        curve[f"{tag}_mean_med"] = np.median(mu, axis=0)
        curve[f"{tag}_mean_5"] = np.percentile(mu, 5, axis=0)
        curve[f"{tag}_mean_95"] = np.percentile(mu, 95, axis=0)
        curve[f"{tag}_pred_5"] = np.percentile(ypp, 5, axis=0)
        curve[f"{tag}_pred_95"] = np.percentile(ypp, 95, axis=0)
    cv = pd.DataFrame(curve)
    cv.round(4).to_csv(TAB / "prob_below_goal_curve.csv", index=False)
    for tag in ["B1", "B2"]:
        for p in [0.5, 0.8]:
            ok = cv[cv[f"{tag}_p_below_goal"] >= p]
            val = ok.no3_kt.max() if len(ok) else np.nan
            print(f"{tag}: highest May load on grid with P(storm-free area < 5,000) >= {p}: {val} kt")

    # storm factor from B3
    s = res["B3 linear + storm factor (all years)"].posterior["s"].values.ravel()
    print(f"B3 storm/wind factor s: median {np.median(s):.2f}, 90% interval "
          f"{np.percentile(s, 5):.2f}-{np.percentile(s, 95):.2f} (fractional area reduction)")

    # thinned draws for figures
    a1, b1_, s1, n1, ns1 = draws["B1 linear (undisturbed)"]
    a2, b2_, s2, n2, ns2 = draws["B2 log load (undisturbed)"]
    k = slice(None, None, 8)
    pd.DataFrame({"B1_a": a1[k], "B1_b": b1_[k], "B1_nstar": ns1[k],
                  "B2_a": a2[k], "B2_b": b2_[k], "B2_nstar": ns2[k],
                  "B3_s": s[k]}).round(3).to_csv(TAB / "bayes_draws_B1_B2.csv", index=False)


if __name__ == "__main__":
    main()
