"""Step 4: merge hypoxic area with USGS May loads, fit regression models, and score
leave-one-year-out (LOYO) cross-validated forecast skill.

Inputs : data/hypoxia_area_lumcon.csv, data/usgs_may_loads_wrtdsk.csv,
         data_raw/nwis/dv_07374000.rdb and dv_07381490.rdb (fallback discharge check)
Outputs: data/model_dataset.csv
         outputs/tables/ols_models.csv, loyo_cv_skill.csv, loyo_predictions.csv,
         storm_year_residuals.csv, nwis_discharge_check.csv
Run:     python -I scripts/04_build_dataset_and_regressions.py
"""
from pathlib import Path
import os

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GHM_RAW_DIR", ROOT / "data_raw"))
TAB = ROOT / "outputs" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

MODELS = {
    "M1 linear NO3":            ("area_km2 ~ no3_kt", None),
    "M2 log NO3":               ("area_km2 ~ np.log(no3_kt)", None),
    "M3 log area ~ log NO3":    ("np.log(area_km2) ~ np.log(no3_kt)", "log"),
    "M4 NO3 + TP":              ("area_km2 ~ no3_kt + tp_kt", None),
    "M5 NO3 + discharge":       ("area_km2 ~ no3_kt + q_kcms", None),
    "M6 NO3 + year":            ("area_km2 ~ no3_kt + I(year - 2000)", None),
    "M7 discharge only":        ("area_km2 ~ q_kcms", None),
    "M8 TP only":               ("area_km2 ~ tp_kt", None),
}


def build():
    h = pd.read_csv(ROOT / "data" / "hypoxia_area_lumcon.csv")
    q = pd.read_csv(ROOT / "data" / "usgs_may_loads_wrtdsk.csv")
    d = h.merge(q, on="year", how="left")
    d["no3_kt"] = d.may_no3_t / 1000
    d["tn_kt"] = d.may_tn_t / 1000
    d["tp_kt"] = d.may_tp_t / 1000
    d["q_kcms"] = d.may_q_m3s / 1000
    keep = ["year", "area_km2", "disturbed", "disturbance_note", "may_no3_t", "may_tn_t",
            "may_tp_t", "may_q_m3s", "no3_kt", "tn_kt", "tp_kt", "q_kcms"]
    d = d[keep]
    d.to_csv(ROOT / "data" / "model_dataset.csv", index=False)
    return d


def nwis_check():
    """May mean discharge from NWIS daily values (Baton Rouge + Simmesport) vs the
    USACE-based flow in the USGS data release. NWIS records start 2004 and 2009."""
    frames = []
    for site in ["07374000", "07381490"]:
        f = RAW / "nwis" / f"dv_{site}.rdb"
        if not f.exists():
            return None
        t = pd.read_csv(f, sep="\t", comment="#", dtype=str).iloc[1:]
        vcol = [c for c in t.columns if c.endswith("_00060_00003")][0]
        t = pd.DataFrame({"date": pd.to_datetime(t.datetime), "q_cfs": pd.to_numeric(t[vcol], errors="coerce")})
        t["site"] = site
        frames.append(t)
    t = pd.concat(frames)
    t = t[t.date.dt.month == 5]
    t["year"] = t.date.dt.year
    g = t.groupby(["year", "site"]).agg(q_cfs=("q_cfs", "mean"), n=("q_cfs", "count")).reset_index()
    g = g[g.n >= 28].pivot(index="year", columns="site", values="q_cfs").dropna()
    g["nwis_may_q_m3s"] = (g.sum(axis=1) * 0.0283168).round()
    rel = pd.read_csv(ROOT / "data" / "usgs_may_loads_wrtdsk.csv")[["year", "may_q_m3s"]]
    out = g[["nwis_may_q_m3s"]].reset_index().merge(rel, on="year", how="left")
    out["ratio_nwis_to_release"] = (out.nwis_may_q_m3s / out.may_q_m3s).round(3)
    out.to_csv(TAB / "nwis_discharge_check.csv", index=False)
    return out


def fit_all(d):
    rows = []
    for name, (f, kind) in MODELS.items():
        m = smf.ols(f, data=d).fit()
        rows.append({"model": name, "formula": f, "n": int(m.nobs), "r2": m.rsquared,
                     "adj_r2": m.rsquared_adj, "aic": m.aic,
                     "resid_sd": np.sqrt(m.scale),
                     "coefs": "; ".join(f"{k}={v:.4g} (p={m.pvalues[k]:.3g})" for k, v in m.params.items())})
    return pd.DataFrame(rows)


def loyo(d, formula, kind):
    preds = []
    for i in d.index:
        tr = d.drop(index=i)
        m = smf.ols(formula, data=tr).fit()
        p = float(m.predict(d.loc[[i]]).iloc[0])
        if kind == "log":  # back-transform with lognormal mean correction
            p = float(np.exp(p + m.scale / 2))
        preds.append(p)
    return np.array(preds)


def skill(obs, pred, train_means):
    err = pred - obs
    mse = np.mean(err ** 2)
    mse_clim = np.mean((train_means - obs) ** 2)  # LOYO climatology (mean of other years)
    return {"rmse": np.sqrt(mse), "mae": np.mean(np.abs(err)), "bias": np.mean(err),
            "skill_vs_climatology": 1 - mse / mse_clim, "r_obs_pred": np.corrcoef(obs, pred)[0, 1]}


def main():
    d = build()
    data = d.dropna(subset=["area_km2", "no3_kt"]).reset_index(drop=True)
    calm = data[data.disturbed == 0].reset_index(drop=True)
    calm_no88 = calm[calm.year != 1988].reset_index(drop=True)
    print(f"years with area and loads: {len(data)} ({data.year.min()}-{data.year.max()}); "
          f"undisturbed: {len(calm)}; disturbed: {int(data.disturbed.sum())}")

    nw = nwis_check()
    if nw is not None:
        print(f"NWIS check: {len(nw)} Mays, NWIS/release May discharge ratio median "
              f"{nw.ratio_nwis_to_release.median():.3f} (range {nw.ratio_nwis_to_release.min():.3f}-"
              f"{nw.ratio_nwis_to_release.max():.3f})")

    fits = []
    for label, sub in [("undisturbed years", calm), ("undisturbed, no 1988", calm_no88), ("all years", data)]:
        f = fit_all(sub)
        f.insert(0, "sample", label)
        fits.append(f)
    # all-years model with a storm/wind indicator
    m = smf.ols("area_km2 ~ no3_kt + disturbed", data=data).fit()
    fits.append(pd.DataFrame([{"sample": "all years", "model": "M9 NO3 + disturbed flag",
                               "formula": "area_km2 ~ no3_kt + disturbed", "n": int(m.nobs),
                               "r2": m.rsquared, "adj_r2": m.rsquared_adj, "aic": m.aic,
                               "resid_sd": np.sqrt(m.scale),
                               "coefs": "; ".join(f"{k}={v:.4g} (p={m.pvalues[k]:.3g})" for k, v in m.params.items())}]))
    fits = pd.concat(fits, ignore_index=True)
    fits.round(4).to_csv(TAB / "ols_models.csv", index=False)
    print(fits[["sample", "model", "n", "r2", "aic", "resid_sd"]].round(3).to_string(index=False))

    # LOYO cross-validation
    cv_rows, pred_tab = [], data[["year", "area_km2", "disturbed", "no3_kt"]].copy()
    for label, sub in [("undisturbed years", calm), ("undisturbed, no 1988", calm_no88)]:
        train_means = np.array([sub.area_km2.drop(index=i).mean() for i in sub.index])
        for name, (f, kind) in MODELS.items():
            p = loyo(sub, f, kind)
            s = skill(sub.area_km2.values, p, train_means)
            cv_rows.append({"sample": label, "model": name, "n": len(sub), **s})
            if label == "undisturbed years":
                pred_tab.loc[pred_tab.year.isin(sub.year), f"loyo_{name.split()[0]}"] = p
    # all years, M9 with the disturbed flag
    train_means = np.array([data.area_km2.drop(index=i).mean() for i in data.index])
    for name, f in [("M1 linear NO3", "area_km2 ~ no3_kt"),
                    ("M9 NO3 + disturbed flag", "area_km2 ~ no3_kt + disturbed")]:
        p = loyo(data, f, None)
        cv_rows.append({"sample": "all years", "model": name, "n": len(data),
                        **skill(data.area_km2.values, p, train_means)})
        if name.startswith("M9"):
            pred_tab["loyo_M9_all"] = p
    cv = pd.DataFrame(cv_rows)
    cv.round(4).to_csv(TAB / "loyo_cv_skill.csv", index=False)
    pred_tab.round(1).to_csv(TAB / "loyo_predictions.csv", index=False)
    print("\nLOYO cross-validated skill:")
    print(cv.round(3).to_string(index=False))

    # Sensitivity: use the LUMCON 2026 web-table values instead of the checked series
    h = pd.read_csv(ROOT / "data" / "hypoxia_area_lumcon.csv")[["year", "lumcon_table_km2"]]
    alt = calm.merge(h, on="year").dropna(subset=["lumcon_table_km2"]).assign(area_km2=lambda x: x.lumcon_table_km2)
    alt = alt.reset_index(drop=True)
    sens = []
    for lab, sub in [("checked series", calm), ("2026 web table", alt)]:
        m2 = smf.ols("area_km2 ~ np.log(no3_kt)", data=sub).fit()
        tm = np.array([sub.area_km2.drop(index=i).mean() for i in sub.index])
        sk = skill(sub.area_km2.values, loyo(sub, "area_km2 ~ np.log(no3_kt)", None), tm)
        sens.append({"area_series": lab, "n": int(m2.nobs), "slope_per_ln_kt": m2.params.iloc[1],
                     "r2": m2.rsquared, "loyo_rmse": sk["rmse"], "loyo_skill": sk["skill_vs_climatology"]})
    sens = pd.DataFrame(sens)
    sens.round(4).to_csv(TAB / "sensitivity_area_series.csv", index=False)
    print("\nSensitivity of M2 (undisturbed years) to the area transcription:")
    print(sens.round(3).to_string(index=False))

    # Storm-year anomalies: predict disturbed years from the model fit to undisturbed years
    m1 = smf.ols("area_km2 ~ no3_kt", data=calm).fit()
    st = data[data.disturbed == 1].copy()
    pr = m1.get_prediction(st).summary_frame(alpha=0.10)
    st["pred_km2"] = pr["mean"].values
    st["pi90_lo"] = pr["obs_ci_lower"].values
    st["pi90_hi"] = pr["obs_ci_upper"].values
    st["obs_minus_pred"] = st.area_km2 - st.pred_km2
    st["obs_over_pred"] = st.area_km2 / st.pred_km2
    st["below_90pi"] = (st.area_km2 < st.pi90_lo).astype(int)
    st = st[["year", "area_km2", "no3_kt", "pred_km2", "pi90_lo", "pi90_hi", "obs_minus_pred",
             "obs_over_pred", "below_90pi", "disturbance_note"]]
    st.round(2).to_csv(TAB / "storm_year_residuals.csv", index=False)
    print("\nDisturbed years vs M1 fitted on undisturbed years:")
    print(st.drop(columns="disturbance_note").round(2).to_string(index=False))
    print(f"median obs/pred in disturbed years: {st.obs_over_pred.median():.2f}; "
          f"{int(st.below_90pi.sum())}/{len(st)} below the 90% prediction interval")


if __name__ == "__main__":
    main()
