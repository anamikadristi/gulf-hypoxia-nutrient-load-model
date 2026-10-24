"""Step 6: figures (PNG, 150 dpi) in outputs/figures/.
Run: python -I scripts/06_figures.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "outputs" / "tables"
FIG = ROOT / "outputs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
GOAL = 5000
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "legend.frameon": False, "figure.facecolor": "white", "axes.titleweight": "bold",
    "axes.titlesize": 11, "axes.titlelocation": "left",
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=150)
    plt.close(fig)


def fig_timeseries(h):
    fig, ax = plt.subplots(figsize=(9, 4))
    ok = h.dropna(subset=["area_km2"])
    calm, dist = ok[ok.disturbed == 0], ok[ok.disturbed == 1]
    ax.bar(calm.year, calm.area_km2, width=0.75, color=BLUE, label="Undisturbed cruise")
    ax.bar(dist.year, dist.area_km2, width=0.75, color=ORANGE, hatch="///", edgecolor="white",
           linewidth=0, label="Storm or unusual wind before/during cruise")
    run = h.set_index("year").area_km2.rolling(5, min_periods=4).mean()
    ax.plot(run.index, run.values, color=INK, lw=2, label="5-year running mean")
    ax.axhline(GOAL, color=INK2, ls="--", lw=1.2)
    ax.text(1984.6, GOAL + 300, "Task Force goal 5,000 km² (5-yr mean)", color=INK2, fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5), zorder=5)
    for y in [1989, 2016]:
        ax.text(y, 400, "no\ncruise", ha="center", fontsize=7, color=INK2)
    ax.set_ylabel("Mid-summer hypoxic area (km²)")
    ax.set_title("Northern Gulf hypoxic zone, 1985-2026 (LUMCON/LSU shelfwide cruises)")
    ax.set_xlim(1984, 2027)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 0.93), fontsize=8)
    save(fig, "fig1_hypoxia_timeseries.png")


def fig_loads(d):
    fig, axes = plt.subplots(3, 1, figsize=(9, 6.5), sharex=True)
    for ax, col, lab in zip(axes, ["no3_kt", "tp_kt", "q_kcms"],
                            ["May NO3+NO2 load (kt N)", "May TP load (kt P)", "May mean discharge (1000 m³/s)"]):
        ax.plot(d.year, d[col], color=BLUE, lw=2, marker="o", ms=4)
        ax.set_ylabel(lab, fontsize=9)
    axes[0].set_title("USGS WRTDS-Kalman May loads and discharge to the Gulf, 1985-2025")
    axes[-1].set_xlabel("Year")
    save(fig, "fig2_may_loads.png")


def fig_scatter(d, curve, draws):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    g = curve.no3_kt.values
    ax.fill_between(g, curve.B2_pred_5, curve.B2_pred_95, color=AQUA, alpha=0.12, lw=0,
                    label="B2 90% predictive interval")
    ax.fill_between(g, curve.B2_mean_5, curve.B2_mean_95, color=AQUA, alpha=0.30, lw=0)
    ax.plot(g, curve.B2_mean_med, color=AQUA, lw=2, label="B2 log-load model (median, 90% band)")
    ax.plot(g, curve.B1_mean_med, color=INK2, lw=1.5, ls=":", label="B1 linear model (median)")
    calm, dist = d[d.disturbed == 0], d[d.disturbed == 1]
    ax.scatter(calm.no3_kt, calm.area_km2, s=36, color=BLUE, edgecolor="white", lw=1.5, zorder=3,
               label="Undisturbed years")
    ax.scatter(dist.no3_kt, dist.area_km2, s=44, marker="D", color=ORANGE, edgecolor="white",
               lw=1.5, zorder=3, label="Storm/wind years (not fitted)")
    for _, r in d.iterrows():
        if r.year in (1988, 2000, 2012, 2017, 2020, 2018, 2002, 2023, 1998, 2013):
            ax.annotate(int(r.year), (r.no3_kt, r.area_km2), xytext=(4, 3), textcoords="offset points",
                        fontsize=7.5, color=INK2)
    ax.axhline(GOAL, color=INK2, ls="--", lw=1.2)
    ax.text(250, GOAL + 300, "5,000 km² goal", ha="right", fontsize=8, color=INK2)
    ax.set_xlim(0, 260)
    ax.set_ylim(-3000, 34000)
    ax.set_xlabel("May NO3+NO2 load to the Gulf (thousand metric tons N)")
    ax.set_ylabel("Mid-summer hypoxic area (km²)")
    ax.set_title("Hypoxic area vs May nitrate load, 1985-2025")
    ax.legend(fontsize=7.5, loc="upper left")
    save(fig, "fig3_area_vs_may_nitrate.png")


def fig_loyo(pred, cv):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), gridspec_kw={"width_ratios": [1, 1.2]})
    ax = axes[0]
    p = pred.dropna(subset=["loyo_M2"])
    ax.scatter(p.area_km2, p.loyo_M2, s=36, color=BLUE, edgecolor="white", lw=1.5, zorder=3)
    lim = [0, 25000]
    ax.plot(lim, lim, color=INK2, lw=1, ls="--")
    for _, r in p.iterrows():
        if abs(r.area_km2 - r.loyo_M2) > 6000:
            ax.annotate(int(r.year), (r.area_km2, r.loyo_M2), xytext=(4, 3), textcoords="offset points",
                        fontsize=7.5, color=INK2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Observed area (km²)")
    ax.set_ylabel("Leave-one-year-out forecast (km²)")
    ax.set_title("M2 (log load), undisturbed years")
    ax = axes[1]
    c = cv[(cv["sample"] == "undisturbed years") & ~cv.model.str.startswith("M3")].sort_values("skill_vs_climatology")
    ax.barh(c.model, c.skill_vs_climatology, color=BLUE, height=0.6)
    for y, v in enumerate(c.skill_vs_climatology):
        ax.text(v + (0.02 if v >= 0 else -0.02), y, f"{v:.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=8, color=INK)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_xlim(min(-0.2, c.skill_vs_climatology.min() - 0.3), 0.7)
    ax.set_xlabel("LOYO skill vs climatology (1 - MSE/MSE_clim)")
    ax.set_title("LOYO skill, 29 undisturbed years")
    ax.grid(axis="y", visible=False)
    save(fig, "fig4_loyo_cv.png")


def fig_reduction(draws, curve, lr):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    bins = np.arange(0, 131, 4)
    ax.hist(draws.B2_nstar.dropna(), bins=bins, color=AQUA, alpha=0.85, label="B2 log-load model")
    ax.hist(draws.B1_nstar.dropna(), bins=bins, histtype="step", color=INK2, lw=1.8, label="B1 linear model (bar at 0 = goal not\nreached even at zero load)")
    ref = lr[lr.reference.str.startswith("1980")].reference_kt.iloc[0]
    rec = lr[lr.reference.str.startswith("2016")].reference_kt.iloc[0]
    ax.axvline(ref, color=BLUE, lw=1.5)
    ax.text(ref + 1, ax.get_ylim()[1] * 0.92, f"1980-96 mean\n{ref:.0f} kt", fontsize=8, color=INK)
    ax.axvline(rec, color=BLUE, lw=1.5, ls="--")
    ax.text(rec - 1, ax.get_ylim()[1] * 0.50, f"2016-25 mean\n{rec:.0f} kt", fontsize=8, color=INK, ha="right")
    ax.set_xlabel("May NO3+NO2 load giving expected area = 5,000 km² (kt N)")
    ax.set_ylabel("Posterior draws")
    ax.set_title("Load needed to reach the goal")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0.08, 1.0))
    ax = axes[1]
    ax.plot(curve.no3_kt, curve.B2_p_below_goal, color=AQUA, lw=2, label="B2 log-load model")
    ax.plot(curve.no3_kt, curve.B1_p_below_goal, color=INK2, lw=1.5, ls=":", label="B1 linear model")
    ax.axvline(ref, color=BLUE, lw=1.5)
    ax.axvline(rec, color=BLUE, lw=1.5, ls="--")
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 260)
    ax.set_xlabel("May NO3+NO2 load (kt N)")
    ax.set_ylabel("P(storm-free year area < 5,000 km²)")
    ax.set_title("Chance a storm-free year is below goal")
    ax.legend(fontsize=8)
    save(fig, "fig5_load_reduction.png")


def fig_storms(st):
    fig, ax = plt.subplots(figsize=(8, 4))
    st = st.sort_values("year")
    x = np.arange(len(st))
    ax.vlines(x, st.pi90_lo, st.pi90_hi, color=GRID, lw=6, label="M1 90% prediction interval")
    ax.scatter(x, st.pred_km2, s=40, color=INK2, marker="_", lw=2, zorder=3, label="Predicted from May NO3 load")
    ax.scatter(x, st.area_km2, s=46, color=ORANGE, marker="D", edgecolor="white", lw=1.5, zorder=4,
               label="Observed")
    ax.set_xticks(x, st.year.astype(int))
    ax.set_ylabel("Hypoxic area (km²)")
    ax.set_title("Storm and wind years vs the load-only model fitted to undisturbed years")
    ax.legend(fontsize=8, loc="upper left", ncol=3)
    ax.set_ylim(0, 33000)
    save(fig, "fig6_storm_years.png")


def main():
    h = pd.read_csv(ROOT / "data" / "hypoxia_area_lumcon.csv")
    d = pd.read_csv(ROOT / "data" / "model_dataset.csv")
    loads = pd.read_csv(ROOT / "data" / "usgs_may_loads_wrtdsk.csv")
    loads = loads[loads.year >= 1985].assign(no3_kt=lambda x: x.may_no3_t / 1000,
                                             tp_kt=lambda x: x.may_tp_t / 1000,
                                             q_kcms=lambda x: x.may_q_m3s / 1000)
    curve = pd.read_csv(TAB / "prob_below_goal_curve.csv")
    draws = pd.read_csv(TAB / "bayes_draws_B1_B2.csv")
    lr = pd.read_csv(TAB / "load_reduction.csv")
    pred = pd.read_csv(TAB / "loyo_predictions.csv")
    cv = pd.read_csv(TAB / "loyo_cv_skill.csv")
    st = pd.read_csv(TAB / "storm_year_residuals.csv")
    fig_timeseries(h)
    fig_loads(loads)
    fig_scatter(d.dropna(subset=["area_km2", "no3_kt"]), curve, draws)
    fig_loyo(pred, cv)
    fig_reduction(draws, curve, lr)
    fig_storms(st)
    print("figures written:", sorted(p.name for p in FIG.glob("*.png")))


if __name__ == "__main__":
    main()
