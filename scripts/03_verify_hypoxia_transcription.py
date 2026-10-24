"""Step 3: double-check the hand-transcribed hypoxic area series (data/hypoxia_area_lumcon.csv).

Checks
  a) Text check: for each year whose value came from a press release, confirm that the
     transcribed number (with thousands comma) appears in the downloaded release text.
  b) Chart check: digitize the bar heights of the LUMCON 1985-2023 histogram
     (bar pixels between the 0 and 25,000 km2 gridlines) and compare with the series and
     with the LUMCON 2026 web table.
  c) Long-term means stated in LSU/LUMCON documents vs means of the series.

Output: outputs/tables/hypoxia_transcription_check.csv and a printed summary.
Run:    python -I scripts/03_verify_hypoxia_transcription.py
"""
from pathlib import Path
import os
import re
import html
import logging

import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

logging.getLogger("pypdf").setLevel(logging.ERROR)

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("GHM_RAW_DIR", ROOT / "data_raw"))
LUM = RAW / "lumcon"


def page_text(path):
    if path.suffix == ".pdf":
        t = " ".join((p.extract_text() or "") for p in PdfReader(str(path)).pages)
    else:
        t = path.read_text(errors="ignore")
        t = re.sub(r"<script.*?</script>|<style.*?</style>", " ", t, flags=re.S)
        t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return re.sub(r"\s+", " ", t)


def digitize_histogram(png):
    im = np.asarray(Image.open(png).convert("RGB")).astype(int)
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    # light-gray horizontal gridlines (0, 5000, ..., 25000 km2)
    gray = (abs(r - g) < 8) & (abs(g - b) < 8) & (r > 190) & (r < 240)
    rows = [y for y in range(im.shape[0]) if gray[y, 180:1200].sum() > 600]
    y_top, y0 = min(rows), max(rows)
    km2_per_px = 25000 / (y0 - y_top)
    bar = (abs(r - 91) < 25) & (abs(g - 155) < 25) & (abs(b - 213) < 25)
    # bar columns: contiguous runs of bar pixels on the row just above the baseline
    cols = np.where(bar[y0 - 3, :])[0]
    runs, s, p = [], cols[0], cols[0]
    for c in cols[1:]:
        if c != p + 1:
            runs.append((s, p))
            s = c
        p = c
    runs.append((s, p))
    runs = [(a, c) for a, c in runs if c - a >= 4]
    out = []
    for a, c in runs:
        x = (a + c) // 2
        y, gap = y0 - 3, 0
        top = y
        while y > y_top - 5 and gap <= 6:  # tolerate thin red/annotation lines across a bar
            if bar[y, x]:
                top, gap = y, 0
            else:
                gap += 1
            y -= 1
        out.append(((a + c) / 2, (y0 - top) * km2_per_px))
    xs = np.array([o[0] for o in out])
    # 1985-2023 = 39 slots, 1988 (40 km2), 1989 and 2016 have no visible bar.
    first, last = xs.min(), xs.max()  # 1985 and 2023 bars
    step = (last - first) / (2023 - 1985)
    years = np.rint(1985 + (xs - first) / step).astype(int)
    return pd.DataFrame({"year": years, "histogram_km2": [o[1] for o in out]}), km2_per_px


def main():
    d = pd.read_csv(ROOT / "data" / "hypoxia_area_lumcon.csv")

    # a) text check
    found = []
    for _, row in d.iterrows():
        y = int(row.year)
        src = str(row.value_source).lower()
        if pd.isna(row.area_km2) or not ("release" in src or "report" in src):
            found.append(np.nan)
            continue
        f = next((LUM / f"pr_{y}{ext}" for ext in (".pdf", ".html") if (LUM / f"pr_{y}{ext}").exists()), None)
        if y == 1999:
            f = LUM / "pr_2000.html"
        if f is None:
            found.append(np.nan)
            continue
        t = page_text(f).replace(" ,", ",")
        v = int(row.area_km2)
        pats = [f"{v:,}", str(v)]
        if y == 2020:  # value derived from the square-mile figure
            pats = ["2,117 square miles"]
        found.append(int(any(p in t for p in pats)))
    d["text_check_found"] = found

    # b) chart check
    h, kpp = digitize_histogram(LUM / "lumcon_histogram_1985_2023.png")
    d = d.merge(h, on="year", how="left")
    text_years = d[(d.text_check_found == 1) & d.histogram_km2.notna()]
    bias = float(np.median(text_years.histogram_km2 - text_years.area_km2))
    d["histogram_bias_corrected_km2"] = (d.histogram_km2 - bias).round(-1)
    d["diff_series_minus_hist"] = d.area_km2 - d.histogram_bias_corrected_km2
    d["diff_table_minus_hist"] = d.lumcon_table_km2 - d.histogram_bias_corrected_km2
    cols = ["year", "area_km2", "lumcon_table_km2", "histogram_bias_corrected_km2",
            "diff_series_minus_hist", "diff_table_minus_hist", "text_check_found", "value_source"]
    out = d[cols]
    out.to_csv(ROOT / "outputs" / "tables" / "hypoxia_transcription_check.csv", index=False)

    print(f"histogram scale {kpp:.1f} km2 per pixel; median digitizing bias {bias:+.0f} km2 "
          f"(from {len(text_years)} press-release years)")
    n_txt = int(np.nansum(d.text_check_found))
    n_try = int(d.text_check_found.notna().sum())
    print(f"text check: {n_txt}/{n_try} press-release values found verbatim in the downloaded text")
    big = out[(out.diff_series_minus_hist.abs() > 300)]
    print("series years differing from the histogram by > 300 km2 (red lines cross 1987 and 2021 bars):")
    print(big.to_string(index=False))
    tab = out[(out.diff_table_minus_hist.abs() > 300)]
    print("LUMCON 2026 web-table years differing from the histogram by > 300 km2:")
    print(tab[["year", "lumcon_table_km2", "histogram_bias_corrected_km2", "area_km2"]].to_string(index=False))

    # c) stated long-term means
    s = d.set_index("year").area_km2
    checks = [("1985-2011 mean, LSU 2012 forecast", 1985, 2011, 14046),
              ("1985-2020 mean (n=34), LSU 2021 forecast", 1985, 2020, 13789),
              ("38-year mean, LSU 2024 release", 1985, 2024, 13675)]
    print("stated long-term means vs this series (and vs the 2026 web table):")
    t = d.set_index("year").lumcon_table_km2
    for name, a, b, stated in checks:
        print(f"  {name}: stated {stated}, series {s.loc[a:b].mean():.0f} (n={s.loc[a:b].notna().sum()}), "
              f"web table {t.loc[a:b].mean():.0f}")


if __name__ == "__main__":
    main()
