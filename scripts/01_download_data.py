"""Step 1: download all public source files into data_raw/ (not tracked by git).

Sources (all public, no account or key needed):
  1. USGS data release, Klager (2026), "Trends in Nutrient and Other Water Quality
     Constituent Loads to the Gulf of America from the Mississippi and Atchafalaya
     Rivers, 1980-2025 (ver. 3.0, May 2026)", ScienceBase item 69fe4aa9b66b01610729c3fd.
     Includes eLists.zip (WRTDS and WRTDS-Kalman daily model output, about 94 MB).
  2. LUMCON / LSU shelfwide hypoxia pages and press releases (gulfhypoxia.net), used
     to transcribe and cross-check the annual hypoxic area series.
  3. USGS NWIS daily discharge, 07374000 (Mississippi R. at Baton Rouge) and
     07381490 (Atchafalaya R. at Simmesport), used only as a fallback cross-check.

Run:  python -I scripts/01_download_data.py
"""
from pathlib import Path
import json
import os
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
# Raw downloads go to data_raw/ (git-ignored). Set GHM_RAW_DIR to keep them elsewhere.
RAW = Path(os.environ.get("GHM_RAW_DIR", ROOT / "data_raw"))
UA = {"User-Agent": "Mozilla/5.0 (public-data research script)"}

SB_ITEM = "69fe4aa9b66b01610729c3fd"

# LUMCON / LSU pages and files used for the hypoxia area transcription.
LUMCON = {
    "shelfwide_cruises.html": "https://gulfhypoxia.net/research/shelfwide-cruises/",
    "lumcon_table_1985_2005.png": "https://gulfhypoxia.net/wp-content/uploads/2026/07/chart-1-2026-edit.png",
    "lumcon_table_2006_2026.png": "https://gulfhypoxia.net/wp-content/uploads/2026/07/chart-2-2026-edit.png",
    "lumcon_histogram_1985_2023.png": "https://gulfhypoxia.net/wp-content/uploads/2024/05/histogram-1985-2023-shelfwide-page-1.png",
    "pr_2000.html": "https://gulfhypoxia.net/mapping-of-the-dead-zone-completed/",
    "pr_2001.html": "https://gulfhypoxia.net/mapping-press-release-july-2001/",
    "pr_2002.html": "https://gulfhypoxia.net/summer-2002-hypoxia/",
    "pr_2003.html": "https://gulfhypoxia.net/tropical-storms-put-the-oxygen-back-into-the-dead-zone/",
    "pr_2004.html": "https://www.lsu.edu/seagrantfish/news/2004/deadzone.htm",
    "pr_2005.html": "https://gulfhypoxia.net/mapping-of-dead-zone-completed/",
    "pr_2006.html": "https://gulfhypoxia.net/2006-press-release/",
    "pr_2007.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/01/PressReleaseCruise2007.pdf",
    "pr_2008.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/01/PressRelease2008.pdf",
    "pr_2009.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/01/2009_Press_Release.pdf",
    "pr_2010.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/01/PressRelease2010.pdf",
    "pr_2011.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/02/PressRelease2011.pdf",
    "pr_2012.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/02/cruise2012-PressRelease.pdf",
    "pr_2013.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/02/cruise2013-PressRelease.pdf",
    "pr_2014.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/02/cruise2014-Hypoxia_Press_Release.pdf",
    "pr_2015.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/02/cruise2015-PRESS_RELEASE_2015.pdf",
    "pr_2017.pdf": "https://gulfhypoxia.net/wp-content/uploads/2017/08/rev-PRESS-RELEASE-2017-LUMCON-LSU.pdf",
    "pr_2018.pdf": "https://gulfhypoxia.net/wp-content/uploads/2018/07/LUMCONLSU-Press-Release-2018-1.pdf",
    "pr_2019.pdf": "https://gulfhypoxia.net/wp-content/uploads/2019/08/Press-release-LSU-LUMCON-2019.pdf",
    "pr_2020.pdf": "https://gulfhypoxia.net/wp-content/uploads/2020/08/2020-press-release.pdf",
    "pr_2021.pdf": "https://gulfhypoxia.net/wp-content/uploads/2021/08/2021-LSULUMCON-press-release.pdf",
    "pr_2022.pdf": "https://gulfhypoxia.net/wp-content/uploads/2022/08/edited-press-release-8-2-22.pdf",
    "pr_2023.pdf": "https://gulfhypoxia.net/wp-content/uploads/2023/08/LSU-LUMCON-hypoxia-press-release-8-2-23-0935.pdf",
    "pr_2024.pdf": "https://gulfhypoxia.net/wp-content/uploads/2024/08/Final-press-release-8-1-2024-NR.pdf",
    "pr_2025.pdf": "https://gulfhypoxia.net/wp-content/uploads/2025/08/LSU-Report-8-6-25.pdf",
    "pr_2026.pdf": "https://gulfhypoxia.net/wp-content/uploads/2026/08/DZ-size-in-2026_CG.pdf",
    "forecast_2022.pdf": "https://gulfhypoxia.net/wp-content/uploads/2022/06/2022-hypoxia-forecast_NR.pdf",
}

NWIS = "https://waterservices.usgs.gov/nwis/dv/?format=rdb&sites={site}&parameterCd=00060&startDT=1985-01-01&endDT=2026-09-30"


def fetch(url, dest, retries=3):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  exists  {dest.name}")
        return
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as f:
                f.write(r.read())
            print(f"  ok      {dest.name}  ({dest.stat().st_size/1e6:.2f} MB)")
            return
        except Exception as e:  # noqa: BLE001
            print(f"  retry {k+1}: {url}: {e}")
            time.sleep(3)
    print(f"  FAILED  {url}")


def main():
    (RAW / "usgs_gulf_loads").mkdir(parents=True, exist_ok=True)
    (RAW / "lumcon").mkdir(parents=True, exist_ok=True)
    (RAW / "nwis").mkdir(parents=True, exist_ok=True)

    print("1. USGS Gulf loads data release (ScienceBase)")
    meta_url = f"https://www.sciencebase.gov/catalog/item/{SB_ITEM}?format=json"
    meta_path = RAW / "usgs_gulf_loads" / "sciencebase_item.json"
    fetch(meta_url, meta_path)
    meta = json.loads(meta_path.read_text())
    print("   title:", meta.get("title"))
    for f in meta.get("files", []):
        fetch(f["url"], RAW / "usgs_gulf_loads" / f["name"])

    print("2. LUMCON / LSU hypoxia pages and press releases")
    for name, url in LUMCON.items():
        fetch(url, RAW / "lumcon" / name)

    print("3. USGS NWIS daily discharge (fallback cross-check)")
    for site in ["07374000", "07381490"]:
        fetch(NWIS.format(site=site), RAW / "nwis" / f"dv_{site}.rdb")

    (RAW / "ACCESS_DATE.txt").write_text(time.strftime("%Y-%m-%d") + "\n")


if __name__ == "__main__":
    main()
