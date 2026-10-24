# Step 2: extract monthly (May) nitrate+nitrite, total nitrogen and total phosphorus
# loads and May mean discharge from the USGS WRTDS-Kalman daily model output.
#
# Input : data_raw/usgs_gulf_loads/eLists.zip (USGS data release, Klager 2026, ver. 3.0)
# Output: data/usgs_monthly_loads_wrtdsk.csv  (all months, each site and the Gulf sum)
#         data/usgs_may_loads_wrtdsk.csv      (May only, Gulf sum = STFR + MELV)
#
# Base R only (readRDS); no packages needed.  Run:  Rscript scripts/02_extract_may_loads.R
#
# EGRET conventions in the Daily table: Q in m3/s; GenFlux (WRTDS-Kalman daily flux) in kg/day.
# STFR = Mississippi River at St. Francisville, LA (QW 07373420; flow from USACE Tarbert Landing)
# MELV = Atchafalaya River at Melville, LA (QW 07381495; flow from USACE Simmesport)

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(f), ".."))
raw <- Sys.getenv("GHM_RAW_DIR", file.path(root, "data_raw"))
zipf <- file.path(raw, "usgs_gulf_loads", "eLists.zip")
tmp <- file.path(tempdir(), "elists")
dir.create(tmp, showWarnings = FALSE)

constits <- c("NO3_NO2", "TN", "TP")
sites <- c("STFR", "MELV")
want <- as.vector(outer(sites, constits, function(s, c) paste0("WRTDS_K_", s, "_", c, ".rds")))
unzip(zipf, files = want, exdir = tmp, overwrite = TRUE)

rows <- list()
for (s in sites) for (cn in constits) {
  e <- readRDS(file.path(tmp, paste0("WRTDS_K_", s, "_", cn, ".rds")))
  d <- e$Daily
  d$year <- as.integer(format(d$Date, "%Y"))
  d$month <- as.integer(format(d$Date, "%m"))
  a <- aggregate(cbind(GenFlux, Q) ~ year + month, data = d,
                 FUN = function(x) c(sum = sum(x), mean = mean(x), n = length(x)))
  out <- data.frame(site = s, constit = cn, year = a$year, month = a$month,
                    load_t = a$GenFlux[, "sum"] / 1000,        # kg -> metric tons
                    q_mean_m3s = a$Q[, "mean"], ndays = a$GenFlux[, "n"])
  rows[[paste(s, cn)]] <- out
}
m <- do.call(rbind, rows)
m <- m[m$year >= 1980, ]

# Gulf total = STFR + MELV (the USGS convention for delivery to the Gulf)
key <- paste(m$constit, m$year, m$month)
g <- aggregate(cbind(load_t, q_mean_m3s) ~ constit + year + month, data = m, FUN = sum)
nsite <- aggregate(site ~ constit + year + month, data = m, FUN = length)
g <- merge(g, nsite)
g <- g[g$site == 2, ]
g$site <- "GULF"
g$ndays <- NA
m <- rbind(m, g[, names(m)])
m <- m[order(m$site, m$constit, m$year, m$month), ]
write.csv(m, file.path(root, "data", "usgs_monthly_loads_wrtdsk.csv"), row.names = FALSE)

may <- m[m$site == "GULF" & m$month == 5, ]
w <- reshape(may[, c("constit", "year", "load_t", "q_mean_m3s")], idvar = "year",
             timevar = "constit", direction = "wide")
out <- data.frame(year = w$year,
                  may_no3_t = round(w$load_t.NO3_NO2),
                  may_tn_t = round(w$load_t.TN),
                  may_tp_t = round(w$load_t.TP),
                  may_q_m3s = round(w$q_mean_m3s.NO3_NO2))
out <- out[order(out$year), ]
write.csv(out, file.path(root, "data", "usgs_may_loads_wrtdsk.csv"), row.names = FALSE)

# QA: annual (water-year) sums from daily output vs the published annual GULF WRTDS_K loads
ann <- m[m$site == "GULF", ]
ann$wy <- ifelse(ann$month >= 10, ann$year + 1, ann$year)
ann <- aggregate(load_t ~ constit + wy, data = ann, FUN = sum)
pub <- read.csv(file.path(raw, "usgs_gulf_loads", "WRTDS_loads.csv"), comment.char = "#",
                colClasses = c(SITE_QW_ID = "character", SITE_FLOW_ID = "character"))
pub <- pub[pub$SITE_ABB == "GULF" & pub$MODTYPE == "WRTDS_K", c("CONSTIT", "WY", "TONS")]
qa <- merge(ann, pub, by.x = c("constit", "wy"), by.y = c("CONSTIT", "WY"))
qa <- qa[qa$wy >= 1981 & qa$wy <= 2025, ]
qa$ratio <- qa$load_t / qa$TONS
write.csv(qa, file.path(root, "outputs", "tables", "qa_daily_sum_vs_published_annual.csv"), row.names = FALSE)
cat("May rows:", nrow(out), " years", min(out$year), "-", max(out$year), "\n")
cat("QA ratio (daily-sum metric t / published TONS), by constituent:\n")
print(aggregate(ratio ~ constit, qa, function(x) round(c(min = min(x), median = median(x), max = max(x)), 4)))
