# Round-trip validation of the EXACT mode: simulate IPD -> render a bold
# survminer-style figure (thick line, "+" censor marks, labels without axis
# tick marks, thick axis lines) -> extract with km_extract_exact -> require
# that the certified outputs COVER the truth:
#   (1) the events count bounds contain the true death count;
#   (2) the median bounds contain the true median (to time-pixel precision);
#   (3) at probe times placed mid-plateau (plateaus wider than 0.8 months, so
#       pixel-level time shifts cannot move a death across a probe), the true
#       S(t) lies within the enumerated solutions' range.
# Exact (d, c) sequence membership is NOT asserted: which side of a riser a
# censor hugging it sits on is pixel-ambiguous by nature; the statistics are
# the contract.
# Run from the repository root: Rscript example/exact_roundtrip.R
source("R/km_extract.R")
source("R/km_extract_exact.R")
suppressPackageStartupMessages(library(survival))

render_bold_km <- function(true, grid, f, lwd_curve = 9, lwd_axis = 7) {
  fit <- survfit(Surv(time, event) ~ 1, data = true)
  png(f, width = 1500, height = 1100)
  par(mar = c(6, 7, 3, 2), lwd = 2)
  plot(fit, conf.int = FALSE, mark.time = TRUE, lwd = lwd_curve,
       cex = 2.4, xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "", xlim = c(0, max(grid)), ylim = c(0, 1))
  usr <- par("usr")
  segments(usr[1], usr[3], usr[2], usr[3], lwd = lwd_axis, xpd = TRUE)
  segments(usr[1], usr[3], usr[1], usr[4], lwd = lwd_axis, xpd = TRUE)
  mtext(grid, side = 1, at = grid, line = 1.2, cex = 2.2)
  yl <- c(1, 0.75, 0.5, 0.25)
  mtext(sprintf("%.2f", yl), side = 2, at = yl, line = 0.8, cex = 2.2, las = 1)
  dev.off()
}

anchors_from <- function(true, grid) {
  data.frame(t = grid,
             n = vapply(grid, function(tt) sum(true$time >= tt), 0L),
             E = vapply(grid, function(tt) sum(true$event == 1 & true$time < tt), 0L))
}

set.seed(2026)
cases <- list(
  list(name = "N=60 mixed censoring",   N = 60, cens = 0.4),
  list(name = "N=20 sparse tail",       N = 20, cens = 0.35),
  list(name = "N=60 heavy ties",        N = 60, cens = 0.5, tie = TRUE),
  list(name = "N=40 censors at deaths", N = 40, cens = 0.45, adjacent = TRUE)
)

fails <- 0L
for (cs in cases) {
  N <- cs$N
  tt <- round(rexp(N, log(2) / 9), 2)
  ev <- rbinom(N, 1, 1 - cs$cens)
  if (isTRUE(cs$tie)) { i <- sample(N, 6); tt[i] <- rep(tt[i[1:2]], 3) }
  if (isTRUE(cs$adjacent)) {
    dd <- which(ev == 1)[1:4]; cc <- which(ev == 0)[1:4]
    tt[cc] <- tt[dd] + c(-0.04, 0.04, -0.04, 0.04)
  }
  tt <- pmin(pmax(tt, 0.3), 30); ev[tt >= 30] <- 0
  true <- data.frame(time = tt, event = ev)
  grid <- c(0, 10, 20, 30)
  f <- tempfile(fileext = ".png")
  render_bold_km(true, grid, f)
  res <- try(suppressMessages(km_extract_exact(
    f, x_ticks = grid, anchors = anchors_from(true, grid),
    max_solutions = 3e5)), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat(sprintf("%-24s FAIL (error: %s)\n", cs$name,
                trimws(conditionMessage(attr(res, "condition")))))
    fails <- fails + 1L; next
  }
  sf <- survfit(Surv(time, event) ~ 1, data = true)
  med_true <- unname(summary(sf)$table[["median"]])
  med_ok <- is.na(med_true) || all(is.na(res$bounds$median)) ||
    (med_true >= res$bounds$median[1] - 0.35 &&
     med_true <= res$bounds$median[2] + 0.35)
  ev_ok <- sum(true$event) >= res$bounds$events[1] &&
           sum(true$event) <= res$bounds$events[2]
  dts <- sort(unique(true$time[true$event == 1]))
  wid <- diff(c(dts, max(true$time)))
  probes <- (dts + wid / 2)[wid > 0.8]
  probes <- head(probes[order(-wid[wid > 0.8])], 4)
  s_ok <- TRUE
  for (tp in probes) {
    S_true <- summary(sf, times = tp, extend = TRUE)$surv
    S_all <- vapply(res$solutions, function(ipd)
      summary(survfit(Surv(time, event) ~ 1, data = ipd),
              times = tp, extend = TRUE)$surv, 0)
    if (S_true < min(S_all) - 1e-9 || S_true > max(S_all) + 1e-9) s_ok <- FALSE
  }
  status <- if (med_ok && ev_ok && s_ok) "PASS" else "FAIL"
  if (status == "FAIL") fails <- fails + 1L
  cat(sprintf("%-24s %s  (sols %d, med %s in [%.2f, %.2f], ev %d in [%d, %d], S-probes %s)\n",
              cs$name, status, res$n_solutions,
              ifelse(is.na(med_true), "NA", sprintf("%.2f", med_true)),
              res$bounds$median[1], res$bounds$median[2],
              sum(true$event), res$bounds$events[1], res$bounds$events[2], s_ok))
}
if (fails == 0L) cat("EXACT ROUND-TRIP PASSED\n") else
  stop(fails, " exact round-trip case(s) failed")
