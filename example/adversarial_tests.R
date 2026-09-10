# Adversarial round-trip tests (added after independent review, 2026-09-10):
# heavy interleaved censoring, curve to zero, large N, and the anchor-free mode.
# Run from the repository root: Rscript example/adversarial_tests.R
source("R/km_extract.R")
library(survival)

run_case <- function(true, grid, label, N = nrow(true), free = FALSE) {
  fit <- survfit(Surv(time, event) ~ 1, data = true)
  f <- tempfile(fileext = ".png")
  png(f, width = 1100, height = 700)
  par(mar = c(5, 5, 2, 2), lwd = 2)
  plot(fit, conf.int = FALSE, mark.time = TRUE, lwd = 3, xaxt = "n", yaxt = "n",
       bty = "l", xlim = c(0, max(grid)))
  axis(1, at = grid, lwd = 2, lwd.ticks = 2)
  axis(2, at = seq(0, 1, 0.25), lwd = 2, lwd.ticks = 2)
  dev.off()
  if (free) {
    res <- suppressMessages(km_extract_free(f, x_ticks = grid, N = N,
                                            y_ticks = seq(1, 0, -0.25)))
  } else {
    anchors <- data.frame(t = grid,
      n = sapply(grid, function(tt) sum(true$time >= tt)),
      E = sapply(grid, function(tt) sum(true$event == 1 & true$time < tt)))
    t_end <- max(true$time)
    E_after <- sum(true$event == 1 & true$time >= max(grid))
    if (E_after > 0) anchors <- rbind(anchors, data.frame(t = t_end, n = 0, E = max(anchors$E) + E_after))
    res <- km_extract(f, x_ticks = grid, anchors = anchors, t_end = t_end,
                      y_ticks = seq(1, 0, -0.25))
  }
  s_t <- summary(fit)$table; s_r <- summary(res$fit)$table
  # primary metric: maximum deviation of the reconstructed KM PATH (the median
  # is a discontinuous functional: when the true curve grazes 0.5, it can jump
  # to an adjacent event time even for a near-perfect reconstruction)
  sf_t <- stepfun(fit$time, c(1, fit$surv)); sf_r <- stepfun(res$fit$time, c(1, res$fit$surv))
  gr <- seq(0.05, min(max(true$time), max(res$ipd$time)) - 0.05, length.out = 400)
  path_err <- max(abs(sf_t(gr) - sf_r(gr)))
  cat(sprintf("%-28s median %6.2f vs %6.2f (d %+5.2f) | ev %d/%d | max path error %.3f\n",
      label, s_r[["median"]], s_t[["median"]], s_r[["median"]] - s_t[["median"]],
      s_r[["events"]], s_t[["events"]], path_err))
  stopifnot(s_r[["records"]] == nrow(true), path_err < 0.06)
  invisible(res)
}

set.seed(7)
# T1: heavy censoring interleaved with deaths early
n <- 40
t1 <- data.frame(time = round(c(runif(13, 0.5, 6), rexp(27, log(2)/9)), 2), event = c(rep(0,13), rep(1,27)))
t1$time <- pmin(t1$time, 24); t1$event[t1$time >= 24] <- 0
run_case(t1, seq(0, 24, 6), "T1 heavy early censoring")

# T2: everyone dies (curve to zero)
t2 <- data.frame(time = round(rexp(25, log(2)/8), 2), event = 1)
run_case(t2, seq(0, ceiling(max(t2$time)/6)*6, 6), "T2 curve to zero")

# T3: large N, small steps
t3 <- data.frame(time = round(rexp(150, log(2)/14), 2), event = rbinom(150, 1, 0.75))
t3$time <- pmin(t3$time, 36); t3$event[t3$time >= 36] <- 0
run_case(t3, seq(0, 36, 6), "T3 large N (150)")

# T5: anchor-free mode with censor ticks
set.seed(11)
t5 <- data.frame(time = round(rexp(20, log(2)/10), 2), event = rbinom(20, 1, 0.65))
t5$time <- pmin(t5$time, 20); t5$event[t5$time >= 20] <- 0
run_case(t5, seq(0, 20, 4), "T5 anchor-free", free = TRUE)

cat("ADVERSARIAL TESTS PASSED\n")
