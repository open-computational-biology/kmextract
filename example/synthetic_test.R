# Self-contained validation: simulate IPD -> draw the KM as a PNG -> re-extract
# from the image -> compare with the known truth. No real data involved.
# run from the repository root: Rscript example/synthetic_test.R
source("R/km_extract.R")
library(survival)

set.seed(42)
n <- 30
true <- data.frame(time = round(rexp(n, log(2) / 9), 2),           # median ~9 mo
                   event = rbinom(n, 1, 0.75))
true$time <- pmin(true$time, 24)                                    # admin cut at 24
true$event[true$time >= 24] <- 0

fit <- survfit(Surv(time, event) ~ 1, data = true)

png("example/synthetic_km.png", width = 1100, height = 700)
par(mar = c(5, 5, 2, 2), lwd = 2)
plot(fit, conf.int = FALSE, mark.time = TRUE, lwd = 3,
     xlab = "Time (months)", ylab = "Survival probability",
     xaxt = "n", yaxt = "n", bty = "l", xlim = c(0, 24))
axis(1, at = seq(0, 24, 6), lwd = 2, lwd.ticks = 2)
axis(2, at = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25)), lwd = 2, lwd.ticks = 2)
dev.off()

# anchors as they would be printed under the plot
grid <- seq(0, 24, 6)
anchors <- data.frame(
  t = grid,
  n = sapply(grid, function(tt) sum(true$time >= tt)),
  E = sapply(grid, function(tt) sum(true$event == 1 & true$time < tt)))
t_end <- max(true$time)
# events after the last anchor? (final drop-to-zero rule)
E_after <- sum(true$event == 1 & true$time >= max(grid))
if (E_after > 0) anchors <- rbind(anchors, data.frame(t = t_end, n = 0, E = max(anchors$E) + E_after))

res <- km_extract("example/synthetic_km.png", x_ticks = seq(0, 24, 6),
                  anchors = anchors, t_end = t_end)

s_true <- summary(survfit(Surv(time, event) ~ 1, data = true))$table
s_rec  <- summary(res$fit)$table
cat(sprintf("TRUE : N=%d ev=%d median %.2f (%.2f-%s)\n", s_true[["records"]], s_true[["events"]],
    s_true[["median"]], s_true[["0.95LCL"]],
    ifelse(is.na(s_true[["0.95UCL"]]), "NR", sprintf("%.2f", s_true[["0.95UCL"]]))))
cat(sprintf("RECON: N=%d ev=%d median %.2f (%.2f-%s)\n", s_rec[["records"]], s_rec[["events"]],
    s_rec[["median"]], s_rec[["0.95LCL"]],
    ifelse(is.na(s_rec[["0.95UCL"]]), "NR", sprintf("%.2f", s_rec[["0.95UCL"]]))))
stopifnot(s_rec[["records"]] == s_true[["records"]],
          s_rec[["events"]] == s_true[["events"]],
          abs(s_rec[["median"]] - s_true[["median"]]) < 0.5)
cat("SELF-TEST PASSED\n")
