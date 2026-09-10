# kmextract

Reconstruct pseudo individual-patient data (IPD) from a published Kaplan–Meier
curve image (PNG/JPEG), constrained by the printed number-at-risk table.

## Method

1. Automated pixel extraction of the KM step function (axis/tick calibration,
   curve tracing with connectivity rules; overlaid colored curves are separated
   by color).
2. Guyot-type reconstruction (Guyot et al., *BMC Med Res Methodol* 2012;12:9):
   deaths and censorings are allocated so that every printed at-risk /
   cumulative-event anchor is honored **exactly**; censoring times within
   intervals are placed uniformly.

## Usage

```r
source("R/km_extract.R")   # requires: png (or jpeg), survival

res <- km_extract(
  "my_km_figure.png",
  x_ticks = c(0, 6, 12, 18, 24),               # x-axis tick values
  anchors = data.frame(                        # the printed at-risk table
    t = c(0, 6, 12, 18, 24),                   #   time
    n = c(30, 21, 9, 4, 1),                    #   number at risk
    E = c(0, 10, 19, 22, 23)),                 #   CUMULATIVE events
  t_end = 24)                                  # where the curve visibly ends

res$ipd         # data.frame(time, event) — the reconstructed pseudo-IPD
res$lifetable   # survfit-style life table
km_check(res, published_median = 9.0)          # compare vs published values
```

For overlaid colored curves, call once per curve with `color = c(r, g, b)`
(0–1, read off the legend) and `seed = TRUE` if curves share the start at
S = 1. If the curve ends with a drop to zero *after* the last printed at-risk
column, append a final anchor row with `n = 0` and the extra event(s).

**No printed at-risk table?** Use `km_extract_free(image, x_ticks, N, y_ticks,
color =)` — deaths are read from the step drops and censoring times from the
"+" tick marks drawn on the curve; only the starting sample size `N` is
needed. This requires the censor ticks to be visible on the figure; without
them, prefer asking the data owner for the at-risk table. In our validation
this mode recovered a 20-patient curve's median within 0.1 months and its
event count exactly (internal validation on a real oncology figure).

Self-contained validation: `Rscript example/synthetic_test.R` simulates a
dataset, draws its KM to PNG, re-extracts it from the image and compares with
the known truth (median recovered within 0.01 months in the shipped example).

## Accuracy, from our validation

What the tool actually reconstructs well is the **survival path**: in our
synthetic round-trip tests (`example/`), the maximum deviation of the
reconstructed KM curve was 0.02–0.05 in survival probability, event counts
were exact, and at-risk anchors are honored exactly by construction. A
synthetic two-arm test recovered the Cox HR within 0.2%.

Point statistics derived from the path behave accordingly:

- **Medians** were recovered within ~0.4 months in all our round-trip tests
  (including a knife-edge case where the true curve grazes S = 0.503, solved
  by the fit-to-trace refinement of the death/censoring allocation). Keep in
  mind the median is a discontinuous functional — always cross-check
  reconstructed medians against the published ones before use.
- Overlaid curves: steps hidden under another curve near the median crossing
  can shift the reconstructed median by ~1 month.
- Reconstructed **median confidence intervals are unreliable** (discrete-jump
  statistic) — use the reconstruction for HR-type analyses, quote medians/CIs
  from the original publication.
- **Anchor-free mode**: each censor tick missed by the detector (overlapping
  ticks merge; low resolution hides them) becomes a spurious death. The
  function reports the detected tick count — verify it against the figure.

Run the shipped validation: `Rscript example/synthetic_test.R` and
`Rscript example/adversarial_tests.R` (heavy censoring, curve-to-zero,
N = 150, anchor-free round trips).

## Caveats

- This is a *reconstruction*, to be labeled as such in any output; it is not a
  substitute for source data. If the data owner can export an aggregated life
  table (`survfit`'s time / n.risk / n.event / n.censor), prefer that: it is
  the same information as the published curve, exact, and needs no image
  processing.
- Manual alternative with the same statistical engine: digitize the curve with
  WebPlotDigitizer and feed the coordinates to the `IPDfromKM` CRAN package
  (Liu N et al., *BMC Med Res Methodol* 2021).

## License

MIT — see LICENSE.
