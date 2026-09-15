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

## Exact mode (`km_extract_exact`)

For figures that print an at-risk (cumulative events) table, exact mode
replaces heuristic reconstruction with **constraint solving**: it extracts
everything the figure determines exactly, and enumerates the complete set of
datasets compatible with the figure where it does not.

```r
source("R/km_extract.R")
source("R/km_extract_exact.R")
res <- km_extract_exact(
  "panel.png", x_ticks = c(0, 10, 20, 30),
  anchors = data.frame(t = c(0, 10, 20, 30), n = c(60, 14, 2, 1),
                       E = c(0, 20, 26, 26)),
  # optional: printed summary values as exact (rounded-interval) constraints
  known = list(surv = data.frame(time = 12, est = 44.0, lcl = 26.1, ucl = 60.5),
               median = c(10.6, 7.2, 18.3), conf.type = "log-log"),
  # optional: a comparator arm (time, event) to bound the Cox HR over the set
  comparator = my_arm)
res$unique        # TRUE when the figure pins the dataset completely
res$n_solutions   # size of the enumerated admissible set
res$bounds        # certified min/max per statistic over the set
res$ipd           # canonical member (median-HR when a comparator is given)
res$certificate   # per-plateau pixel deviations of the canonical solution
km_overlay_exact(res, "panel.png", "overlay.png")   # human-readable certificate
```

How it works: (1) plateaus are re-measured on clean line columns and risers
located on pure stroke columns; (2) censor marks are censused by three
complementary channels — symmetric protrusion on plateaus, horizontal arm
signatures at plateau edges, floating-bar residues at risers — each mark
assigned to its inter-death gap by its LEVEL (which also resolves marks merged
into a riser stroke); (3) a depth-first exact solver enumerates every
(deaths-per-riser, censors-per-gap) sequence satisfying the anchors exactly,
every measured plateau to pixel tolerance (with an affine allowance for
label-centroid calibration drift), and at least the visible marks per gap;
(4) printed summary values, when supplied, further filter the set as
rounded-interval constraints.

The contract, honestly stated: counts and orderings — which determine every
KM statistic — are recovered exactly whenever the figure determines them
(`unique = TRUE`); otherwise `bounds` gives the exact range of any statistic
over all admissible datasets. Times are pixel-limited (~half a line width),
which no statistic feels. Validation: `Rscript example/exact_roundtrip.R`
(bold curves, no axis ticks, heavy ties, censors hugging deaths; the truth
must be covered by the certified bounds — and is, including a case where the
solution is unique and equals the truth).

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
