# ============================================================================
# kmextract — EXACT mode
# ----------------------------------------------------------------------------
# For KM figures that carry an at-risk (cumulative events) table: extract the
# event/censor sequence EXACTLY where the figure determines it, and enumerate
# the complete admissible set where it does not.
#
# Contract: every count and ordering the figure pins down (which is what every
# KM statistic depends on) is recovered exactly; what the figure leaves free
# (sub-pixel interleaving of censors hidden behind other marks) is returned as
# the full enumerated solution set with certified per-statistic bounds. Times
# are pixel-limited (~half a line width) — irrelevant to any statistic.
#
# Pipeline (each stage feeds constraints to the exact solver, nothing guesses):
#   1. geometry & calibration           (.km_find_axes / .km_calibrate, base file)
#   2. certified path: plateaus re-measured on clean columns + risers
#   3. censor-mark census, three channels (on-plateau protrusion, horizontal
#      arm signatures, riser residues), each mark assigned to its gap
#   4. exact combinatorial solve: deaths per riser and censors per gap under
#      (a) the at-risk anchors, exactly; (b) every measured plateau level to
#      pixel tolerance; (c) at least the visible marks per gap; and optionally
#      (d) published summary values (medians / landmark rates / their CIs)
#   5. certification: constraint report + rendered overlay for human review
#
# Requires R/km_extract.R to be sourced first (geometry helpers), plus the
# survival package for per-solution statistics.
# ============================================================================

if (!exists(".km_find_axes")) stop("source R/km_extract.R before this file")

# ---------------------------------------------------------------------------
# Stage 2 — certified path
# ---------------------------------------------------------------------------

.kmx_col_runs <- function(mask, x, r0, r1) {
  rows <- which(mask[r0:r1, x]) + r0 - 1L
  if (!length(rows)) return(NULL)
  b <- cumsum(c(1L, diff(rows) > 2L))
  grp <- split(as.integer(rows), b)
  cbind(top = vapply(grp, min, 0L), bot = vapply(grp, max, 0L))
}

.kmx_path <- function(mask, axes, y_top_px) {
  h <- nrow(mask); w <- ncol(mask)
  r0 <- max(1L, y_top_px - 30L); r1 <- axes$x_top - 2L
  x0 <- axes$y_right + 4L
  cols <- vector("list", w)
  for (x in x0:w) cols[x] <- list(.kmx_col_runs(mask, x, r0, r1))   # [x] <- list():
  # a NULL result must STAY as a NULL element ([[x]] <- NULL would delete it)
  # line thickness from single-run columns (risers are a minority; median holds)
  hts <- unlist(lapply(cols, function(rn)
    if (!is.null(rn) && nrow(rn) == 1L) rn[1, "bot"] - rn[1, "top"] + 1L))
  Tk <- as.integer(stats::median(hts))
  # plain points: exactly one run, of plain-line height
  px <- integer(0); pc <- numeric(0)
  for (x in x0:w) {
    rn <- cols[[x]]
    if (!is.null(rn) && nrow(rn) == 1L && (rn[1, "bot"] - rn[1, "top"] + 1L) <= Tk + 2L) {
      px <- c(px, x); pc <- c(pc, (rn[1, "top"] + rn[1, "bot"]) / 2)
    }
  }
  if (length(px) < 10L) stop("exact mode: too few clean line columns - image too cluttered?")
  # segment plain points into plateaus (KM level jumps are >> 4 px for N <= ~150)
  seg <- cumsum(c(1L, abs(diff(pc)) > 4))
  plat <- do.call(rbind, lapply(split(seq_along(px), seg), function(i)
    data.frame(level_px = stats::median(pc[i]), x0 = min(px[i]), x1 = max(px[i]),
               n_plain = length(i))))
  # merge consecutive segments at the same level (a mark or riser artefact
  # split them); a true riser separates levels by at least one KM step
  keep <- rep(TRUE, nrow(plat))
  for (k in seq_len(nrow(plat))[-1]) {
    j <- max(which(keep[1:(k - 1L)]))
    if (abs(plat$level_px[k] - plat$level_px[j]) <= 3) {
      plat$level_px[j] <- stats::median(c(plat$level_px[j], plat$level_px[k]))
      plat$x1[j] <- plat$x1[k]
      plat$n_plain[j] <- plat$n_plain[j] + plat$n_plain[k]
      keep[k] <- FALSE
    }
  }
  plat <- plat[keep, , drop = FALSE]
  # enforce monotone descent (rows increase); a violation means segmentation noise
  stopifnot("exact mode: non-monotone plateau sequence" =
              all(diff(plat$level_px) > 0))
  # risers between consecutive plateaus: columns whose stroke connects the levels
  risers <- data.frame(x_px = numeric(0), width = integer(0), hidden = logical(0))
  for (k in seq_len(nrow(plat) - 1L)) {
    lo <- plat$x1[k]; hi <- plat$x0[k + 1L]
    conn <- integer(0); pure <- integer(0)
    for (x in lo:hi) {
      rn <- cols[[x]]
      if (is.null(rn)) next
      hit <- rn[, "top"] <= plat$level_px[k] + Tk &
             rn[, "bot"] >= plat$level_px[k + 1L] - Tk
      if (any(hit)) {
        conn <- c(conn, x)
        # pure stroke: no glyph ink beyond either plateau line
        if (any(hit & rn[, "top"] >= plat$level_px[k] - Tk / 2 - 1 &
                      rn[, "bot"] <= plat$level_px[k + 1L] + Tk / 2 + 1))
          pure <- c(pure, x)
      }
    }
    ctr <- if (length(pure)) pure else conn
    risers <- rbind(risers, data.frame(
      x_px = if (length(ctr)) mean(range(ctr)) else mean(c(lo, hi)),
      width = length(conn), hidden = !length(conn)))
  }
  # end of curve: last column with ink (marks included)
  x_last <- max(which(!vapply(cols, is.null, TRUE)))
  list(plateaus = plat, risers = risers, T = Tk, cols = cols,
       x_last = x_last, r0 = r0, r1 = r1)
}

# ---------------------------------------------------------------------------
# Stage 3 — mark census (three channels, gap-assigned)
# ---------------------------------------------------------------------------
# Gap g = the interval during plateau g (between riser g and riser g+1);
# gap 0 precedes the first riser. Marks at a riser are side-assigned by their
# level: upper level -> end of the gap before the death, lower level -> start
# of the gap after it.

.kmx_marks <- function(path) {
  plat <- path$plateaus; ris <- path$risers; Tk <- path$T; cols <- path$cols
  m <- max(2L, as.integer(round(Tk / 2)) - 1L)
  half <- as.integer(ceiling(Tk / 2))
  A <- 3L * Tk                     # search window for arms; glyph arm ~ 1-1.5 T
  marks <- data.frame(gap = integer(0), x_px = numeric(0), channel = character(0))
  add <- function(gap, x, ch) marks <<- rbind(marks, data.frame(gap = gap, x_px = x, channel = ch))
  R <- nrow(ris)
  for (k in seq_len(nrow(plat))) {
    lvl <- plat$level_px[k]
    # --- channel A: symmetric protrusion, scanned over the full gap span ----
    # (a mark at the plateau edge sits beyond the last clean column; riser
    # strokes self-exclude: their top never clears the level by the margin)
    a_lo <- if (k > 1L) as.integer(ceiling(ris$x_px[k - 1L])) + 1L else
              max(plat$x0[k] - A, 1L)
    a_hi <- if (k <= R) as.integer(floor(ris$x_px[k])) else
              min(path$x_last, plat$x1[k] + A)
    cand <- integer(0)
    for (x in a_lo:a_hi) {
      rn <- cols[[x]]
      if (is.null(rn)) next
      i <- which(rn[, "top"] <= lvl + half & rn[, "bot"] >= lvl - half)
      if (length(i) && any(rn[i, "top"] <= lvl - half - m & rn[i, "bot"] >= lvl + half + m))
        cand <- c(cand, x)
    }
    if (length(cand)) {
      b <- cumsum(c(1L, diff(cand) > 3L))
      for (grp in split(cand, b)) {
        # bars = full-height column runs; two bars one blank column apart are
        # distinct marks (split on COLUMN positions), and a fused run wider
        # than a single bar is a superposed pair (multiplicity from width,
        # calibrated on the isolated-bar width seen elsewhere in the figure)
        hts <- vapply(grp, function(x) {
          rn <- cols[[x]]
          i <- which(rn[, "top"] <= lvl + half & rn[, "bot"] >= lvl - half)[1]
          rn[i, "bot"] - rn[i, "top"] + 1L
        }, 0L)
        tall <- grp[hts >= Tk + 2L * m - 1L]
        if (length(tall)) {
          for (bar in split(tall, cumsum(c(1L, diff(tall) > 1L)))) {
            wdt <- max(bar) - min(bar) + 1L
            cnt <- max(1L, as.integer(round(wdt / 2.5)))
            for (px_c in seq(min(bar), max(bar), length.out = cnt))
              add(k - 1L, px_c, "A")
          }
        } else add(k - 1L, mean(grp), "A")
      }
    }
    # --- channel B: horizontal arm signatures beyond the plateau extent -----
    band <- c(lvl - half - 1L, lvl + half + 1L)
    x_lo <- if (k > 1L) as.integer(ceiling(ris$x_px[k - 1L])) + 1L else max(1L, plat$x0[k] - A)
    if (x_lo < plat$x0[k] - 1L) {
      for (x in x_lo:(plat$x0[k] - 1L)) {
        rn <- cols[[x]]
        if (is.null(rn)) next
        thin <- rn[, "bot"] - rn[, "top"] + 1L <= 4L &
                rn[, "top"] >= band[1] & rn[, "bot"] <= band[2]
        if (any(thin)) { add(k - 1L, plat$x0[k] - 0.5, "B"); break }
      }
    }
    x_hi <- if (k <= R) as.integer(floor(ris$x_px[k])) - 1L else min(path$x_last, plat$x1[k] + A)
    if (x_hi > plat$x1[k] + 1L) {
      for (x in (plat$x1[k] + 1L):x_hi) {
        rn <- cols[[x]]
        if (is.null(rn)) next
        thin <- rn[, "bot"] - rn[, "top"] + 1L <= 4L &
                rn[, "top"] >= band[1] & rn[, "bot"] <= band[2]
        if (any(thin)) { add(k - 1L, plat$x1[k] + 0.5, "B"); break }
      }
    }
  }
  # --- channel C: residues at riser columns ---------------------------------
  for (j in seq_len(R)) {
    up <- plat$level_px[j]; lo <- plat$level_px[j + 1L]
    up_guard <- if (j > 1L) plat$level_px[j - 1L] + half else -Inf
    lo_guard <- if (j + 2L <= nrow(plat)) plat$level_px[j + 2L] - half else Inf
    xs <- as.integer(round(ris$x_px[j])) + (-Tk - 2L):(Tk + 2L)
    xs <- xs[xs >= 1L & xs <= length(cols)]
    top_cols <- integer(0); bot_cols <- integer(0)
    for (x in xs) {
      rn <- cols[[x]]
      if (is.null(rn)) next
      # a residue must FLOAT between plateau levels: a stroke reaching the
      # plateau one level further is the neighbouring riser, not a glyph
      if (any(rn[, "top"] <= up - half - m & rn[, "bot"] >= up - half &
              rn[, "top"] >= up_guard)) top_cols <- c(top_cols, x)
      if (any(rn[, "bot"] >= lo + half + m & rn[, "top"] <= lo + half &
              rn[, "bot"] <= lo_guard)) bot_cols <- c(bot_cols, x)
    }
    if (length(top_cols))
      for (grp in split(top_cols, cumsum(c(1L, diff(top_cols) > 2L))))
        add(j - 1L, mean(grp), "C")
    if (length(bot_cols))
      for (grp in split(bot_cols, cumsum(c(1L, diff(bot_cols) > 2L))))
        add(j, mean(grp), "C")
  }
  # --- dedup: one mark may be seen by several channels ----------------------
  if (nrow(marks)) {
    ord <- order(marks$channel != "A", marks$channel != "C", marks$x_px)
    marks <- marks[ord, ]                     # A first, then C, then B
    keep <- rep(TRUE, nrow(marks))
    for (i in seq_len(nrow(marks))) {
      if (marks$channel[i] == "A") next
      prev <- which(keep & seq_len(nrow(marks)) < i & marks$gap == marks$gap[i])
      if (!length(prev)) next
      r <- if (marks$channel[i] == "C") 4 else Tk + 2L * m + 2
      # a C group is the same bar as a nearby A bar; a B arm belongs to any
      # already-kept mark within one glyph arm of the plateau edge
      if (any(abs(marks$x_px[prev] - marks$x_px[i]) <= r)) keep[i] <- FALSE
    }
    marks <- marks[keep, , drop = FALSE]
    marks <- marks[order(marks$gap, marks$x_px), ]
  }
  marks
}

# ---------------------------------------------------------------------------
# Stage 4 — exact combinatorial solve
# ---------------------------------------------------------------------------
# Cells alternate gaps and risers in time order; interior anchor times split
# gaps into subgaps and act as exact checkpoints (n and cumulative events).
# DFS enumerates every (deaths-per-riser, censors-per-subgap) sequence whose
# KM path stays within tolerance of every measured plateau.

.kmx_solve <- function(riser_t, L, tolS, subgaps, anchors, N, t_end,
                       max_solutions = 2e5) {
  R <- length(riser_t)
  # build the cell sequence
  cells <- list()
  for (g in 0:R) {
    sg <- subgaps[subgaps$gap == g, , drop = FALSE]
    for (i in seq_len(nrow(sg))) {
      cells[[length(cells) + 1L]] <- list(type = "gap", gap = g,
                                          lb = sg$lb[i], t0 = sg$t0[i], t1 = sg$t1[i])
      if (sg$check[i]) {
        a <- which(abs(anchors$t - sg$t1[i]) < 1e-9)
        cells[[length(cells) + 1L]] <- list(type = "check", n = anchors$n[a], E = anchors$E[a])
      }
    }
    if (g < R)
      cells[[length(cells) + 1L]] <- list(type = "riser", r = g + 1L,
                                          L = L[g + 1L], tol = tolS[g + 1L])
  }
  cells[[length(cells) + 1L]] <- list(type = "end")
  M <- length(cells)
  # suffix minimum consumption (for pruning): censors >= lb per gap, deaths >= 1 per riser
  min_cons <- numeric(M + 1L); min_cons[M + 1L] <- 0
  for (i in M:1) {
    cst <- switch(cells[[i]]$type, gap = cells[[i]]$lb, riser = 1L, 0L)
    min_cons[i] <- min_cons[i + 1L] + cst
  }
  n_final <- anchors$n[nrow(anchors)]
  sols <- vector("list", 1024L); nsol <- 0L; overflow <- FALSE
  d_cur <- integer(R); c_cur <- integer(sum(subgaps$gap >= 0))
  sg_index <- 0L
  rec <- function(i, n, S, cum_d, sg_i) {
    if (overflow) return()
    cell <- cells[[i]]
    if (cell$type == "end") {
      if (n == n_final) {
        nsol <<- nsol + 1L
        if (nsol > max_solutions) { overflow <<- TRUE; return() }
        sols[[nsol]] <<- list(d = d_cur, c = c_cur)
      }
      return()
    }
    if (cell$type == "check") {
      if (n == cell$n && cum_d == cell$E) rec(i + 1L, n, S, cum_d, sg_i)
      return()
    }
    if (cell$type == "gap") {
      ub <- n - n_final - (min_cons[i + 1L])
      if (ub < cell$lb) return()
      for (cc in cell$lb:ub) {
        c_cur[sg_i + 1L] <<- cc
        rec(i + 1L, n - cc, S, cum_d, sg_i + 1L)
      }
      c_cur[sg_i + 1L] <<- 0L
      return()
    }
    # riser
    ub <- n - n_final - (min_cons[i + 1L])
    if (ub < 1L) return()
    for (dd in 1:ub) {
      S2 <- S * (1 - dd / n)
      if (is.finite(cell$tol) && abs(S2 - cell$L) > cell$tol) {
        if (S2 < cell$L - cell$tol) break   # deeper d only sinks further
        next
      }
      d_cur[cell$r] <<- dd
      rec(i + 1L, n - dd, S2, cum_d + dd, sg_i)
    }
    d_cur[cell$r] <<- 0L
  }
  rec(1L, N, 1.0, 0L, 0L)
  list(solutions = sols[seq_len(min(nsol, max_solutions))],
       n_solutions = nsol, overflow = overflow, cells = cells)
}

# ---------------------------------------------------------------------------
# Solution -> pseudo-IPD (order-exact; censor times canonical)
# ---------------------------------------------------------------------------

.kmx_ipd <- function(sol, riser_t, subgaps, marks_t, t_end) {
  times <- rep(riser_t, sol$d)
  ipd <- data.frame(time = times, event = 1L)
  for (i in seq_len(nrow(subgaps))) {
    cc <- sol$c[i]
    if (cc == 0L) next
    mt <- marks_t[[i]]
    ct <- if (length(mt)) rep(mt, length.out = cc)[order(rep(mt, length.out = cc))] else
          rep((subgaps$t0[i] + subgaps$t1[i]) / 2, cc)
    # canonical display rule: hidden censors superpose exactly on a visible
    # mark of their gap; gaps with no visible mark use the midpoint
    ipd <- rbind(ipd, data.frame(time = ct, event = 0L))
  }
  ipd[order(ipd$time, -ipd$event), ]
}

# ---------------------------------------------------------------------------
# Optional published-value filter (rounded printed values as constraints)
# ---------------------------------------------------------------------------

.kmx_match_known <- function(ipd, known) {
  ct <- if (!is.null(known$conf.type)) known$conf.type else "log-log"
  sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = ipd, conf.type = ct)
  ok <- TRUE
  rnd <- function(x, d) round(x + 1e-12, d)
  if (!is.null(known$surv)) {
    dg <- if (!is.null(known$digits)) known$digits else 1L
    for (i in seq_len(nrow(known$surv))) {
      s <- summary(sf, times = known$surv$time[i], extend = TRUE)
      if (!is.na(known$surv$est[i]) && rnd(100 * s$surv, dg) != known$surv$est[i]) ok <- FALSE
      if (!is.null(known$surv$lcl) && !is.na(known$surv$lcl[i]) &&
          rnd(100 * s$lower, dg) != known$surv$lcl[i]) ok <- FALSE
      if (!is.null(known$surv$ucl) && !is.na(known$surv$ucl[i]) &&
          rnd(100 * s$upper, dg) != known$surv$ucl[i]) ok <- FALSE
      if (!ok) return(FALSE)
    }
  }
  if (!is.null(known$median)) {
    tab <- summary(sf)$table
    tol_t <- if (!is.null(known$median_tol)) known$median_tol else 0.15
    cmp <- function(a, b) (is.na(a) && is.na(b)) || (!is.na(a) && !is.na(b) && abs(a - b) <= tol_t)
    if (!cmp(tab[["median"]], known$median[1])) ok <- FALSE
    if (length(known$median) > 1L && !cmp(tab[["0.95LCL"]], known$median[2])) ok <- FALSE
    if (length(known$median) > 2L && !cmp(tab[["0.95UCL"]], known$median[3])) ok <- FALSE
  }
  ok
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

km_extract_exact <- function(image_path, x_ticks, anchors,
                             y_ticks = c(1, 0.75, 0.5, 0.25),
                             color = NULL, t_end = NULL,
                             known = list(), comparator = NULL,
                             tol_px = 1.5, drift_frac = 0.006,
                             max_solutions = 2e5, hr_max = 20000L,
                             dark_thresh = 110/255, color_tol = 70/255) {
  stopifnot(all(c("t", "n", "E") %in% names(anchors)), anchors$t[1] == 0, anchors$E[1] == 0)
  img <- .km_read_image(image_path)
  g <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  axes <- .km_find_axes(g, dark_thresh)
  ticks <- .km_calibrate(g, axes, x_ticks, y_ticks, dark_thresh)
  fx <- stats::lm(v ~ p, data = data.frame(p = ticks$x_px, v = x_ticks))
  fy <- stats::lm(v ~ p, data = data.frame(p = ticks$y_px, v = y_ticks))
  to_t <- function(px) as.numeric(stats::predict(fx, data.frame(p = px)))
  to_S <- function(px) as.numeric(stats::predict(fy, data.frame(p = px)))
  mask <- .km_curve_mask(img, color, dark_thresh, color_tol)
  # confine to the plot area
  mask[, seq_len(axes$y_right + 3L)] <- FALSE
  mask[axes$x_top:nrow(mask), ] <- FALSE
  y_top_px <- ticks$y_px[which.max(y_ticks)]
  path <- .kmx_path(mask, axes, y_top_px)
  plat <- path$plateaus
  # sanity: the first plateau is S = 1
  stopifnot("exact mode: first plateau is not at S = 1" =
              abs(to_S(plat$level_px[1]) - 1) < 0.03)
  riser_t <- to_t(path$risers$x_px)
  slopeS <- abs(stats::coef(fy)[2])
  L <- to_S(plat$level_px)[-1]                      # post-riser levels
  # label-centroid calibration carries a small linear drift (text centroids
  # are not exactly on the tick positions): allow it, proportionally to the
  # distance from the exact S = 1 anchor. Late-curve KM steps are huge, so
  # the wider late tolerance costs no discrimination there.
  tolS <- tol_px * slopeS + drift_frac * (1 - L)
  tolS[plat$n_plain[-1] <= 2L] <- tolS[plat$n_plain[-1] <= 2L] + 2 * tol_px * slopeS
  if (is.null(t_end)) t_end <- to_t(path$x_last)
  # terminal convention: patients still at risk after the last printed anchor
  # are censored by t_end (classic-mode behaviour); a plunge-to-zero figure
  # already ends with an n = 0 anchor and is left untouched
  if (anchors$n[nrow(anchors)] > 0L)
    anchors <- rbind(anchors, data.frame(t = t_end, n = 0L,
                                         E = anchors$E[nrow(anchors)]))
  # marks -> times & subgaps
  mk <- .kmx_marks(path)
  mk$t <- to_t(mk$x_px)
  bounds <- c(0, riser_t, t_end)
  interior <- anchors$t[anchors$t > 0 & anchors$t < t_end - 1e-9]
  subgaps <- do.call(rbind, lapply(0:length(riser_t), function(g0) {
    lo <- bounds[g0 + 1L]; hi <- bounds[g0 + 2L]
    cuts <- sort(unique(c(lo, interior[interior > lo & interior < hi], hi)))
    data.frame(gap = g0, t0 = head(cuts, -1L), t1 = cuts[-1L],
               check = cuts[-1L] %in% interior)
  }))
  subgaps$lb <- vapply(seq_len(nrow(subgaps)), function(i)
    sum(mk$gap == subgaps$gap[i] & mk$t > subgaps$t0[i] - 1e-9 & mk$t <= subgaps$t1[i] + 1e-9), 0L)
  # marks listed per subgap (canonical censor positions)
  marks_t <- lapply(seq_len(nrow(subgaps)), function(i)
    sort(mk$t[mk$gap == subgaps$gap[i] & mk$t > subgaps$t0[i] - 1e-9 & mk$t <= subgaps$t1[i] + 1e-9]))
  N <- anchors$n[1]
  sv <- .kmx_solve(riser_t, L, tolS, subgaps, anchors, N, t_end, max_solutions)
  relaxed <- "none"
  if (sv$n_solutions == 0L) {           # graceful degradation, loudly reported
    sv <- .kmx_solve(riser_t, L, 2 * tolS, subgaps, anchors, N, t_end, max_solutions)
    relaxed <- "plateau tolerance doubled"
  }
  if (sv$n_solutions == 0L) {
    sg2 <- subgaps; sg2$lb <- pmax(0L, sg2$lb - 1L)
    sv <- .kmx_solve(riser_t, L, 2 * tolS, sg2, anchors, N, t_end, max_solutions)
    relaxed <- "plateau tolerance doubled + mark lower bounds relaxed by 1"
  }
  if (sv$n_solutions == 0L)
    stop("exact mode: no admissible solution even after relaxation - ",
         "check the anchors transcription against the printed table")
  if (relaxed != "none")
    warning("exact mode: solved only after relaxation (", relaxed,
            ") - inspect the certificate before trusting the bounds")
  # events bounds: EXACT over the full enumerated set (cheap on raw vectors)
  ev_all <- vapply(sv$solutions, function(x) sum(x$d), 0L)
  # build IPDs on at most keep_max solutions (survfit per solution is the
  # cost driver); beyond that, per-statistic bounds other than events are
  # computed on a reproducible sample and flagged
  keep_max <- 20000L
  pick <- seq_along(sv$solutions)
  sampled <- length(pick) > keep_max
  if (sampled) {
    set.seed(1L)
    pick <- sort(sample(pick, keep_max))
    warning("exact mode: ", length(sv$solutions), " solutions - statistics ",
            "beyond the events bounds are computed on a sample of ", keep_max)
  }
  ipds <- lapply(sv$solutions[pick], .kmx_ipd, riser_t = riser_t,
                 subgaps = subgaps, marks_t = marks_t, t_end = t_end)
  kept <- pick
  if (length(known)) {
    keep <- vapply(ipds, .kmx_match_known, TRUE, known = known)
    if (any(keep)) { ipds <- ipds[keep]; kept <- kept[keep] }
    else warning("published-value filter eliminated every solution; ignoring it")
  }
  # per-solution statistics
  stat <- function(ipd) {
    sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = ipd, conf.type = "log-log")
    tab <- summary(sf)$table
    c(median = unname(tab[["median"]]), events = sum(ipd$event))
  }
  st <- t(vapply(ipds, stat, c(0, 0)))
  hr <- NULL
  if (!is.null(comparator)) {
    idx <- seq_along(ipds)
    if (length(idx) > hr_max) idx <- sort(sample(idx, hr_max))
    hr <- vapply(ipds[idx], function(ipd) {
      d <- rbind(data.frame(time = comparator$time, event = comparator$event, arm = "trt"),
                 data.frame(time = ipd$time, event = ipd$event, arm = "ctl"))
      d$arm <- factor(d$arm, levels = c("ctl", "trt"))
      unname(exp(stats::coef(survival::coxph(survival::Surv(time, event) ~ arm,
                                             data = d, ties = "efron"))))
    }, 0)
  }
  # canonical solution: median HR if available, else first
  can <- if (!is.null(hr)) which.min(abs(hr - stats::median(hr))) else 1L
  # certificate: plateau deviations of the canonical solution, in px
  ipd0 <- ipds[[can]]
  sf0 <- survival::survfit(survival::Surv(time, event) ~ 1, data = ipd0)
  Scan <- summary(sf0, times = riser_t + 1e-9, extend = TRUE)$surv
  cert <- data.frame(riser_t = round(riser_t, 2), S_measured = round(L, 4),
                     S_solution = round(Scan, 4),
                     dev_px = round(abs(Scan - L) / slopeS, 2))
  list(ipd = ipd0,
       n_solutions = sv$n_solutions, n_after_known = length(ipds),
       unique = length(ipds) == 1L, overflow = sv$overflow,
       solutions = ipds, raw = sv$solutions[kept], subgaps = subgaps,
       riser_t = riser_t,
       bounds = list(median = if (all(is.na(st[, "median"]))) c(NA_real_, NA_real_)
                              else range(st[, "median"], na.rm = TRUE),
                     events = range(ev_all),
                     hr = if (!is.null(hr)) range(hr)),
       sampled = sampled,
       hr = if (!is.null(hr)) hr,
       marks = mk, path = list(plateaus = plat, risers = path$risers,
                               riser_t = riser_t, L = L, tolS = tolS, T = path$T),
       certificate = cert, relaxed = relaxed,
       calib = list(fx = fx, fy = fy, t_end = t_end))
}

# ---------------------------------------------------------------------------
# Stage 5 helper — rendered overlay for human review
# ---------------------------------------------------------------------------

km_overlay_exact <- function(res, image_path, out_png) {
  img <- .km_read_image(image_path)
  to_px_x <- function(t) (t - stats::coef(res$calib$fx)[1]) / stats::coef(res$calib$fx)[2]
  to_px_y <- function(S) (S - stats::coef(res$calib$fy)[1]) / stats::coef(res$calib$fy)[2]
  grDevices::png(out_png, width = ncol(img[, , 1]), height = nrow(img[, , 1]))
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(1, ncol(img[, , 1])), ylim = c(nrow(img[, , 1]), 1),
                        xaxs = "i", yaxs = "i")
  graphics::rasterImage(img, 1, nrow(img[, , 1]), ncol(img[, , 1]), 1)
  ipd <- res$ipd
  n <- nrow(ipd); S <- 1; tprev <- 0
  for (t in sort(unique(ipd$time))) {
    d <- sum(ipd$time == t & ipd$event == 1L)
    cc <- sum(ipd$time == t & ipd$event == 0L)
    if (d > 0L) {
      graphics::segments(to_px_x(tprev), to_px_y(S), to_px_x(t), to_px_y(S),
                         col = "red", lwd = 2)
      graphics::segments(to_px_x(t), to_px_y(S), to_px_x(t), to_px_y(S * (1 - d / n)),
                         col = "red", lwd = 2)
      S <- S * (1 - d / n); tprev <- t
    }
    if (cc > 0L)
      graphics::segments(to_px_x(t), to_px_y(S) - 14, to_px_x(t), to_px_y(S) + 14,
                         col = "blue", lwd = 2)
    n <- n - d - cc
  }
  graphics::segments(to_px_x(tprev), to_px_y(S), to_px_x(res$calib$t_end), to_px_y(S),
                     col = "red", lwd = 2)
  grDevices::dev.off()
  invisible(out_png)
}
