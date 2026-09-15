# kmextract — reconstruct pseudo individual patient data (IPD) from a Kaplan-Meier
# curve image (PNG/JPEG), constrained by the printed at-risk table.
#
# Method: automated pixel extraction of the KM step function + Guyot-type
# reconstruction (Guyot et al., BMC Med Res Methodol 2012;12:9), with the
# number-at-risk / cumulative-events anchors honored exactly.
#
# Main entry point:
#   km_extract(image_path, x_ticks, anchors, ...)
# returns list(ipd, lifetable, curve, fit) — see README for arguments and examples.
#
# Requires: png (and/or jpeg), survival.
#
# Pixel-geometry constants (tick-label grouping 15px, drop-vs-line threshold
# 8px, reappearance jump cap 150px, tick bands ~25px) are tuned for figures
# roughly 700-2000px wide; validated down to 450x300. Below that, expect loud
# failures (tick-count mismatch), not silent ones.

.km_read_image <- function(path) {
  ext <- tolower(tools::file_ext(path))
  img <- switch(ext,
    png  = png::readPNG(path),
    jpg  = jpeg::readJPEG(path),
    jpeg = jpeg::readJPEG(path),
    stop("Unsupported image format '", ext, "' (use png or jpeg; export PDFs to PNG first)"))
  if (length(dim(img)) == 2) img <- array(rep(img, 3), dim = c(dim(img), 3))
  if (dim(img)[3] == 2) img <- array(rep(img[, , 1], 3), dim = c(dim(img)[1:2], 3))  # gray+alpha
  img[, , 1:3, drop = FALSE]
}

.km_find_axes <- function(g, dark_thresh) {
  dark <- g < dark_thresh
  if (!any(dark)) stop("no dark pixels found - could not detect plot axes ",
                       "(is this a plot image? try raising dark_thresh)")
  h <- nrow(g); w <- ncol(g)
  longest_run <- function(v) {
    if (!any(v)) return(0L)
    r <- rle(v); max(r$lengths[r$values])
  }
  rows <- floor(h * 0.3):h
  x_row <- rows[which.max(vapply(rows, function(r) longest_run(dark[r, ]), 0L))]
  cols <- 1:floor(w * 0.4)
  y_col <- cols[which.max(vapply(cols, function(c) longest_run(dark[, c]), 0L))]
  # extend to the full thickness of the axis lines, on both sides: x_row/y_col
  # are the OUTER edges (bottom / left), x_top/y_right the INNER edges that the
  # plot area starts from — bold-axis themes make these bands several px thick.
  while (x_row + 1L <= h && sum(dark[x_row + 1L, ]) > 0.5 * sum(dark[x_row, ])) x_row <- x_row + 1L
  while (y_col - 1L >= 1L && sum(dark[, y_col - 1L]) > 0.5 * sum(dark[, y_col])) y_col <- y_col - 1L
  x_top <- x_row
  while (x_top - 1L >= 1L && sum(dark[x_top - 1L, ]) > 0.5 * sum(dark[x_row, ])) x_top <- x_top - 1L
  y_right <- y_col
  while (y_right + 1L <= w && sum(dark[, y_right + 1L]) > 0.5 * sum(dark[, y_col])) y_right <- y_right + 1L
  list(x_row = x_row, y_col = y_col, x_top = x_top, y_right = y_right)
}


.km_find_labels <- function(g, axes, dark_thresh) {
  # Fallback calibration for figures WITHOUT axis tick marks (labels only):
  # ggplot centers each tick label under its tick position, so the centroid of
  # each label's text blob calibrates the axis. Digit blobs of one label are
  # merged with a gap scaled to the image size (a fixed small gap would split
  # "20" into "2" and "0").
  dark <- g < dark_thresh
  h <- nrow(g); w <- ncol(g)
  gap_x <- max(30L, as.integer(round(w / 40)))
  gap_y <- max(20L, as.integer(round(h / 40)))
  # x labels: first contiguous text row-block below the x-axis band
  rows <- (axes$x_row + 4L):min(h, axes$x_row + 120L)
  has_txt <- vapply(rows, function(r) any(dark[r, min(axes$y_right + 1L, w):w]), TRUE)
  x_px <- integer(0)
  if (any(has_txt)) {
    idx <- rows[has_txt]
    blocks <- split(idx, cumsum(c(1L, diff(idx) > 4L)))
    b <- blocks[[1]]
    xb <- dark[b, , drop = FALSE]
    x_px <- .km_group_centers(setdiff(which(colSums(xb) > 0), 1:axes$y_right),
                              gap = gap_x)
  }
  # y labels: text in the band left of the y-axis, above the x-axis band
  # (excludes the rotated axis title, which sits further left)
  c1 <- max(1L, axes$y_col - 6L); c0 <- max(1L, axes$y_col - as.integer(round(w / 15)))
  y_px <- integer(0)
  if (c1 > c0) {
    yb <- dark[1:max(1L, axes$x_top - 10L), c0:c1, drop = FALSE]
    y_px <- .km_group_centers(which(rowSums(yb) > 0), gap = gap_y)
  }
  list(x_px = x_px, y_px = y_px)
}


.km_calibrate <- function(g, axes, x_ticks, y_ticks, dark_thresh) {
  # Tick marks first; when their counts do not match the expected values, fall
  # back to label centroids axis by axis.
  ticks <- .km_find_ticks(g, axes, dark_thresh)
  labs <- NULL
  if (length(ticks$x_px) != length(x_ticks)) {
    labs <- .km_find_labels(g, axes, dark_thresh)
    if (length(labs$x_px) == length(x_ticks)) {
      message("x-axis: no usable tick marks (found ", length(ticks$x_px),
              ", expected ", length(x_ticks), ") - calibrated on label centroids")
      ticks$x_px <- labs$x_px
    } else {
      stop("found ", length(ticks$x_px), " x-ticks and ", length(labs$x_px),
           " x-labels, expected ", length(x_ticks))
    }
  }
  if (length(ticks$y_px) != length(y_ticks)) {
    if (is.null(labs)) labs <- .km_find_labels(g, axes, dark_thresh)
    if (length(labs$y_px) == length(y_ticks)) {
      message("y-axis: no usable tick marks (found ", length(ticks$y_px),
              ", expected ", length(y_ticks), ") - calibrated on label centroids")
      ticks$y_px <- labs$y_px
    } else {
      stop("found ", length(ticks$y_px), " y-ticks and ", length(labs$y_px),
           " y-labels, expected ", length(y_ticks))
    }
  }
  ticks
}

.km_group_centers <- function(idx, gap = 15L) {
  if (!length(idx)) return(integer(0))
  b <- cumsum(c(1L, diff(idx) > gap))
  as.integer(round(tapply(idx, b, mean)))
}

.km_find_ticks <- function(g, axes, dark_thresh) {
  dark <- g < dark_thresh
  h <- nrow(g); w <- ncol(g)
  # adaptive: find the first row below the x-axis containing dark pixels
  # (tick marks or the tick-label line), then use a band from there
  r0 <- axes$x_row + 3L
  while (r0 <= min(h, axes$x_row + 60L) && !any(dark[r0, (axes$y_col + 1L):w])) r0 <- r0 + 1L
  xb <- dark[r0:min(h, r0 + 20L), , drop = FALSE]
  x_px <- .km_group_centers(setdiff(which(colSums(xb) > 0), 1:axes$y_col))
  c0 <- axes$y_col - 3L
  while (c0 >= max(1L, axes$y_col - 60L) && !any(dark[1:(axes$x_row - 1L), c0])) c0 <- c0 - 1L
  yb <- dark[, max(1L, c0 - 25L):c0, drop = FALSE]
  ax <- rle(dark[, axes$y_col])                  # vertical span of the y-axis LINE:
  ends <- cumsum(ax$lengths)                     # longest contiguous dark run
  i <- which(ax$values)[which.max(ax$lengths[ax$values])]
  y_lo <- max(1L, ends[i] - ax$lengths[i] + 1L - 8L)
  y_px <- .km_group_centers(intersect(which(rowSums(yb) > 0), y_lo:(axes$x_row - 1L)))
  list(x_px = x_px, y_px = y_px)
}

.km_curve_mask <- function(img, color, dark_thresh, color_tol) {
  if (is.null(color)) {
    g <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
    return(g < dark_thresh)
  }
  d2 <- (img[, , 1] - color[1])^2 + (img[, , 2] - color[2])^2 + (img[, , 3] - color[3])^2
  mx <- pmax(img[, , 1], img[, , 2], img[, , 3])
  mn <- pmin(img[, , 1], img[, , 2], img[, , 3])
  # require some saturation so light-gray gridlines never pass the tolerance
  sqrt(d2) < color_tol & (mx - mn) > 0.07
}

.km_trace <- function(mask, axes, y_top_px, seed) {
  h <- nrow(mask); w <- ncol(mask)
  prev <- if (seed) y_top_px else NA_integer_
  xs <- (axes$y_right + 6L):w
  ys <- rep(NA_integer_, length(xs))
  hh <- rep(0L, length(xs))
  tt <- rep(NA_integer_, length(xs))
  bb <- rep(NA_integer_, length(xs))
  last_seen <- NA_integer_
  for (k in seq_along(xs)) {
    x <- xs[k]
    col <- which(mask[1:(axes$x_top - 2L), x]); col <- col[col > 5L]
    if (!length(col)) { ys[k] <- prev; next }
    last_seen <- x
    if (is.na(prev)) {
      cand <- col[which.min(abs(col - y_top_px))]
      if (abs(cand - y_top_px) <= 12L) { prev <- cand; ys[k] <- prev }
      next
    }
    b <- cumsum(c(1L, diff(col) > 3L))
    runs <- split(col, b)
    joined <- Filter(function(r) min(r) - 4L <= prev && prev <= max(r) + 4L, runs)
    if (length(joined)) {
      r <- joined[[1]]
      hh[k] <- max(r) - min(r) + 1L
      tt[k] <- min(r); bb[k] <- max(r)
      # A tall run SYMMETRIC around the current level is a censor tick ("+")
      # straddling the line, not a drop: keep the trace on the line. A true
      # drop is asymmetric (the run hangs below the upper plateau the trace
      # is on), so it fails the >= 7 px upper-overhang test.
      tick_like <- (max(r) - min(r) > 8L) &&
                   (prev - min(r) >= 7L) && (max(r) - prev >= 7L)
      newy <- if (tick_like) prev
              else if (max(r) - min(r) > 8L) max(prev, max(r) - 2L)   # vertical drop
              else max(prev, as.integer(round((min(r) + max(r)) / 2)))  # flat line: center
      prev <- newy
    } else {
      below <- Filter(function(r) prev < min(r) && min(r) <= prev + 150L, runs)
      if (length(below)) {
        r <- below[[which.min(vapply(below, min, 0L))]]
        hh[k] <- max(r) - min(r) + 1L
        tt[k] <- min(r); bb[k] <- max(r)
        prev <- if (max(r) - min(r) <= 8L) as.integer(round((min(r) + max(r)) / 2)) else max(r) - 2L
      }
    }
    ys[k] <- prev
  }
  keep <- !is.na(ys)
  list(x = xs[keep], y = ys[keep], h = hh[keep],
       top = tt[keep], bot = bb[keep], last_seen = last_seen)
}


.km_censor_ticks <- function(tr, t, S, min_drop) {
  # A censor tick ("+") is a stroke centered on the line: it protrudes ABOVE the
  # line's top edge AND BELOW its bottom edge, roughly symmetrically. A step
  # drop only produces stroke on one side of the traced level (from the upper
  # plateau down to it), so requiring BOTH protrusions rejects drops by
  # construction (the upward overhang alone would still be fooled at the drop
  # column, whose top edge is the upper plateau).
  flat_cols <- tr$h > 0L & tr$h <= 8L          # plain-line columns (cf. drop rule)
  thick <- stats::median(tr$h[flat_cols])
  if (is.na(thick)) thick <- stats::median(tr$h[tr$h > 0L])
  if (is.na(thick)) return(numeric(0))
  m <- max(2L, as.integer(ceiling(thick / 2)))  # required overhang on each side
  half <- as.integer(ceiling(thick / 2))
  cand <- which(tr$h > 0L &
                tr$h <= 6L * max(thick, 3L) &   # a "+" is ~3-4x the line width;
                                                # excludes panel borders/axes
                tr$top <= tr$y - half - m &     # pokes above the line top edge
                tr$bot >= tr$y + half + m)      # AND below the line bottom edge
  out <- numeric(0)
  if (length(cand)) {
    b <- cumsum(c(1L, diff(cand) > 3L))
    for (grp in split(cand, b)) out <- c(out, t[as.integer(round(mean(grp)))])
  }
  out
}

.km_reconstruct <- function(steps, anchors, t_end, cens_hint = numeric(0)) {
  ipd_t <- numeric(0); ipd_e <- integer(0)
  S_prev <- 1
  for (k in seq_len(nrow(anchors) - 1L)) {
    t0 <- anchors$t[k];  n0 <- anchors$n[k];  E0 <- anchors$E[k]
    t1 <- anchors$t[k+1]; n1 <- anchors$n[k+1]; E1 <- anchors$E[k+1]
    D <- E1 - E0; C <- (n0 - n1) - D
    stopifnot("inconsistent anchors (D<0 or C<0)" = D >= 0 && C >= 0)
    seg <- steps[steps$t >= t0 & steps$t < t1, , drop = FALSE]
    if (!nrow(seg) && D > 0L) seg <- data.frame(t = (t0 + t1) / 2, S = NA)
    m <- nrow(seg)
    if (m == 0L && C == 0L) next
    # decision variables: d[j] deaths at step j; cg[g] censors in gap g,
    # where gap g = 0..m lies before step g+1 (gap m = after the last step)
    walk <- function(d, cg) {
      # returns per-step reconstructed S and end-of-interval S, or NULL if infeasible
      if (any(cg < 0L) || any(d < 0L)) return(NULL)
      n <- n0 - cg[1L]
      Sp <- S_prev; Ss <- numeric(m)
      for (j in seq_len(m)) {
        if (d[j] > n) return(NULL)
        if (n > 0L) Sp <- Sp * (1 - d[j] / n)
        Ss[j] <- Sp
        n <- n - d[j] - cg[j + 1L]
      }
      if (n < 0L) return(NULL)
      list(Ss = Ss, S_end = Sp)
    }
    err <- function(d, cg) {
      w <- walk(d, cg)
      if (is.null(w)) return(Inf)
      ok <- !is.na(seg$S)
      if (!any(ok)) return(0)
      sum((w$Ss[ok] - seg$S[ok])^2)
    }
    # ---- initial allocation ----
    gaps_lo <- c(t0, seg$t); gaps_hi <- c(seg$t, t1)
    cg <- integer(m + 1L)
    if (C > 0L) {
      hints <- cens_hint[cens_hint >= t0 & cens_hint < t1]
      for (h in hints[seq_len(min(length(hints), C))]) {
        g <- findInterval(h, c(t0, seg$t))          # 1..m+1
        cg[g] <- cg[g] + 1L
      }
      left <- C - sum(cg)
      if (left > 0L) {                               # spread remainder evenly
        pos <- t0 + (t1 - t0) * seq_len(left) / (left + 1L)
        for (h in pos) { g <- findInterval(h, c(t0, seg$t)); cg[g] <- cg[g] + 1L }
      }
    }
    d <- integer(m)
    if (m > 0L) {
      n <- n0 - cg[1L]; Sp <- S_prev
      for (j in seq_len(m)) {
        d[j] <- if (!is.na(seg$S[j]) && Sp > 0) min(n, max(0L, as.integer(round(n * (1 - seg$S[j] / Sp))))) else 0L
        if (n > 0L && d[j] > 0L) Sp <- Sp * (1 - d[j] / n)
        n <- n - d[j] - cg[j + 1L]
      }
      # exact death total, best-fit adjustments
      diff_d <- D - sum(d); guard <- 0L
      while (diff_d != 0L) {
        best <- NULL; best_e <- Inf
        for (j in seq_len(m)) {
          trial <- d
          if (diff_d > 0L) trial[j] <- trial[j] + 1L
          else if (trial[j] > 0L) trial[j] <- trial[j] - 1L else next
          e <- err(trial, cg)
          if (e < best_e) { best_e <- e; best <- trial }
        }
        stopifnot("no feasible allocation adjustment" = !is.null(best))
        d <- best; diff_d <- diff_d + if (diff_d > 0L) -1L else 1L
        guard <- guard + 1L; stopifnot("allocation loop stuck" = guard < 10000L)
      }
      # ---- greedy refinement (Guyot-style iteration): move single censors
      # between gaps, or single deaths between steps, while the fit improves
      cur <- err(d, cg); guard <- 0L
      repeat {
        improved <- FALSE
        if (C > 0L) for (a in seq_len(m + 1L)) for (b in seq_len(m + 1L)) {
          if (a == b || cg[a] <= 0L) next   # test the CURRENT vector: cg
          trial <- cg                        # mutates inside the loop
          trial[a] <- trial[a] - 1L; trial[b] <- trial[b] + 1L
          e <- err(d, trial)
          if (e < cur - 1e-12) { cg <- trial; cur <- e; improved <- TRUE }
        }
        if (m > 1L) for (a in seq_len(m)) for (b in seq_len(m)) {
          if (a == b || d[a] <= 0L) next
          trial <- d; trial[a] <- trial[a] - 1L; trial[b] <- trial[b] + 1L
          e <- err(trial, cg)
          if (e < cur - 1e-12) { d <- trial; cur <- e; improved <- TRUE }
        }
        guard <- guard + 1L
        if (!improved || guard > 100L) break
      }
    }
    # ---- emit ----
    for (j in seq_len(m)) if (d[j] > 0L) {
      ipd_t <- c(ipd_t, rep(seg$t[j], d[j])); ipd_e <- c(ipd_e, rep(1L, d[j]))
    }
    for (g in seq_len(m + 1L)) if (cg[g] > 0L) {
      lo <- gaps_lo[g]; hi <- gaps_hi[g]
      ipd_t <- c(ipd_t, lo + (hi - lo) * seq_len(cg[g]) / (cg[g] + 1L))
      ipd_e <- c(ipd_e, rep(0L, cg[g]))
    }
    w <- walk(d, cg)
    S_prev <- if (!is.null(w)) w$S_end else S_prev
  }
  nl <- anchors$n[nrow(anchors)]; tl <- anchors$t[nrow(anchors)]
  if (nl > 0L) { ipd_t <- c(ipd_t, rep(max(tl, t_end), nl)); ipd_e <- c(ipd_e, rep(0L, nl)) }
  ipd <- data.frame(time = ipd_t, event = ipd_e)[order(ipd_t), ]
  stopifnot("event total != last anchor" = sum(ipd$event) == anchors$E[nrow(anchors)],
            "N != initial at-risk" = nrow(ipd) == anchors$n[1])
  rownames(ipd) <- NULL
  ipd
}

#' Extract pseudo-IPD from a KM curve image.
#'
#' @param image_path PNG/JPEG of the KM plot (one curve per call; for overlaid
#'   colored curves, call once per curve with its `color`).
#' @param x_ticks numeric values of the x-axis ticks, left to right.
#' @param anchors data.frame(t, n, E): the printed at-risk table — time, number
#'   at risk, CUMULATIVE events. Must start at t = 0 with E = 0.
#' @param t_end time where the curve visibly ends. Patients still at risk at the
#'   last anchor are censored there. If the curve ends with a drop to zero AFTER
#'   the last printed anchor, append a final anchor row with n = 0 and the
#'   corresponding extra event count instead.
#' @param y_ticks y-axis tick values top to bottom (default 1, .75, .5, .25, 0).
#' @param color NULL for a black/dark curve, else c(r, g, b) in 0-1 for a
#'   colored curve (read it off the legend, e.g. with a color picker).
#' @param seed TRUE when overlaid curves share the start at S = 1 and this
#'   curve may be hidden at the beginning (forces the trace to start at S = 1).
#' @return list(ipd = data.frame(time, event),
#'              lifetable = survfit-style table,
#'              curve = extracted step coordinates,
#'              fit = survfit object)
km_extract <- function(image_path, x_ticks, anchors, t_end,
                       y_ticks = c(1, 0.75, 0.5, 0.25, 0),
                       color = NULL, seed = FALSE,
                       dark_thresh = 110/255, color_tol = 70/255) {
  stopifnot(all(c("t", "n", "E") %in% names(anchors)))
  if (is.unsorted(anchors$t)) stop("anchors must be sorted by t")
  if (anchors$E[1] != 0) stop("anchors must start with E = 0 (cumulative events)")
  if (length(x_ticks) < 2) stop("need at least 2 x_ticks for calibration")
  if (t_end < max(anchors$t[anchors$n > 0], 0))
    stop("t_end must be >= the last anchor time with patients still at risk")
  img <- .km_read_image(image_path)
  g <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  axes <- .km_find_axes(g, dark_thresh)
  ticks <- .km_calibrate(g, axes, x_ticks, y_ticks, dark_thresh)
  mask <- .km_curve_mask(img, color, dark_thresh, color_tol)
  y_top <- ticks$y_px[which.max(y_ticks)]
  tr <- .km_trace(mask, axes, y_top, seed)
  if (!length(tr$x) || all(tr$h == 0))
    stop("no curve pixels matched the mask - check `color` (and `seed` for overlaid curves)")
  fx <- lm(v ~ p, data = data.frame(p = ticks$x_px, v = x_ticks))
  fy <- lm(v ~ p, data = data.frame(p = ticks$y_px, v = y_ticks))
  t <- as.numeric(predict(fx, data.frame(p = tr$x)))
  S <- as.numeric(predict(fy, data.frame(p = tr$y)))
  S <- cummin(pmin(pmax(S, 0), 1))
  # step detection
  n0 <- anchors$n[1]; min_drop <- 0.5 / n0
  st_t <- numeric(0); st_S <- numeric(0); S_run <- S[1]
  for (i in seq_along(t)[-1]) if (S_run - S[i] >= min_drop) {
    st_t <- c(st_t, t[i]); st_S <- c(st_S, S[i]); S_run <- S[i]
  }
  cens_hint <- .km_censor_ticks(tr, t, S, min_drop)
  ipd <- .km_reconstruct(data.frame(t = st_t, S = st_S), anchors, t_end,
                         cens_hint = cens_hint)
  fit <- survival::survfit(survival::Surv(time, event) ~ 1, data = ipd)
  lifetable <- data.frame(time = fit$time, n.risk = fit$n.risk,
                          n.event = fit$n.event, n.censor = fit$n.censor)
  list(ipd = ipd, lifetable = lifetable,
       curve = data.frame(t = t, S = S), fit = fit)
}

#' Compare a reconstruction against published summary statistics.
km_check <- function(res, published_median = NA, published_lcl = NA,
                     published_ucl = NA, published_events = NA) {
  s <- summary(res$fit)$table
  data.frame(
    N = s[["records"]], events = s[["events"]],
    median = round(s[["median"]], 2),
    lcl = round(s[["0.95LCL"]], 2), ucl = round(s[["0.95UCL"]], 2),
    pub_median = published_median, pub_lcl = published_lcl,
    pub_ucl = published_ucl, pub_events = published_events)
}


#' Extract pseudo-IPD from a KM image WITHOUT a printed at-risk table.
#'
#' Deaths are read from the step drops; censoring times are read from the
#' censor tick marks ("+") drawn on the curve. Use when the figure has no
#' number-at-risk table (then only N, the starting sample size, is needed).
#' Requires the tick marks to be visible; if they are not, results are
#' unreliable — say so rather than trusting them.
km_extract_free <- function(image_path, x_ticks, N,
                            y_ticks = c(1, 0.75, 0.5, 0.25, 0),
                            color = NULL, seed = FALSE, t_end = NULL,
                            dark_thresh = 110/255, color_tol = 70/255) {
  stopifnot(N >= 1, length(x_ticks) >= 2)
  img <- .km_read_image(image_path)
  g <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  axes <- .km_find_axes(g, dark_thresh)
  ticks <- .km_calibrate(g, axes, x_ticks, y_ticks, dark_thresh)
  mask <- .km_curve_mask(img, color, dark_thresh, color_tol)
  y_top <- ticks$y_px[which.max(y_ticks)]
  tr <- .km_trace(mask, axes, y_top, seed)
  if (!length(tr$x) || all(tr$h == 0))
    stop("no curve pixels matched the mask - check `color` (and `seed` for overlaid curves)")
  fx <- lm(v ~ p, data = data.frame(p = ticks$x_px, v = x_ticks))
  fy <- lm(v ~ p, data = data.frame(p = ticks$y_px, v = y_ticks))
  t <- as.numeric(predict(fx, data.frame(p = tr$x)))
  S <- as.numeric(predict(fy, data.frame(p = tr$y)))
  S <- cummin(pmin(pmax(S, 0), 1))
  if (is.null(t_end)) t_end <- as.numeric(predict(fx, data.frame(p = tr$last_seen)))
  min_drop <- 0.5 / N
  # steps
  st_t <- numeric(0); st_S <- numeric(0); S_run <- S[1]
  for (i in seq_along(t)[-1]) if (S_run - S[i] >= min_drop) {
    st_t <- c(st_t, t[i]); st_S <- c(st_S, S[i]); S_run <- S[i]
  }
  # censor ticks: tall runs on a locally flat stretch of the curve
  if (all(tr$h == 0)) stop("curve trace produced no measurable line - image too degraded?")
  cens_t <- .km_censor_ticks(tr, t, S, min_drop)
  # sequential reconstruction: deaths from drop ratios, censors from ticks
  if (!length(st_t) && !length(cens_t))
    stop("no steps and no censor ticks detected - nothing to reconstruct")
  ev <- rbind(if (length(st_t)) data.frame(t = st_t, S = st_S, type = "d"),
              if (length(cens_t)) data.frame(t = cens_t, S = NA, type = "c"))
  ev <- ev[order(ev$t), ]
  message(length(st_t), " step(s) and ", length(cens_t),
          " censor tick(s) detected - verify the tick count against the figure: ",
          "each missed tick can turn a censoring into a spurious death")
  n <- N; Sp <- 1
  ipd_t <- numeric(0); ipd_e <- integer(0)
  for (i in seq_len(nrow(ev))) {
    if (ev$type[i] == "d") {
      d <- max(1L, as.integer(round(n * (1 - ev$S[i] / Sp))))
      d <- min(d, n)
      ipd_t <- c(ipd_t, rep(ev$t[i], d)); ipd_e <- c(ipd_e, rep(1L, d))
      Sp <- Sp * (1 - d / n); n <- n - d
    } else if (n > 0L) {
      ipd_t <- c(ipd_t, ev$t[i]); ipd_e <- c(ipd_e, 0L); n <- n - 1L
    }
  }
  if (n > 0L) { ipd_t <- c(ipd_t, rep(t_end, n)); ipd_e <- c(ipd_e, rep(0L, n)) }
  ipd <- data.frame(time = ipd_t, event = ipd_e)[order(ipd_t), ]
  rownames(ipd) <- NULL
  fit <- survival::survfit(survival::Surv(time, event) ~ 1, data = ipd)
  list(ipd = ipd,
       lifetable = data.frame(time = fit$time, n.risk = fit$n.risk,
                              n.event = fit$n.event, n.censor = fit$n.censor),
       curve = data.frame(t = t, S = S),
       censor_ticks = cens_t, fit = fit)
}
