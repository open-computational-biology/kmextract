# kmextract — reconstruct pseudo individual patient data (IPD) from a Kaplan-Meier
# curve image (PNG/JPEG), constrained by the printed at-risk table.
#
# Method: automated pixel extraction of the KM step function + Guyot-type
# reconstruction (Guyot et al., BMC Med Res Methodol 2012;12:9), with the
# number-at-risk / cumulative-events anchors honored exactly.
#
# Main entry point:
#   km_extract(image_path, x_ticks, anchors, ...)
# returns list(ipd, lifetable, curve) — see README for arguments and examples.
#
# Requires: png (and/or jpeg), survival.

.km_read_image <- function(path) {
  ext <- tolower(tools::file_ext(path))
  img <- switch(ext,
    png  = png::readPNG(path),
    jpg  = jpeg::readJPEG(path),
    jpeg = jpeg::readJPEG(path),
    stop("Unsupported image format '", ext, "' (use png or jpeg; export PDFs to PNG first)"))
  if (length(dim(img)) == 2) img <- array(rep(img, 3), dim = c(dim(img), 3))
  img[, , 1:3, drop = FALSE]
}

.km_find_axes <- function(g, dark_thresh) {
  dark <- g < dark_thresh
  h <- nrow(g); w <- ncol(g)
  longest_run <- function(v) {
    if (!any(v)) return(0L)
    r <- rle(v); max(r$lengths[r$values])
  }
  rows <- floor(h * 0.3):h
  x_row <- rows[which.max(vapply(rows, function(r) longest_run(dark[r, ]), 0L))]
  cols <- 1:floor(w * 0.4)
  y_col <- cols[which.max(vapply(cols, function(c) longest_run(dark[, c]), 0L))]
  # extend to the full thickness of the axis lines
  while (x_row + 1L <= h && sum(dark[x_row + 1L, ]) > 0.5 * sum(dark[x_row, ])) x_row <- x_row + 1L
  while (y_col - 1L >= 1L && sum(dark[, y_col - 1L]) > 0.5 * sum(dark[, y_col])) y_col <- y_col - 1L
  list(x_row = x_row, y_col = y_col)
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
  xs <- (axes$y_col + 6L):w
  ys <- rep(NA_integer_, length(xs))
  hh <- rep(0L, length(xs))
  last_seen <- NA_integer_
  for (k in seq_along(xs)) {
    x <- xs[k]
    col <- which(mask[1:(axes$x_row - 2L), x]); col <- col[col > 5L]
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
      newy <- if (max(r) - min(r) > 8L) max(prev, max(r) - 2L)   # vertical drop
              else max(prev, as.integer(round((min(r) + max(r)) / 2)))  # flat line: center
      prev <- newy
    } else {
      below <- Filter(function(r) prev < min(r) && min(r) <= prev + 150L, runs)
      if (length(below)) {
        r <- below[[which.min(vapply(below, min, 0L))]]
        prev <- if (max(r) - min(r) <= 8L) as.integer(round((min(r) + max(r)) / 2)) else max(r) - 2L
      }
    }
    ys[k] <- prev
  }
  keep <- !is.na(ys)
  list(x = xs[keep], y = ys[keep], h = hh[keep], last_seen = last_seen)
}

.km_reconstruct <- function(steps, anchors, t_end) {
  ipd_t <- numeric(0); ipd_e <- integer(0)
  S_prev <- 1
  for (k in seq_len(nrow(anchors) - 1L)) {
    t0 <- anchors$t[k];  n0 <- anchors$n[k];  E0 <- anchors$E[k]
    t1 <- anchors$t[k+1]; n1 <- anchors$n[k+1]; E1 <- anchors$E[k+1]
    D <- E1 - E0; C <- (n0 - n1) - D
    stopifnot("inconsistent anchors (D<0 or C<0)" = D >= 0 && C >= 0)
    seg <- steps[steps$t >= t0 & steps$t < t1, , drop = FALSE]
    n <- n0; Sp <- S_prev
    d <- integer(nrow(seg))
    if (nrow(seg)) for (i in seq_len(nrow(seg))) {
      d[i] <- if (Sp > 0) max(0L, as.integer(round(n * (1 - seg$S[i] / Sp)))) else 0L
      n <- n - d[i]; Sp <- seg$S[i]
    }
    diff_d <- D - sum(d)
    if (nrow(seg) && diff_d != 0L) {
      ord <- order(d, decreasing = TRUE); i <- 0L; guard <- 0L
      while (diff_d != 0L) {
        j <- ord[(i %% length(ord)) + 1L]
        if (diff_d > 0L) { d[j] <- d[j] + 1L; diff_d <- diff_d - 1L }
        else if (d[j] > 0L) { d[j] <- d[j] - 1L; diff_d <- diff_d + 1L }
        i <- i + 1L; guard <- guard + 1L
        stopifnot("allocation loop stuck" = guard < 10000L)
      }
    } else if (!nrow(seg) && D > 0L) {
      seg <- data.frame(t = (t0 + t1) / 2, S = NA); d <- D
    }
    if (nrow(seg)) for (i in seq_len(nrow(seg))) if (d[i] > 0L) {
      ipd_t <- c(ipd_t, rep(seg$t[i], d[i])); ipd_e <- c(ipd_e, rep(1L, d[i]))
    }
    if (C > 0L) {
      ipd_t <- c(ipd_t, t0 + (t1 - t0) * seq_len(C) / (C + 1))
      ipd_e <- c(ipd_e, rep(0L, C))
    }
    n <- n0; Sp <- S_prev
    if (nrow(seg)) for (i in seq_len(nrow(seg))) { if (n > 0) Sp <- Sp * (1 - d[i] / n); n <- n - d[i] }
    S_prev <- Sp
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
  img <- .km_read_image(image_path)
  g <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  axes <- .km_find_axes(g, dark_thresh)
  ticks <- .km_find_ticks(g, axes, dark_thresh)
  if (length(ticks$x_px) != length(x_ticks))
    stop("found ", length(ticks$x_px), " x-ticks, expected ", length(x_ticks),
         " (at px ", paste(ticks$x_px, collapse = ","), ")")
  if (length(ticks$y_px) != length(y_ticks))
    stop("found ", length(ticks$y_px), " y-ticks, expected ", length(y_ticks))
  mask <- .km_curve_mask(img, color, dark_thresh, color_tol)
  y_top <- ticks$y_px[which.max(y_ticks)]
  tr <- .km_trace(mask, axes, y_top, seed)
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
  ipd <- .km_reconstruct(data.frame(t = st_t, S = st_S), anchors, t_end)
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
  img <- .km_read_image(image_path)
  g <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
  axes <- .km_find_axes(g, dark_thresh)
  ticks <- .km_find_ticks(g, axes, dark_thresh)
  if (length(ticks$x_px) != length(x_ticks))
    stop("found ", length(ticks$x_px), " x-ticks, expected ", length(x_ticks))
  if (length(ticks$y_px) != length(y_ticks))
    stop("found ", length(ticks$y_px), " y-ticks, expected ", length(y_ticks))
  mask <- .km_curve_mask(img, color, dark_thresh, color_tol)
  y_top <- ticks$y_px[which.max(y_ticks)]
  tr <- .km_trace(mask, axes, y_top, seed)
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
  thick <- stats::median(tr$h[tr$h > 0])
  flat_at <- function(k) {
    lo <- max(1L, k - 8L); hi <- min(length(S), k + 8L)
    (S[lo] - S[hi]) < min_drop / 2
  }
  cand <- which(tr$h >= thick + 5L & vapply(seq_along(tr$h), flat_at, TRUE))
  cens_t <- numeric(0)
  if (length(cand)) {
    b <- cumsum(c(1L, diff(cand) > 3L))
    for (grp in split(cand, b)) cens_t <- c(cens_t, t[as.integer(round(mean(grp)))])
  }
  # sequential reconstruction: deaths from drop ratios, censors from ticks
  ev <- rbind(if (length(st_t)) data.frame(t = st_t, S = st_S, type = "d"),
              if (length(cens_t)) data.frame(t = cens_t, S = NA, type = "c"))
  ev <- ev[order(ev$t), ]
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
