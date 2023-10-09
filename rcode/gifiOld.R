source("splinebasis.R")
source("matrix.R")
library (nnls)
dyn.load("matrix.so")
dyn.load("splinebasis.so")

gifiEngine <-
  function (data,
            knots,
            degrees,
            ordinal,
            sets = 1:ncol(data),
            copies = rep (1, ncol(data)),
            ndim = 2,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = TRUE) {
    nvars <- ncol (data)
    nobs <- nrow (data)
    nsets <- max (sets)
    ncops <- sum (copies)
    lmax <- ndim * nsets
    indices <- gifiIndices (nvars, nsets, ordinal, sets, copies)
    setup <-
      gifiSetup (data, knots, degrees, ordinal, sets , copies, indices$expand)
    initials <-
      gifiInitials (
        setup$q,
        setup$r,
        nobs,
        nsets,
        ncops,
        ndim,
        indices$expand,
        indices$nvarsets,
        indices$ordicops,
        seed
      )
    fold <- 0
    itel <- 1
    setcopies <- indices$setcopies
    h <- initials$h
    a <- initials$a
    x <- initials$x
    for (s in 1:nsets) {
      hh <- h[1:nobs, which (setcopies == s), drop = FALSE]
      aa <- a[which (setcopies == s), 1:ndim, drop = FALSE]
      fold <- fold + sum ((x - hh %*% aa) ^ 2)
    }
    fold <- fold / lmax
    repeat {
      xz <- matrix(0, nobs, ndim)
      fnew <- fmid <- 0
      for (s in 1:nsets) {
        id <- which (setcopies == s)
        hh <- h[1:nobs, id, drop = FALSE]
        lf <- lsRC (hh, x)
        aa <- lf$solution
        rs <- lf$residuals
        kappa <- max (eigen (crossprod (aa))$values)
        fmid <- fmid + sum (rs ^ 2)
        target <- hh + tcrossprod (rs, aa) / kappa
        k <- 1
        for (l in id) {
          ql <- setup$q[[l]]
          vl <- setup$v[[l]]
          if (indices$ordicops[l]) {
            hz <- drop (crossprod (ql, target[, k]))
            ns <- nnls (vl, hz)
            rz <- coefficients (ns)
            hk <- drop (ql %*% (hz - drop (vl %*% rz)))
          }
          else {
            hk <- ql %*% crossprod (ql, target[, k])
          }
          hh[, k] <- hk / sqrt (sum (hk ^ 2))
          k <- k + 1
        }
        ha <- hh %*% aa
        xz <- xz + ha
        fnew <- fnew + sum ((x - ha) ^ 2)
        h[1:nobs, id] <- hh
        a[id, 1:ndim] <- aa
      }
      fmid <- fmid / lmax
      fnew <- fnew / lmax
      if (verbose)
        cat(
          "Iteration: ",
          formatC (itel, width = 3, format = "d"),
          "fold: ",
          formatC (
            fold,
            digits = 8,
            width = 12,
            format = "f"
          ),
          "fmid: ",
          formatC (
            fmid,
            digits = 8,
            width = 12,
            format = "f"
          ),
          "fnew: ",
          formatC (
            fnew,
            digits = 8,
            width = 12,
            format = "f"
          ),
          "\n"
        )
      if ((itel == itmax) || ((fold - fnew) < eps))
        break
      itel <- itel + 1
      fold <- fnew
      x <- gsRC (center (xz))$q
    }
    return (list (
      f = fnew,
      ntel = itel,
      x = x,
      a = a,
      h = h,
      q = setup$q,
      r = setup$r
    ))
  }

gifiIndices <- function (nvars, nsets, ordinal, sets, copies) {
  ordicops <- logical (0)
  expand <- numeric (0)
  for (j in 1:nvars) {
    ordicops <- c (ordicops, ordinal [j], rep (FALSE, copies[j] - 1))
    expand <- c (expand, rep (j, copies [j]))
  }
  setcopies <- sets[expand]
  nvarsets <- numeric(0)
  for (s in 1:nsets) {
    nvarsets <- c(nvarsets, sum (copies[which (sets == s)]))
  }
  return (list (
    ordicops = ordicops,
    setcopies = setcopies,
    expand = expand,
    nvarsets = nvarsets
  ))
}

gifiSetup <-
  function (data,
            knots,
            degrees,
            ordinal,
            sets,
            copies,
            expand) {
    nvars <- ncol (data)
    nobs <- nrow (data)
    ncops <- sum (copies)
    g <- b <- q <- r <- v <- list()
    for (l in 1:ncops) {
      j <- expand[l]
      mm <- is.na(data[, j])
      nm <- !mm
      dm <- data[nm, j]
      if (degrees[j] < 0) {
        gn <- ifelse (outer (dm, unique (dm), "=="), 1, 0)
      }
      else {
        gn <- bsplineBasis (dm, degrees[j], knots[[j]])
      }
      nr <- ncol (gn)
      ns <- sum (mm)
      gg <- matrix (0, nobs, nr + ns)
      gg[!mm, 1:nr] <- gn
      if (ns > 0) {
        gg[mm, nr + (1:ns)] <- diag (ns)
      }
      gg <- center (gg)[,-1, drop = FALSE]
      g <- c (g, list (gg))
      bb <- gsRC (g[[l]])
      q <- c(q, list (bb$q))
      r <- c (r, list (bb$r))
      vv <- makeV (data[, j], q[[l]])
      v <- c (v, list (vv))
    }
    return (list (q = q, r = r, v = v))
  }

gifiInitials <-
  function (q,
            r,
            nobs,
            nsets,
            ncops,
            ndim,
            expand,
            nvarsets,
            ordicops,
            seed) {
    set.seed (seed)
    x <- matrix (rnorm (nobs * ndim), nobs, ndim)
    x <- gsRC (center (x))$q
    a <- matrix (rnorm (ncops * ndim), ncops, ndim)
    h <- matrix (0, nobs, ncops)
    for (l in 1:ncops) {
      ql <- q[[l]]
      rl <- r[[l]]
      if (ordicops[l])
        cf <- 1:ncol (ql)
      else
        cf <- rnorm (ncol (ql))
      h[, l] <- drop (ql %*% rl %*% cf)
    }
    h <- center (h)
    h <- apply (h, 2, function (z)
      z / sqrt (sum (z ^ 2)))
    return (list (x = x, a = a, h = h))
  }

makeV <- function (x, g) {
  r <- order (x)[1:length(which(!is.na(x)))]
  n <- length (r)
  m <- ncol (g)
  v <- numeric (0)
  k <- 0
  for (i in 1:(n - 1)) {
    ri <- r[i]
    rj <- r[i + 1]
    if (is.na(x[ri]) || is.na(x[rj]))
      next
    if (x[rj] > x[ri]) {
      v <- c(v, g[ri,] - g[rj,])
      k <- k + 1
    }
  }
  return (t(matrix (v, k, m, byrow = TRUE)))
}

