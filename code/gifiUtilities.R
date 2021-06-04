
makeNumeric <- function (x) {
  do <- function (y) {
    u <- unique (y)
    return (drop (ifelse (outer (y, u, "=="), 1, 0) %*% (1:length (u))))
  }
  if (is.vector (x)) {
    return (do (x))
  }
  else {
    return (apply (x, 2, do))
  }
}

center <- function (x) {
  do <- function (z) {
    z - mean (z)
  }
  if (is.matrix (x))
    return (apply (x, 2, do))
  else
    return (do (x))
}

normalize <- function (x) {
  do <- function (z) {
    z / sqrt (sum (z ^ 2))
  }
  if (is.matrix (x))
    return (apply (x, 2, do))
  else
    return (do (x))
}

makeMissing <- function (data, basis, missing) {
  there <- which (!is.na (data))
  notthere <- which (is.na (data))
  nmis <- length (notthere)
  nobs <- length (data)
  ndim <- ncol (basis)
  if (missing == "m") {
    abasis <- matrix (0, nobs, ndim + nmis)
    abasis [there, 1:ndim] <- basis
    abasis [notthere, ndim + 1:nmis] <- diag(nmis)
    basis <- abasis
  }
  if (missing == "a") {
    abasis <- matrix (0, nobs, ndim)
    abasis [there,] <- basis
    abasis [notthere,] <- 1 / ndim
    basis <- abasis
  }
  if (missing == "s") {
    abasis <- matrix (0, nobs, ndim + 1)
    abasis [there, 1:ndim] <- basis
    abasis [notthere, ndim + 1] <- 1
    basis <- abasis
  }
  return (basis)
}

makeIndicator <- function (x) {
  return (as.matrix(ifelse(outer(
    x, sort(unique(x)), "=="
  ), 1, 0)))
}

reshape <- function (x, n) {
  if (length (x) == 1)
    return (rep (x, n))
  else
    return (x)
}

aline <- function (a) {
  abline (0, a[2] / a[1])
}

aperp <- function (a, x) {
  abline (x * (sum (a ^ 2) / a[2]),-a[1] / a[2])
}

aproj <- function (a, h, x) {
  mu <- (h - sum (a * x)) / (sum (a ^ 2))
  return (x + mu * a)
}

stepPlotter <- function (x, y, knots, xlab) {
  y <- as.matrix (y)
  plot (x,
        y[, 1],
        type = "n",
        xlab = xlab,
        ylab = "Transform")
  nknots <- length (knots)
  knots <- c(min(x) - 1, knots, max(x) + 1)
  for (i in 1:(nknots + 1)) {
    ind <- which ((x >= knots [i]) & (x < knots[i + 1]))
    lev <- median (y [ind, 1])
    lines (rbind (c(knots[i], lev), c (knots[i + 1], lev)), col = "RED", lwd = 3)
    if (ncol (y) == 2) {
      lev <- median (y [ind, 2])
      lines (rbind (c(knots[i], lev), c (knots[i + 1], lev)), col = "BLUE", lwd = 3)
    }
  }
}

starPlotter <- function (x, y, main = "") {
  plot(
    x,
    xlab = "dimension 1",
    ylab = "dimension 2",
    col = "RED",
    cex = .5,
    main = main
  )
  points(y, col = "BLUE", cex = .5)
  for (i in 1:nrow(x))
    lines (rbind (x[i, ], y[i, ]))
}

regressionPlotter <-
  function (table,
            x,
            y,
            xname = "Columns",
            yname = "Rows",
            main = "",
            lines = TRUE,
            cex = 1.0,
            ticks = "n") {
    if (ticks != "n") {
      ticks <- NULL
    }
    nr <- nrow (table)
    nc <- ncol (table)
    sr <- rowSums (table)
    sc <- colSums (table)
    rc <- sum (table)
    x <- x - sum (sr * x) / rc
    y <- y - sum (sc * y) / rc
    x <- x / sqrt (sum (sr * (x ^ 2)) / rc)
    y <- y / sqrt (sum (sc * (y ^ 2)) / rc)
    ar <- drop ((table %*% y) / sr)
    ac <- drop ((x %*% table) / sc)
    plot (
      0,
      xlim = c (min(y), max(y)),
      ylim = c (min(x), max(x)),
      xlab = xname,
      ylab = yname,
      main = main,
      xaxt = ticks,
      yaxt = ticks,
      type = "n"
    )
    if (lines) {
      for (i in 1:nr)
        abline (h = x[i])
      for (j in 1:nc)
        abline (v = y[j])
    }
    for (i in 1:nr) {
      for (j in 1:nc) {
        text(y[j],
             x[nr - i + 1],
             as.character(table[i, j]),
             cex = cex,
             col = "RED")
      }
    }
    lines (y, ac, col = "BLUE")
    lines (ar, x, col = "BLUE")
  }

corList <- function (x) {
  m <- length (x)
  n <- nrow (x[[1]])
  h <- matrix (0, n, 0)
  for (i in 1:m) {
    h <- cbind (h, x[[i]])
  }
  return (cor (h))
}

preCorals <- function (x) {
  n <- sum (x)
  r <- nrow (x)
  s <- ncol (x)
  v <- numeric (0)
  for (i in 1:r)
    for (j in 1:s)
      v <- c(v, rep(c(i, j), x[i, j]))
  return (matrix (v, n, 2, byrow = TRUE))
}

postCorals <- function (ff, x) {
  y <- matrix(0, max(ff), ncol (x))
  for (i in 1:nrow (x))
    y[ff[i],] <- x[i,]
  return (y)
}

preCoranals <- function (x, y) {
  n <- sum (x)
  m <- ncol (y)
  r <- nrow (x)
  s <- ncol (x)
  v <- numeric (0)
  for (i in 1:r)
    for (j in 1:s)
      v <- c(v, rep(c(y[i,], j), x[i, j]))
  return (matrix (v, n, m + 1, byrow = TRUE))
}

mprint <- function (x, d = 2, w = 5) {
  print (noquote (formatC (
    x, di = d, wi = w, fo = "f"
  )))
}

burtTable <- function (gifi) {
  nsets <- length (gifi)
  nobs <- length(gifi[[1]][[1]]$data)
  hh <- matrix (0, nobs, 0)
  hl <- list ()
  for (i in 1:nsets) {
    gifiSet <- gifi[[i]]
    nvars <- length (gifiSet)
    hi <- matrix(0, nobs, 0)
    for (j in 1:nvars) {
      gifiVariable <- gifiSet[[j]]
      hi <- cbind (hi, gifiVariable$basis)
    }
    hl <- c (hl, list (crossprod (hi)))
    hh <- cbind (hh, hi)
  }
  return (list (cc = crossprod (hh), dd = directSum (hl)))
}

interactiveCoding <- function (data) {
  cmin <- apply (data, 2, min)
  cmax <- apply (data, 2, max)
  if (!all(cmin == 1))
    stop ("data must be start at 1")
  nobs <- nrow(data)
  h <- numeric(0)
  for (i in 1:nobs)
    h <- c(h, decode (data[i, ], cmax))
  return (h)
}

makeColumnProduct <- function (x) {
  makeTwoColumnProduct <- function (a, b) {
    n <- nrow (a)
    ma <- ncol (a)
    mb <- ncol (b)
    ab <- matrix (0, n, ma * mb)
    k <- 1
    for (i in 1:ma) {
      for (j in 1:mb) {
        ab[, k] <- a[, i] * b[, j]
        k <- k + 1
      }
    }
    return (ab)
  }
  if (!is.list(x)) {
    x <- list (x)
  }
  m <- length (x)
  z <- matrix (1, nrow(x[[1]]), 1)
  for (k in 1:m)
    z <- makeTwoColumnProduct (z, x[[k]])
  return (z)
}


profileFrequencies <- function (data) {
  h <- interactiveCoding (data)
  cmax <- apply (data, 2, max)
  u <- unique (h)
  m <- length (u)
  g <- ifelse (outer (h, u, "=="), 1, 0)
  n <- colSums (g)
  h <- matrix (0, m, ncol (data))
  for (j in 1:m)
    h[j, ] <- encode (u[j], cmax)
  return (list (h = h, n = n))
}


directSum <- function (x) {
  m <- length (x)
  nr <- sum (sapply (x, nrow))
  nc <- sum (sapply (x, ncol))
  z <- matrix (0, nr, nc)
  kr <- 0
  kc <- 0
  for (i in 1:m) {
    ir <- nrow (x[[i]])
    ic <- ncol (x[[i]])
    z[kr + (1:ir), kc + (1:ic)] <- x[[i]]
    kr <- kr + ir
    kc <- kc + ic
  }
  return (z)
}
