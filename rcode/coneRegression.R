
amalgm <- function (x, w = rep (1, length (x))) {
  n <- length (x)
  a <- rep (0, n)
  b <- rep (0, n)
  y <- rep (0, n)
  lf <-
    .Fortran (
      "AMALGM",
      n = as.integer (n),
      x = as.double (x),
      w = as.double (w),
      a = as.double (a),
      b = as.double (b),
      y = as.double (y),
      tol = as.double(1e-15),
      ifault = as.integer(0)
    )
  return (lf$y)
}

isotone <-
  function (x,
            y,
            w = rep (1, length (x)),
            ties = "s") {
    there <- which (!is.na (x))
    notthere <- which (is.na (x))
    xthere <- x[there]
    f <- sort(unique(xthere))
    g <- lapply(f, function (z)
      which(x == z))
    n <- length (x)
    k <- length (f)
    if (ties == "s") {
      w <- sapply (g, length)
      h <- lapply (g, function (z)
        y[z])
      m <- sapply (h, sum) / w
      r <- amalgm (m, w)
      s <- rep (0, n)
      for (i in 1:k)
        s[g[[i]]] <- r[i]
    s[notthere] <- y[notthere]
    }
    if (ties == "p") {
      h <- lapply (g, function (z)
        y[z])
      m <- rep (0, n)
      s <- rep (0, n)
      for (i in 1:k) {
        ii <- order (h[[i]])
        g[[i]] <- g[[i]][ii]
        h[[i]] <- h[[i]][ii]
      }
      m <- unlist (h)
      r <- amalgm (m, w)
      s[there] <- r[order (unlist (g))]
      s[notthere] <- y[notthere]
    }
    if (ties == "t") {
      w <- sapply (g, length)
      h <- lapply (g, function (x)
        y[x])
      m <- sapply (h, sum) / w
      r <- amalgm (m, w)
      s <- rep (0, n)
      for (i in 1:k)
        s[g[[i]]] <- y[g[[i]]] + (r[i] - m[i])
      s[notthere] <- y[notthere]
    }
    return (s)
  }

coneRegression <-
  function (data,
            target,
            basis = matrix (data, length(data), 1),
            type = "i",
            ties = "s",
            missing = "m",
            itmax = 1000,
            eps = 1e-6) {
    itel <- 1
    there <- which (!is.na (data))
    notthere <- which (is.na (data))
    nmis <- length (notthere)
    solution <- rep(0, length (data))
    wdata <- data[there]
    wtarget <- target[there]
    wbasis <- basis [there, ]
    if (type == "s")  {
      solution  <- drop (basis %*% lm.fit (basis, target)$coefficients)
    }
     if ((type == "c") && (missing != "a")) {
      solution[there] <- isotone (x = wdata, y = wtarget, ties = ties)
      if (nmis > 0) {
        if (missing == "m")
          solution[notthere] <- target[notthere]
        if (missing == "s")
          solution[notthere] <- mean (target[notthere])
      }
     }
    if ((type == "i")  || ((type == "c") && (missing == "a"))) {
      solution <-
        dykstra (
          target = target,
          basis = basis,
          data = data,
          ties = ties,
          itmax = itmax,
          eps = eps
        )
    }
    return (solution)
  }

dykstra <- function (target, basis, data, ties, itmax, eps) {
  x0 <- target
  itel <- 1
  a <- b <- rep (0, length (target))
  fold <- Inf
  repeat {
    x1 <- drop (basis %*% lm.fit (basis, x0 - a)$coefficients)
    a <- a + x1 - x0
    x2 <- isotone (data, x1 - b, ties = ties)
    b <- b + x2 - x1
    fnew <- sum ((target - (x1 + x2) / 2) ^ 2)
    xdif <- max (abs (x1 - x2))
    if ((itel == itmax) || (xdif < eps))
      break
    itel <- itel + 1
    x0 <- x2
    fold <- fnew
  }
  return ((x1 + x2) / 2)
}