
aspectEngine <-
  function (gifi,
            afunc,
            eps = 1e-6,
            itmax = 100,
            verbose = 1,
            monotone = FALSE,
            ...) {
    nsets <- length (gifi)
    for (i in 1:nsets) {
      gifiSet <- gifi[[i]]
      nvars <- length (gifiSet)
      for (j in 1:nvars) {
        gifiVar <- gifiSet[[j]]
        q <- gifiVar$qr$q
      }
    }
    itel <- 1
    tdata <- matrix (0, n, m)
    for (j in 1:m) {
      tdata[, j] <- bd$x[[j]]
    }
    tdata <- apply (tdata, 2, function (z)
      z - mean (z))
    tdata <- apply (tdata, 2, function (z)
      z / sqrt (sum (z ^ 2)))
    corr <- crossprod (tdata)
    af <- afunc(corr, ...)
    fold <- af$f
    g <- af$g
    repeat {
      for (j in 1:m) {
        target <- drop (tdata[, -j] %*% g[-j, j])
        k <- bd$b[[j]]
        v <- bd$v[[j]]
        u <- drop (crossprod(k, target))
        s0 <- sum(target * tdata[, j])
        if (ordinal[j]) {
          ns <- nnls(v, u)
          rr <- residuals(ns)
          tt <- drop(k %*% rr)
        } else
          tt <- drop (k %*% u)
        tt <- tt - mean (tt)
        sq <- sum(tt ^ 2)
        if (sq > 1e-15) {
          tt <- tt / sqrt (sq)
          tdata[, j] <- tt
        }
        s1 <- sum(target * tdata[, j])
        if (verbose > 1)
          cat (
            "**** Variable",
            formatC(j, width = 3, format = "d"),
            "Before",
            formatC(
              s0,
              digits = 8,
              width = 12,
              format = "f"
            ),
            "After",
            formatC(
              s1,
              digits = 8,
              width = 12,
              format = "f"
            ),
            "\n"
          )
        if (!monotone) {
          corr <- cor (tdata)
          af <- afunc (corr, ...)
          fnew <- af$f
          g <- af$g
        }
      }
      if (monotone) {
        corr <- cor (tdata)
        af <- afunc (corr, ...)
        fnew <- af$f
        g <- af$g
      }
      if (verbose > 0)
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
          "fnew: ",
          formatC (
            fnew,
            digits = 8,
            width = 12,
            format = "f"
          ),
          "\n"
        )
      if ((itel == itmax) || ((fnew - fold) < eps))
        break
      itel <- itel + 1
      fold <- fnew
    }
    return (list (
      tdata = tdata,
      f = fnew,
      r = corr,
      g = g
    ))
  }
