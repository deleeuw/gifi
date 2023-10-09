
gifiEngine <-
  function (gifi,
            ndim,
            itmax,
            eps,
            seed,
            verbose) {
    set.seed (seed)
    nobs <- nrow (as.matrix (gifi[[1]][[1]]$data))
    nsets <- length (gifi)
    nvars <- sum (sapply (gifi, length))
    itel <- 1
    if (nvars == 1)
      stop("a gifiAnalysis needs more than one variable")
    x <- matrix (rnorm (nobs * ndim), nobs, ndim)
    x <- gsRC (center (x))$q
    xGifi <- xGifi (gifi, x)
    fold <- 0
    asets <- 0
    for (i in 1:nsets) {
      gifiSet <- gifi[[i]]
      xGifiSet <- xGifi[[i]]
      nvarset <- length (gifiSet)
      ha <- matrix (0, nobs, ndim)
      activeCount <- 0
      for (j in 1:nvarset) {
        if (gifiSet[[j]]$active) {
          activeCount <- activeCount + 1
          ha <- ha + xGifiSet[[j]]$scores
        }
      }
      if (activeCount > 0) {
        asets <- asets + 1
        fold <- fold + sum ((x - ha) ^ 2)
      }
    }
    fold <- fold / (asets * ndim)
    repeat {
      xz <- matrix(0, nobs, ndim)
      fnew <- fmid <- 0
      for (i in 1:nsets) {
        gifiSet <- gifi[[i]]
        xGifiSet <- xGifi[[i]]
        nvarset <- length (gifiSet)
        hh <- matrix (0, nobs, 0)
        activeCount <- 0
        for (j in 1:nvarset) {
          if (gifiSet[[j]]$active) {
            activeCount <- activeCount + 1
            hh <- cbind (hh, xGifiSet[[j]]$transform)
          }
        }
        if (activeCount == 0)
          next
        lf <- lsRC (hh, x)
        aa <- lf$solution
        rs <- lf$residuals
        kappa <- max (eigen (crossprod (aa))$values)
        fmid <- fmid + sum (rs ^ 2)
        target <- hh + tcrossprod (rs, aa) / kappa
        hh <- matrix (0, nobs, 0)
        scopies <- 0
        for (j in 1:nvarset) {
          gifiVar <- gifiSet[[j]]
          jdata <- gifiVar$data
          jbasis <- gifiVar$basis
          jcopies <- gifiVar$copies
          jdegree <- gifiVar$degree
          jties <- gifiVar$ties
          jmissing <- gifiVar$missing
          jordinal <- gifiVar$ordinal
          ja <- aa[scopies + 1:jcopies, , drop = FALSE]
          jtarget <- target[, scopies + 1:jcopies, drop = FALSE]
          hj <-
            gifiTransform (
              data = jdata,
              target = jtarget,
              basis = jbasis,
              copies = jcopies,
              degree = jdegree,
              ordinal = jordinal,
              ties = jties,
              missing = jmissing
            )
          hj <- gsRC(normalize (center (hj)))$q
          sc <- hj %*% ja
          xGifi[[i]][[j]]$transform <- hj
          xGifi[[i]][[j]]$weights <- ja
          xGifi[[i]][[j]]$scores <- sc
          xGifi[[i]][[j]]$quantifications <-
            lsRC(jbasis, sc)$solution
          activeCount <- 0
          if (gifiSet[[j]]$active) {
            activeCount <- activeCount + 1
            hh <- cbind (hh, hj)
          }
          scopies <- scopies + jcopies
        }
        if (activeCount > 0) {
          ha <- hh %*% aa
          xz <- xz + ha
          fnew <- fnew + sum ((x - ha) ^ 2)
        }
      }
      fmid <- fmid / (asets * ndim)
      fnew <- fnew / (asets * ndim)
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
       if (((itel == itmax) || ((fold - fnew) < eps)) && (itel > 1))
        break
      itel <- itel + 1
      fold <- fnew
      x <- gsRC (center (xz))$q
    }
    return (list (
      f = fnew,
      ntel = itel,
      x = x,
      xGifi = xGifi
    ))
  }

gifiTransform <-
  function (data,
            target,
            basis,
            copies,
            degree,
            ordinal,
            ties,
            missing) {
    nobs <- nrow (as.matrix (data))
    h <- matrix (0, nobs, copies)
    if (degree == -1) {
      if (ordinal) {
        h[, 1] <-
          coneRegression (
            data = data,
            target = target[, 1],
            type = "c",
            ties = ties,
            missing = missing
          )
      }
      else {
        h[, 1] <-
          coneRegression (
            data = data,
            target = target[, 1],
            basis = basis,
            type = "s",
            missing = missing
          )
      }
    }
    if (degree >= 0) {
      if (ordinal) {
        h[, 1] <-
          coneRegression (
            data = data,
            target = target[, 1],
            basis = basis,
            type = "i",
            ties = ties,
            missing = missing
          )
      }
      else {
        h[, 1] <-
          coneRegression (
            data = data,
            target = target[, 1],
            basis = basis,
            type = "s",
            ties = ties,
            missing = missing
          )
      }
    }
    if (copies > 1) {
      for (l in 2:copies)
        h[, l] <-
          coneRegression (
            data = data,
            target = target[, l],
            basis = basis,
            type = "s",
            ties = ties,
            missing = missing
          )
    }
    return (h)
  }
