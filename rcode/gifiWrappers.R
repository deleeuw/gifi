homals <-
  function (data,
            knots = knotsD (data),
            degrees = -1,
            ordinal = FALSE,
            ndim = 2,
            ties = "s",
            missing = "m",
            names = colnames (data, do.NULL = FALSE),
            active = TRUE,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE)  {
    nvars <- ncol (data)
    g <- makeGifi (
      data = data,
      knots = knots,
      degrees = reshape (degrees, nvars),
      ordinal = reshape (ordinal, nvars),
      ties = reshape (ties, nvars),
      copies = rep (ndim, ncol (data)),
      missing = reshape (missing, nvars),
      active = reshape (active, nvars),
      names = names,
      sets = 1:nvars
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    a <- v <- z <- d <- y <- o <- as.list (1:ncol(data))
    dsum <- matrix (0, ndim, ndim)
    nact <- 0
    for (j in 1:nvars) {
      jgifi <- h$xGifi[[j]][[1]]
      v[[j]] <- jgifi$transform
      a[[j]] <- jgifi$weights
      y[[j]] <- jgifi$scores
      z[[j]] <- jgifi$quantifications
      cy <- crossprod (y[[j]])
      if (g[[j]][[1]]$active) {
        dsum <- dsum + cy
        nact <- nact + 1
      }
      d[[j]] <- cy
      o[[j]] <- crossprod (h$x, v[[j]])
    }
    return (structure (
      list (
        transform = v,
        rhat = corList (v),
        objectscores = h$x,
        scores = y,
        quantifications = z,
        dmeasures = d,
        lambda = dsum / nact,
        weights = a,
        loadings = o,
        ntel = h$ntel,
        f = h$f
      ),
      class = "homals"
    ))
  }

corals <-
  function (data,
            ftype = TRUE,
            xknots = NULL,
            yknots = NULL,
            xdegree = -1,
            ydegree = -1,
            xordinal = FALSE,
            yordinal = FALSE,
            xties = "s",
            yties = "s",
            xmissing = "m",
            ymissing = "m",
            xname = "X",
            yname = "Y",
            ndim = 2,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE) {
    if (ftype) {
      xy <- preCorals (as.matrix(data))
      x <- xy[, 1, drop = FALSE]
      y <- xy[, 2, drop = FALSE]
    } else {
      x <- data[, 1, drop = FALSE]
      y <- data[, 2, drop = FALSE]
    }
    if (is.null(xknots))
      xknots <- knotsD(x)
    if (is.null(yknots))
      yknots <- knotsD(y)
    g <- makeGifi (
      data = cbind (x, y),
      knots = c(xknots, yknots),
      degrees = c(xdegree, ydegree),
      ordinal = c(xordinal, yordinal),
      ties = c(xties, yties),
      copies = c(ndim, ndim),
      missing = c(xmissing, ymissing),
      active = c(TRUE, TRUE),
      names = c(xname, yname),
      sets = c(1, 2)
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    xg <- h$xGifi[[1]][[1]]
    yg <- h$xGifi[[2]][[1]]
    return (structure (
      list(
        burt = crossprod (cbind(g[[1]][[1]]$basis, g[[2]][[1]]$basis)),
        objectscores = h$x,
        xtransform = postCorals (x, xg$transform),
        ytransform = postCorals (y, yg$transform),
        rhat = cor (cbind (xg$transform, yg$transform)),
        xweights = xg$weights,
        yweights = yg$weights,
        xscores = postCorals (x, xg$scores),
        yscores = postCorals (y, yg$scores),
        xdmeasure = crossprod (xg$scores),
        ydmeasure = crossprod (yg$scores),
        xquantifications = xg$quantifications,
        yquantifications = yg$quantifications,
        xloadings = crossprod (xg$transform, h$x),
        yloadings = crossprod (yg$transform, h$x),
        lambda = (crossprod (xg$scores) + crossprod (yg$scores)) / 2,
        ntel = h$ntel,
        f = h$f
      ),
      class = "corals"
    ))
  }

coranals <- function () {

}

morals <-
  function (x,
            y,
            xknots = knotsQ(x),
            yknots = knotsQ(y),
            xdegrees = 2,
            ydegrees = 2,
            xordinal = TRUE,
            yordinal = TRUE,
            xties = "s",
            yties = "s",
            xmissing = "m",
            ymissing = "m",
            xnames = colnames (x, do.NULL = FALSE),
            ynames = "Y",
            xactive = TRUE,
            xcopies = 1,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE) {
    npred <- ncol (x)
    nobs <- nrow (x)
    xdegrees <- reshape (xdegrees, npred)
    xordinal <- reshape (xordinal, npred)
    xties <- reshape (xties, npred)
    xmissing <- reshape (xmissing, npred)
    xactive <- reshape (xactive, npred)
    xcopies <- reshape (xcopies, npred)
    g <- makeGifi (
      data = cbind (x, y),
      knots = c (xknots, yknots),
      degrees = c (xdegrees, ydegrees),
      ordinal = c (xordinal, yordinal),
      sets =  c (rep(1, npred), 2),
      copies = c (xcopies, 1),
      ties = c (xties, yties),
      missing = c (xmissing, ymissing),
      active = c (xactive, TRUE),
      names = c (xnames, ynames)
    )
    h <- gifiEngine(
      gifi = g,
      ndim = 1,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    xhat <- matrix (0, nobs, 0)
    for (j in 1:npred)
      xhat <- cbind (xhat, h$xGifi[[1]][[j]]$transform)
    yhat <- h$xGifi[[2]][[1]]$transform
    rhat <- cor (cbind (xhat, yhat))
    qxy <- lsRC(xhat, yhat)$solution
    ypred <- xhat %*% qxy
    yres <- yhat - ypred
    smc <- sum (yhat * ypred)
    return (structure (
      list (
        objscores = h$x,
        xhat = xhat,
        yhat = yhat,
        rhat = rhat,
        beta = qxy,
        ypred = ypred,
        yres = yres,
        smc = smc,
        ntel = h$ntel,
        f = h$f
      ),
      class = "morals"
    ))
  }

princals <-
  function (data,
            knots = knotsQ (data),
            degrees = 2,
            ordinal = TRUE,
            copies = 1,
            ndim = 2,
            ties = "s",
            missing = "m",
            names = colnames (data, do.NULL = FALSE),
            active = TRUE,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE)  {
    aname <- deparse (substitute (data))
    nvars <- ncol (data)
    nobs <- nrow (data)
    g <- makeGifi (
      data = data,
      knots = knots,
      degrees = reshape (degrees, nvars),
      ordinal = reshape (ordinal, nvars),
      sets =  1:nvars,
      copies = reshape (copies, nvars),
      ties = reshape (ties, nvars),
      missing = reshape (missing, nvars),
      active = reshape (active, nvars),
      names = names
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    a <- v <- z <- d <- y <- o <- as.list (1:nvars)
    dsum <- matrix (0, ndim, ndim)
    for (j in 1:nvars) {
      jgifi <- h$xGifi[[j]][[1]]
      v[[j]] <- jgifi$transform
      a[[j]] <- jgifi$weights
      y[[j]] <- jgifi$scores
      z[[j]] <- jgifi$quantifications
      cy <- crossprod (y[[j]])
      dsum <- dsum + cy
      d[[j]] <- cy
      o[[j]] <- crossprod (h$x, v[[j]])
    }
    return (structure (
      list (
        transform = v,
        rhat = corList (v),
        objectscores = h$x,
        scores = y,
        quantifications = z,
        dmeasures = d,
        lambda = dsum / ncol (data),
        weights = a,
        loadings = o,
        ntel = h$ntel,
        f = h$f
      ),
      class = "princals"
    ))
  }

criminals <-
  function (x,
            y,
            xknots = knotsQ(x),
            yknots = knotsD(y),
            xdegrees = 2,
            ydegrees = -1,
            xordinal = TRUE,
            yordinal = FALSE,
            xcopies = 1,
            xties = "s",
            yties = "s",
            xmissing = "m",
            ymissing = "m",
            xactive = TRUE,
            xnames = colnames (x, do.NULL = FALSE),
            ynames = "Y",
            ndim = 2,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE) {
    aname <- deparse (substitute (data))
    npred <- ncol (x)
    nobs <- nrow (x)
    g <- makeGifi (
      data = cbind (x, y),
      knots = c(xknots, yknots),
      degrees = c (reshape (xdegrees, npred), ydegrees),
      ordinal = c (reshape (xordinal, npred), yordinal),
      sets =  c (rep(1, npred), 2),
      copies = c (reshape (xcopies, npred), length (unique (y))),
      ties = c (reshape (xties, npred), yties),
      missing = c (reshape (xmissing, npred), ymissing),
      active = c (reshape (xactive, npred), TRUE),
      names = c(xnames, ynames)
    )
    h <- gifiEngine(
      gifi = g,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    x <- matrix (0, nobs, 0)
    for (j in 1:npred)
      x <- cbind (x, h$xGifi[[1]][[j]]$transform)
    y <- as.vector(y)
    g <- ifelse (outer (y, unique (y), "=="), 1, 0)
    d <- colSums (g)
    v <- crossprod (x)
    u <- crossprod (g, x)
    b <- crossprod (u, (1 / d) * u)
    w <- v - b
    e <- eigen (v)
    k <- e$vectors
    l <- sqrt (abs (e$values))
    l <- ifelse (l < 1e-7, 0, 1 / l)
    f <- eigen (outer(l, l) * crossprod (k, b %*% k))
    a <- k %*% (l * f$vectors)
    z <- x %*% a
    u <- (1 / d) * crossprod (g, x)
    return (structure (
      list(
        objectscores = h$x,
        xhat = x,
        loadings = a,
        scores = z,
        groupmeans = u,
        bwratios = f$values,
        ntel = h$ntel,
        f = h$f
      ),
      class = "criminals"
    ))
  }

canals <-
  function (x,
            y,
            xknots = knotsQ(x),
            yknots = knotsQ(y),
            xdegrees = rep(2, ncol(x)),
            ydegrees = rep(2, ncol(y)),
            xordinal = rep (TRUE, ncol (x)),
            yordinal = rep (TRUE, ncol (y)),
            xcopies = rep (1, ncol (x)),
            ycopies = rep (1, ncol (y)),
            ndim = 2,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE) {
    h <- gifiEngine(
      data = cbind (x, y),
      knots = c(xknots, yknots),
      degrees = c(xdegrees, ydegrees),
      ordinal = c(xordinal, yordinal),
      sets = c(rep(1, ncol(x)), rep(2, ncol (y))),
      copies = c(xcopies, ycopies),
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    x <- h$h[, 1:sum(xcopies)]
    y <- h$h[, -(1:sum(xcopies))]
    u <- crossprod (x)
    v <- crossprod (y)
    w <- crossprod (x, y)
    a <- solve (chol (u))
    b <- solve (chol (v))
    s <- crossprod (a, w %*% b)
    r <- svd (s)
    xw <- a %*% (r$u)
    yw <- b %*% (r$v)
    xs <- x %*% xw
    ys <- y %*% yw
    xl <- crossprod (x, xs)
    yl <- crossprod (y, ys)
    return (structure (
      list(
        xhat = x,
        yhat = y,
        rhat = cor (cbind (x, y)),
        cancors = r$d,
        xweights = xw,
        yweights = yw,
        xscores = xs,
        yscores = ys,
        xloadings = xl,
        yloadings = yl,
        ntel = h$ntel,
        f = h$f
      ),
      class = "canals"
    ))
  }


overals <-
  function (data,
            sets,
            copies,
            knots = knotsQ (data),
            degrees = rep (2, ncol (data)),
            ordinal = rep (TRUE, ncol (data)),
            ndim = 2,
            itmax = 1000,
            eps = 1e-6,
            seed = 123,
            verbose = FALSE)  {
    h <- gifiEngine(
      data = data,
      knots = knots,
      degrees = degrees,
      ordinal = ordinal,
      sets = sets,
      copies = copies,
      ndim = ndim,
      itmax = itmax,
      eps = eps,
      seed = seed,
      verbose = verbose
    )
    xhat <- h$h
    rhat <- cor (xhat)
    a <- h$a
    y <- xhat
    for (j in 1:ncol(data)) {
      k <- (1:ndim) + (j - 1) * ndim
      y[, k] <- xhat[, k] %*% a[k,]
    }
    return (structure (
      list (
        xhat = xhat,
        rhat = rhat,
        objscores = h$x,
        quantifications = y,
        ntel = h$ntel,
        f = h$f
      ),
      class = "overals"
    ))
  }


primals <- function () {

}

addals <- function () {

}

pathals <- function () {

}