maxeig <- function (r, p) {
  e <- eigen (r)
  f <- sum (e$values[1:p])
  g <- tcrossprod(e$vectors[,1:p])
  return (list (f = f, g = g))
}

maxcor <- function (r, p) {
  f <- sum (r ^ p)
  g <- p * (r ^ (p - 1))
  return (list (f = f, g = g))
}

maxabs <- function (r, p) {
  f <- sum (abs(r) ^ p)
  g <- p * (abs(r) ^ (p - 1)) * sign(r)
  return (list (f = f, g = g))
}

maxdet <- function (r) {
  f <- -log(det (r))
  g <- -solve(r)
  return (list (f = f, g = g))
}

maxsmc <- function (r, p) {
  beta <- solve (r[-p,-p], r[p,-p])
  f <- sum (beta * r[p,-p])
  h <- rep (1, nrow (r))
  h[-p] <- -beta
  g <- -outer (h, h)
  return (list (f = f, g = g))
}

maxsum <- function (r, p) {
  f <- sum (sqrt (r ^ 2 + p))
  g <- r / sqrt (r ^ 2 + p)
  return (list (f = f, g = g))
}

maximage <- function (r) {
  n <- nrow(r)
  f <- 0
  g <- matrix (0, n, n)
  for (p in 1:n) {
    beta <- solve (r[-p,-p], r[p,-p])
    f <- f + sum (beta * r[p,-p])
    h <- rep (1, nrow (r))
    h[-p] <- -beta
    g <- g - outer (h, h)
  }
  return (list (f = f, g = g))
}

maxfac <- function (r, p) {
  fa <- factanal (NULL, p, covmat = r, rotation = "none")
  s <- tcrossprod (fa$loadings) + diag (fa$unique)
  g <- - solve (s)
  f <- -log(det (s)) + sum (g * r)
  return (list (f = f, g = g))
}
