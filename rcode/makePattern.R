makeGifiPattern <- function (gifiFrame, nequations, atype) {
  nvars <- length (gifiFrame)
  neqs <- length (nequations)
  ncopies <- sapply (1:nvars, function (j) gifiFrame[[j]]$copies)
  a <- matrixList (nvars, neqs, ncopies, nequations)
  for (i in 1:nvars) {
    for (j in 1:nac) {
      if ((atype[i, j] == "0") || (atype[i, j] == "1")) {
        a[[i]][[j]] <- NULL
      }
      else {
        a[[i]][[j]] <- matrix (rnorm (hrow[i] * hcol[j]), hrow[i], hcol[j])
      }
    }
  }
  return (a)
}

matrixList <- function (n, m, nr, nc) {
  a <- list ()
  ai <- as.list (1:m)
  for (i in 1:n) {
    a <- c (a, list (ai))
    for (j in 1:m) {
      a[[i]][[j]] <- matrix (0, nr[i], nc[j])
    }
  }
  return (a)
}
