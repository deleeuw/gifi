# an object of class gifiVariable contains information about the variable that
# does not change during computation -- it stores the input data and parameters

makeGifiVariable <-
  function (data,
            weights,
            knots,
            degree,
            ordinal,
            ties,
            copies,
            missing,
            active,
            name) {
    there <- which (!is.na (data))
    notthere <- which (is.na (data))
    nmis <- length (notthere)
    nobs <- length (data)
    if (length (there) == 0)
      stop ("a gifiVariable cannot be completely missing")
# make a basis matrix for the nonmissing data
    work <- data[there]
    if (degree == -2) {
      type <- "orthoblock"
      basis <- NULL
    }
    if (degree == -1) {
      type <- "categorical"
      basis <- makeIndicator (work)
      if (ncol (basis) == 1) {
        stop ("a gifiVariable must have more than one category")
      }
      if (ncol (basis) == 2) {
        type <- "binary"
      }
    }
    if (degree >= 0) {
      if (length (knots) == 0)
        type <- "polynomial"
      else
        type <- "splinical"
      basis <- bsplineBasis (work, degree, knots)
    }
# make the basis complete by adding rows for missing data
    if (nmis > 0)
      basis <- makeMissing (data, basis, missing)
# correct for too many copies
    copies <- min (copies, ncol (basis) - 1)
# QR decomposition of basis
    qr <- gsRC (center (basis))
    if (qr$rank == 0)
      stop ("a gifiVariable cannot be completely zero")
    return (structure (
      list (
        data = data,
        basis = basis,
        qr = qr,
        copies = copies,
        degree = degree,
        ties = ties,
        missing = missing,
        ordinal = ordinal,
        active = active,
        name = name,
        type = type
      ),
      class = "gifiVariable"
    ))
  }

# an object of class gifiSet is a list of gifiVariable objects

makeGifiSet <-
  function (data,
            weights,
            knots,
            degrees,
            ordinal,
            ties,
            copies,
            missing,
            active,
            names) {
    nvars <- ncol (data)
    varList <- as.list (1:nvars)
    for (i in 1:nvars) {
      varList[[i]] <-
        makeGifiVariable (
          data = data[, i],
          weights = weights[, i],
          knots = knots[[i]],
          degree = degrees[i],
          ordinal = ordinal[i],
          ties = ties[i],
          copies = copies[i],
          missing = missing[i],
          active = active[i],
          name = names[i]
        )
    }
    return (structure (varList, class = "gifiSet"))
  }

# an object of class gifi is a list of objects of class gifiSet

makeGifi <-
  function (data,
            weights,
            knots,
            degrees,
            ordinal,
            ties,
            copies,
            missing,
            active,
            names,
            sets) {
    nsets <- max (sets)
    setList <- as.list (1:nsets)
    for (i in 1:nsets) {
      k <- which (sets == i)
      setList [[i]] <-
        makeGifiSet (
          data = data[, k, drop = FALSE],
          weights = weights[, k],
          knots = knots[k],
          degrees = degrees[k],
          ordinal = ordinal[k],
          ties = ties[k],
          copies = copies[k],
          missing = missing[k],
          active = active[k],
          names = names[k]
        )
    }
    return (structure (setList, class = "gifi"))
  }

# an object of class xGifiVariable contains information about the variable that
# changes during computation -- it stores the initial estimates, which will
# become the eventual output

xGifiVariable <- function (gifiVariable, x) {
  ndim <- ncol (x)
  basis <- gifiVariable$basis
  nbas <- ncol (basis)
  nobs <- length (gifiVariable$data)
  copies <- gifiVariable$copies
  transform <- matrix (0, nobs, copies)
  transform[, 1] <- drop(basis %*% (1:nbas))
  if (copies > 1) {
    for (i in 2:copies)
      transform[, i] <- drop (basis %*% rnorm (nbas))
  }
  transform <- gsRC (normalize (center (transform)))$q
  weights <- lsRC (transform, x)$solution
  scores <- transform %*% weights
  quantifications <- lsRC (basis, scores)$solution
  return (structure (
    list(
      transform = transform,
      weights = weights,
      scores = scores,
      quantifications = quantifications
    ),
    class = "xGifiVariable"
  ))
}

# an object of class xGifiSet is a list of objects of class xGifiVariable

xGifiSet <- function (gifiSet, x) {
  nvars <- length (gifiSet)
  varList <- as.list (1:nvars)
  for (i in 1:nvars) {
    varList[[i]] <- xGifiVariable (gifiSet[[i]], x)
  }
  return (structure (varList, class = "xGifiSet"))
}

# an object of class xGifi is a list of objects of class xGifiSet

xGifi <- function (gifi, x) {
  nsets <- length (gifi)
  setList <- as.list (1:nsets)
  for (i in 1:nsets) {
    setList[[i]] <- xGifiSet (gifi[[i]], x)
  }
  return (structure (setList, class = "xGifi"))
}
