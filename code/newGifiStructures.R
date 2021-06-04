# an object of class gifiVariable contains information about the variable that
# does not change during computation -- it stores the input data and parameters

makeGifiVariable <-
  function (data,
            knots,
            degree,
            ordinal,
            ties,
            copies,
            missing,
            active,
            name) {
    if (is.NULL (data)) {
      type <- "latent"
      basis <- diag (nobs)
      qr = gsRC (center (basis))
    }
    else {
      there <- which (!is.na (data))
      notthere <- which (is.na (data))
      nmis <- length (notthere)
      nobs <- length (data)
      work <- data[there]
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
      if ((nmis > 0) && (degree > -2))
        basis <- makeMissing (data, basis, missing)
      copies <- min (copies, ncol (basis) - 1)
      qr <- gsRC (center (basis))
      if (qr$rank == 0)
        stop ("a gifiVariable cannot be completely zero")
    }
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

# an object of class gifiFrame is a list of objects of class gifiVariable

makeGifiFrame <-
  function (data,
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
      varList [[i]] <-
        makeGifiVariable (
          data = data[, i, drop = FALSE],
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
    return (structure (varList, class = "gifiFrame"))
  }

# an object of class xGifiVariable contains information about the variable that
# changes during computation -- it stores the initial estimates, which will
# become the eventual output

xGifiVariable <- function (gifiVariable) {
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
  quantifications <- lsRC (basis, transform)$solution
  return (structure (
    list(transform = transform,
         quantifications = quantifications),
    class = "xGifiVariable"
  ))
}


# an object of class xGifiFrame is a list of objects of class xGifiVariable

xGifiFrame <- function (gifiFrame) {
  nvars <- length (gifiFrame)
  varList <- as.list (1:nvars)
  for (i in 1:nvars) {
    varList[[i]] <- xGifiVariable (gifiFrame[[i]])
  }
  return (structure (varList, class = "xGifiFrame"))
}
