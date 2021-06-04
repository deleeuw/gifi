dyn.load("code/coding.so")

decode <- function(cell, dims) {
  if (length(cell) != length(dims)) {
    stop("Dimension error")
  }
  if (any(cell > dims) || any (cell < 1)) {
    stop("No such cell")
  }
  .Call("DECODE", as.integer(cell), as.integer(dims))
}


encode <- function(ind, dims) {
  if (length(ind) > 1) {
    stop ("Dimension error")
  }
  if ((ind < 1) || (ind > prod(dims))) {
    stop ("No such cell")
  }
  .Call("ENCODE", as.integer(ind), as.integer(dims))
}