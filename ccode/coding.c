#include <R.h>
#include <Rinternals.h>

SEXP DECODE( SEXP, SEXP );
SEXP ENCODE( SEXP, SEXP );

SEXP
  DECODE( SEXP cell, SEXP dims )
  {
    int             aux = 1, n = length(dims);
    SEXP            ind;
    PROTECT( ind = allocVector( INTSXP, 1 ) );
    INTEGER( ind )[0] = 1;
    for( int i = 0; i < n; i++ ) {
      INTEGER( ind )[0] += aux * ( INTEGER( cell )[i] - 1 );
      aux *= INTEGER(dims)[i];
    }
    UNPROTECT( 1 );
    return (ind);
  }

SEXP
  ENCODE( SEXP ind, SEXP dims )
  {
    int             n = length(dims), aux = INTEGER(ind)[0], pdim = 1;
    SEXP            cell;
    PROTECT( cell = allocVector( INTSXP, n ) );
    for ( int i = 0; i < n - 1; i++ )
      pdim *= INTEGER( dims )[i];
    for ( int i = n - 1; i > 0; i-- ){
      INTEGER( cell )[i] = ( aux - 1 ) / pdim;
      aux -= pdim * INTEGER( cell )[i];
      pdim /= INTEGER( dims )[i - 1];
      INTEGER( cell )[i] += 1;
    }
    INTEGER( cell )[0] = aux;
    UNPROTECT( 1 );
    return cell;
  }