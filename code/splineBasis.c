#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

double bs (int nknots, int nspline, int degree, double x, double * knots);
int mindex (int i, int j, int nrow);

void splinebasis (int *d, int *n, int *m, double * x, double * knots, double * basis) {
    int mm = *m, dd = *d, nn = *n;
    int k = mm - dd - 1, i , j, ir, jr;
    for (i = 0; i < nn; i++) {
        ir = i + 1;
        if (x[i] == knots[mm - 1]) {
           basis [mindex (ir, k, nn) - 1] =  1.0;
            for (j = 0; j <  (k - 1);  j++) {
                jr = j + 1;
                basis [mindex (ir, jr, nn) - 1] = 0.0;
            }
        } else {
            for (j = 0; j < k ; j++) {
                jr = j + 1;
                basis [mindex (ir, jr, nn) - 1] = bs (mm, jr, dd + 1, x[i], knots);
            }
        }
    }
}

int mindex (int i, int j, int nrow) {
    return (j - 1) * nrow + i;
}

double bs (int nknots, int nspline, int updegree, double x, double * knots) {
    double y, y1, y2, temp1, temp2;
     if (updegree == 1) {
        if ((x >= knots[nspline - 1]) && (x < knots[nspline]))
            y = 1.0;
        else
            y = 0.0;
     }
     else {
        temp1 = 0.0;
        if ((knots[nspline + updegree - 2] - knots[nspline - 1]) > 0)
            temp1 = (x - knots[nspline - 1]) / (knots[nspline + updegree - 2] - knots[nspline - 1]);
        temp2 = 0.0;
        if ((knots[nspline + updegree - 1] - knots[nspline]) > 0)
            temp2 = (knots[nspline + updegree - 1] - x) / (knots[nspline + updegree - 1] - knots[nspline]);
        y1 = bs(nknots, nspline, updegree - 1, x, knots);
        y2 = bs(nknots, nspline + 1, updegree - 1, x, knots);
        y =  temp1 * y1 + temp2 * y2;
    }
    return y;
}
