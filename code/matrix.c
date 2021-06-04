#include <math.h>

#define MINDEX(i, j, n) (((j)-1) * (n) + (i)-1)

void gsC(double *, double *, int *, int *, int *, int *, double *);

void gsC(double *x, double *r, int *n, int *m, int *rank, int *pivot,
         double *eps) {
  int i, j, ip, nn = *n, mm = *m, rk = *m, jwork = 1;
  double s = 0.0, p;
  for (j = 1; j <= mm; j++) {
    pivot[j - 1] = j;
  }
  while (jwork <= rk) {
    for (j = 1; j < jwork; j++) {
      s = 0.0;
      for (i = 1; i <= nn; i++) {
        s += x[MINDEX(i, jwork, nn)] * x[MINDEX(i, j, nn)];
      }
      r[MINDEX(j, jwork, mm)] = s;
      for (i = 1; i <= nn; i++) {
        x[MINDEX(i, jwork, nn)] -= s * x[MINDEX(i, j, nn)];
      }
    }
    s = 0.0;
    for (i = 1; i <= nn; i++) {
      s += x[MINDEX(i, jwork, nn)] * x[MINDEX(i, jwork, nn)];
    }
    if (s > *eps) {
      s = sqrt(s);
      r[MINDEX(jwork, jwork, mm)] = s;
      for (i = 1; i <= nn; i++) {
        x[MINDEX(i, jwork, nn)] /= s;
      }
      jwork += 1;
    }
    if (s <= *eps) {
      ip = pivot [rk - 1];
      pivot[rk - 1] = pivot[jwork - 1];
      pivot[jwork - 1] = ip;
      for (i = 1; i <= nn; i++) {
        p = x[MINDEX(i, rk, nn)];
        x[MINDEX(i, rk, nn)] = x[MINDEX(i, jwork, nn)];
        x[MINDEX(i, jwork, nn)] = p;
      }
      for (j = 1; j <= mm; j++) {
        p = r[MINDEX(j, rk, mm)];
        r[MINDEX(j, rk, mm)] = r[MINDEX(j, jwork, mm)];
        r[MINDEX(j, jwork, mm)] = p;
      }
      rk -= 1;
    }
  }
  *rank = rk;
}
