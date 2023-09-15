// An example to illustrate the use of LINCOA.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>

static void fun(const double x[], double *f, const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;
  (void)data;
}

static void callback(int n, const double x[], double f, int nf, int tr, double cstrv, int m_nlcon, const double nlconstr[], _Bool *terminate)
{
  (void)n;
  printf("progress: x=[%g;%g] f=%g cstrv=%g nf=%d tr=%d\n", x[0], x[1], f, cstrv, nf, tr);
  (void)m_nlcon;
  (void)nlconstr;
  *terminate = 0;
}

int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;
  const int n = 2;
  double x[2] = {0.0, 0.0};
  double f = 0.0;
  double cstrv = 0.0;
  // x1<=4, x2<=3, x1+x2<=10
  const int m_ineq = 3;
  double Aineq[3*2] = {1.0, 0.0,
                       0.0, 1.0,
                       1.0, 1.0};
  double bineq[3] = {4.0,
                     3.0,
                     10.0};
  const int m_eq = 0;
  double *Aeq = NULL;
  double *beq = NULL;
  double xl[2] = {-6.0, -6.0};
  double xu[2] = {6.0, 6.0};
  const double rhobeg = 1.0;
  const double rhoend = 1e-3;
  const double ftarget = -INFINITY;
  const int iprint = PRIMA_MSG_EXIT;
  const int maxfun = 200*n;
  const int npt = 2*n+1;
  int nf = 0;
  void *data = NULL;
  const int rc = prima_lincoa(&fun, data, n, x, &f, &cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, &callback);
  const char *msg = prima_get_rc_string(rc);
  printf("x*={%g, %g} f*=%g cstrv=%g rc=%d msg='%s' evals=%d\n", x[0], x[1], f, cstrv, rc, msg, nf);
  return (fabs(x[0]-3)>2e-2 || fabs(x[1]-2)>2e-2);
}
