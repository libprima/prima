// An example to illustrate the use of COBYLA.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdint.h>

static void fun(const double x[], double *f, double constr[], const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;
  constr[0] = x1*x1 + x2*x2 - 13;// ||x||^2<=13
  (void)data;
}

int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;
  const int m_nlcon = 1;
  const int n = 2;
  double x[2] = {0.0, 0.0};
  double f = 0.0;
  double cstrv = 0.0;
  double nlconstr[m_nlcon];
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
  int nf = 0;
  void *data = NULL;
  const int rc = prima_cobyla(m_nlcon, &fun, data, n, x, &f, &cstrv, nlconstr, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, iprint);
  const char *msg = prima_get_rc_string(rc);
  printf("x*={%g, %g} f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", x[0], x[1], f, cstrv, nlconstr[0], rc, msg, nf);
  return (fabs(x[0]-3)>2e-2 || fabs(x[1]-2)>2e-2);
}
