// Test for data callback argument.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

const int n = 2;
#define m_nlcon 1
int debug = 0;
static int int_data = 0xff;
const void * data_ref = &int_data;

static void fun(const double x[], double *f, const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;

  static int count = 0;
  if (debug)
  {
    ++ count;
    printf("count=%d\n", count);
  }

  // check data is ok
  if (data != data_ref)
  {
    printf("invalid data\n");
    *f = NAN;
  }
}

static void fun_con(const double x[], double *f, double constr[], const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;
  constr[0] = x1*x1 + x2*x2 - 13;// ||x||^2<=13

  static int count = 0;
  if (debug)
  {
    ++ count;
    printf("count=%d\n", count);
  }

  // check data is ok
  if (data != data_ref)
  {
    printf("invalid data\n");
    *f = NAN;
  }
}

int main(int argc, char * argv[])
{
  char *algo = "bobyqa";
  if (argc > 1)
    algo = argv[1];
  printf("algo=%s\n", algo);

  if (argc > 2)
    debug = (strcmp(argv[2], "debug") == 0);
  printf("debug=%d\n", debug);

  double x[] = {0, 0};
  double xl[] = {-6.0, -6.0};
  double xu[] = {6.0, 6.0};
  double f = 0.0;
  double cstrv = 0.0;
  double nlconstr[m_nlcon];
  for (int j = 0; j < m_nlcon; ++ j)
    nlconstr[j] = 0.0;
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
  const double rhobeg = 1.0;
  const double rhoend = 1e-6;
  const double ftarget = -INFINITY;
  const int iprint = PRIMA_MSG_RHO;
  const int maxfun = 500*n;
  int nf = 0;
  int rc = 0;
  if(strcmp(algo, "bobyqa") == 0)
  {
    rc = prima_bobyqa(&fun, data_ref, n, x, &f, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "cobyla") == 0)
  {
    rc = prima_cobyla(m_nlcon, &fun_con, data_ref, n, x, &f, &cstrv, nlconstr, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, iprint);
  }
  else if(strcmp(algo, "lincoa") == 0)
  {
    rc = prima_lincoa(&fun, data_ref, n, x, &f, &cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "newuoa") == 0)
  {
    rc = prima_newuoa(&fun, data_ref, n, x, &f, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "uobyqa") == 0)
  {
    rc = prima_uobyqa(&fun, data_ref, n, x, &f, &nf, rhobeg, rhoend, ftarget, maxfun, iprint);
  }
  else
  {
    printf("incorrect algo\n");
    return 1;
  }
  const char *msg = prima_get_rc_string(rc);

  printf("f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", f, cstrv, nlconstr[0], rc, msg, nf);
  return (fabs(x[0]-3)>2e-2 || fabs(x[1]-2)>2e-2);
}
