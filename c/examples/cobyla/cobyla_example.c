// An example to illustrate the use of COBYLA.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdint.h>

#define M_NLCON 1

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
  const int n = 2;
  double x0[2] = {0.0, 0.0};
  prima_problem problem;
  prima_init_problem(&problem, n);
  problem.calcfc = &fun;
  problem.x0 = x0;
  prima_options options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_EXIT;
  options.rhoend= 1e-3;
  options.maxfun = 200*n;
  problem.m_nlcon = M_NLCON;
  // x1<=4, x2<=3, x1+x2<=10
  problem.m_ineq = 3;
  double Aineq[3*2] = {1.0, 0.0,
                       0.0, 1.0,
                       1.0, 1.0};
  double bineq[3] = {4.0,
                     3.0,
                     10.0};
  problem.Aineq = Aineq;
  problem.bineq = bineq;
  double xl[2] = {-6.0, -6.0};
  double xu[2] = {6.0, 6.0};
  problem.xl = xl;
  problem.xu = xu;
  prima_result result;
  const int rc = prima_minimize(PRIMA_COBYLA, &problem, &options, &result);
  printf("x*={%g, %g} f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", result.x[0], result.x[1], result.f, result.cstrv, result.nlconstr[0], rc, result.message, result.nf);
  prima_free_problem(&problem);
  prima_free_result(&result);
  return (fabs(result.x[0]-3)>2e-2 || fabs(result.x[1]-2)>2e-2);
}
