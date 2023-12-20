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
  constr[0] = x1*x1 + x2*x2 - 12;// ||x||^2<=12
  (void)data;
}


static void callback(const int n, const double x[], const double f, const int nf, const int tr, const double cstrv, const int m_nlcon, const double nlconstr[], bool *terminate)
{
  (void)n;
  (void)m_nlcon;
  printf("best point so far: x=[%g;%g] f=%g cstrv=%g nlconstr=%g nf=%d tr=%d\n", x[0], x[1], f, cstrv, nlconstr[0], nf, tr);
  *terminate = 0;
}


int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;
  const int n = 2;
  double x0[2] = {0.0, 0.0};
  // set up the problem
  prima_problem_t problem;
  prima_init_problem(&problem, n);
  problem.calcfc = &fun;
  problem.x0 = x0;
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
  // set up the options
  prima_options_t options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_EXIT;
  options.rhoend= 1e-3;
  options.maxfun = 200*n;
  options.callback = &callback;
  // initialize the result
  prima_result_t result;
  // run the solver
  const int rc = prima_minimize(PRIMA_COBYLA, &problem, &options, &result);
  printf("x*={%g, %g} f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", result.x[0], result.x[1], result.f, result.cstrv, result.nlconstr[0], rc, result.message, result.nf);
  prima_free_problem(&problem);
  prima_free_result(&result);
  return (fabs(result.x[0]-2.86)>2e-2 || fabs(result.x[1]-1.94)>2e-2);
}
