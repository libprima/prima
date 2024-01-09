// An example to illustrate the use of LINCOA.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>

static void fun(const double x[], double *f, const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = pow(x1-5, 2) + pow(x2-4, 2);
  (void)data;
}

static void callback(int n, const double x[], double f, int nf, int tr, double cstrv, int m_nlcon, const double nlconstr[], bool *terminate)
{
  (void)n;
  (void)m_nlcon;
  (void)nlconstr;
  printf("best point so far: x=[%g;%g] f=%g cstrv=%g nf=%d tr=%d\n", x[0], x[1], f, cstrv, nf, tr);
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
  problem.calfun = &fun;
  problem.x0 = x0;
  // Define the constraints. We define constraints that will be active
  // in order to demonstrate their usage. The constraint is x1 + x2 < 8
  problem.m_ineq = 1;
  double Aineq[1*2] = {1.0, 1.0};
  double bineq[1] = {8.0};
  problem.Aineq = Aineq;
  problem.bineq = bineq;
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
  const int rc = prima_minimize(PRIMA_LINCOA, &problem, &options, &result);
  printf("x*={%g, %g} f*=%g cstrv=%g rc=%d msg='%s' evals=%d\n", result.x[0], result.x[1], result.f, result.cstrv, rc, result.message, result.nf);
  prima_free_problem(&problem);
  prima_free_result(&result);
  return (fabs(result.x[0]-4.5)>2e-2 || fabs(result.x[1]-3.5)>2e-2);
}
