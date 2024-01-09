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
  *f = pow(x1-5, 2) + pow(x2-4, 2);
  // We add a constraint we know will be active in order to demonstrate usage
  // The constraint is x(1)**2 - 9 <= 0, meaning |x1| <= 3.
  constr[0] = pow(x1, 2) - 9;
  (void)data;
}


static void callback(const int n, const double x[], const double f, const int nf, const int tr, const double cstrv, const int m_nlcon, const double nlconstr[], bool *terminate)
{
  (void)n;
  (void)cstrv;
  (void)m_nlcon;
  printf("best point so far: x=[%g;%g] f=%g nlconstr=%g nf=%d tr=%d\n", x[0], x[1], f, nlconstr[0], nf, tr);
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
  return (fabs(result.x[0]-3)>2e-2 || fabs(result.x[1]-4)>2e-2);
}
