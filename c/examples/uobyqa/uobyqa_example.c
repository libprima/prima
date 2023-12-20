// An example to illustrate the use of UOBYQA.

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

static void callback(int n, const double x[], double f, int nf, int tr, double cstrv, const int m_nlcon, const double nlconstr[], bool *terminate)
{
  (void)n;
  (void)cstrv;
  (void)m_nlcon;
  (void)nlconstr;
  printf("best point so far: x=[%g;%g] f=%g nf=%d tr=%d\n", x[0], x[1], f, nf, tr);
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
  const int rc = prima_minimize(PRIMA_UOBYQA, &problem, &options, &result);
  printf("x*={%g, %g} rc=%d msg='%s' evals=%d\n", result.x[0], result.x[1], rc, result.message, result.nf);
  prima_free_problem(&problem);
  prima_free_result(&result);
  return (fabs(result.x[0]-3)>2e-2 || fabs(result.x[1]-2)>2e-2);
}
