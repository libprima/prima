// An example to illustrate the use of NEWUOA.

#include <math.h>
#include <stdio.h>

// Make PRIMA available
#include "prima/prima.h"

// Objective function
static void fun(const double x[], double *f, const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = pow(x1-5, 2) + pow(x2-4, 2);
  (void)data;
}

// Callback function
static void callback(int n, const double x[], double f, int nf, int tr, double cstrv, int m_nlcon, const double nlconstr[], bool *terminate)
{
  (void)n;
  (void)cstrv;
  (void)m_nlcon;
  (void)nlconstr;
  printf("Best point so far: x = {%g, %g}, f = %g, nf = %d, tr = %d\n", x[0], x[1], f, nf, tr);
  *terminate = 0;
}

// Main function
int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;
  const int n = 2;
  double x0[2] = {0.0, 0.0};

  // Set up the problem
  prima_problem_t problem;
  prima_init_problem(&problem, n);
  problem.calfun = &fun;
  problem.x0 = x0;

  // Set up the options
  prima_options_t options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_EXIT;
  options.rhoend= 1e-3;
  options.maxfun = 200*n;
  options.callback = &callback;

  // Call the solver
  prima_result_t result;
  const int rc = prima_minimize(PRIMA_NEWUOA, &problem, &options, &result);

  // Print the result
  printf("x* = {%g, %g}, rc = %d, msg = '%s', evals = %d\n", result.x[0], result.x[1], rc, result.message, result.nf);

  // Check the result
  int success = (fabs(result.x[0] - 5) > 2e-2 || fabs(result.x[1] - 4) > 2e-2);

  // Free the result
  prima_free_result(&result);

  return success;
}
