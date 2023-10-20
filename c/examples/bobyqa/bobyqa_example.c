// An example to illustrate the use of BOBYQA.

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

int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;
  const int n = 2;
  double x[2] = {0.0, 0.0};
  double f = 0.0;
  double xl[2] = {-6.0, -6.0};
  double xu[2] = {6.0, 6.0};
  
  prima_options options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_EXIT;
  options.rhoend= 1e-3;
  options.maxfun = 200*n;

  int nf = 0;
  const int rc = prima_bobyqa(&fun, n, x, &f, xl, xu, &nf, &options);
  const char *msg = prima_get_rc_string(rc);
  printf("x*={%g, %g} rc=%d msg='%s' evals=%d\n", x[0], x[1], rc, msg, nf);
  return (fabs(x[0]-3)>2e-2 || fabs(x[1]-2)>2e-2);
}
