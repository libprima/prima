// A stress test on excessively large problems.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define n_max 2000
int n = 0;
#define m_ineq_max 1000
int m_ineq = 0;
#define m_nlcon 200
const double alpha = 4.0;
int debug = 0;

static double random_gen(double a, double b)
{
  return a + rand() * (b - a) / RAND_MAX;
}

static void fun(const double x[], double *f, const void *data)
{
  // Rosenbrock function
  *f = 0.0;
  for (int i = 0; i < n-1; ++ i)
    *f += (x[i] - 1.0) * (x[i] - 1.0) + alpha * (x[i+1] - x[i]*x[i]) * (x[i+1] - x[i]*x[i]);

  static int count = 0;
  if (debug)
  {
    ++ count;
    printf("count=%d\n", count);
  }
  (void)data;
}

static void fun_con(const double x[], double *f, double constr[], const void *data)
{
  // Rosenbrock function
  *f = 0.0;
  for (int i = 0; i < n-1; ++ i)
    *f += (x[i] - 1.0) * (x[i] - 1.0) + alpha * (x[i+1] - x[i]*x[i]) * (x[i+1] - x[i]*x[i]);
  // x_{i+1} <= x_i^2
  for (int i = 0; i < MIN(m_nlcon, n-1); ++ i)
    constr[i] = x[i+1] - x[i] * x[i];

  static int count = 0;
  if (debug)
  {
    ++ count;
    printf("count=%d\n", count);
  }
  (void)data;
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

  // set seed to year/week
  char buf[10] = {0};
  time_t t = time(NULL);
  struct tm *tmp = localtime(&t);
  int rc = strftime(buf, 10, "%y%W", tmp);
  if (!rc)
    return 1;
  unsigned seed = atoi(buf);
  printf("seed=%d\n", seed);
  srand(seed);

  double x[n_max];
  double xl[n_max];
  double xu[n_max];
  prima_options options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_RHO;
  options.maxfun = 500*n_max;
  double *Aineq = malloc(n_max*m_ineq_max*sizeof(double));
  double bineq[m_ineq_max];
  options.Aineq = Aineq;
  options.bineq = bineq;
  options.xl = xl;
  options.xu = xu;
  for (int j = 0; j < m_ineq_max; ++ j)
    bineq[j] = random_gen(-1.0, 1.0);

  for (int i = 0; i < n_max; ++ i)
  {
    for (int j = 0; j < m_ineq; ++ j)
      Aineq[j*n_max+i] = random_gen(-1.0, 1.0);
    x[i] = random_gen(-1.0, 1.0);
    xl[i] = -1.0;
    xu[i] = 1.0;
  }
  prima_results results;
  if(strcmp(algo, "bobyqa") == 0)
  {
    n = 1600;
    rc = prima_bobyqa(&fun, n, x, &options, &results);
  }
  else if(strcmp(algo, "cobyla") == 0)
  {
    n = 800;
    options.m_nlcon = m_nlcon;
    options.m_ineq = 600;
    rc = prima_cobyla(&fun_con, n, x, &options, &results);
  }
  else if(strcmp(algo, "lincoa") == 0)
  {
    n = 1000;
    options.m_ineq = 1000;
    rc = prima_lincoa(&fun, n, x, &options, &results);
  }
  else if(strcmp(algo, "newuoa") == 0)
  {
    n = 1600;
    rc = prima_newuoa(&fun, n, x, &options, &results);
  }
  else if(strcmp(algo, "uobyqa") == 0)
  {
    n = 100;
    rc = prima_uobyqa(&fun, n, x, &options, &results);
  }
  else
  {
    printf("incorrect algo\n");
    return 1;
  }
  const char *msg = prima_get_rc_string(rc);

  printf("f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", results.f, results.cstrv, results.nlconstr ? results.nlconstr[0] : 0.0, rc, msg, results.nf);
  prima_free_options(&options);
  prima_free_results(&results);
  return 0;
}
