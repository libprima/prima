// A stress test on excessively large problems.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int n = 10;
const int m_nlcon = 200;
const double alpha = 4.0;
int debug = 0;
char *algo = "bobyqa";

static double random_gen(double a, double b)
{
  return a + rand() * (b - a) / RAND_MAX;
}

static void fun_inner(const double x[], double *f)
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
}

static void fun_outer(const double x_global[], double *f_global)
{
  int rc = 0;
  double f = 0.0;
  double cstrv = 0.0;
  double x[n];
  double xl[n];
  double xu[n];
  int nf = 0;
  const int m_ineq = 5;
  double *Aineq = malloc(n*m_ineq*sizeof(double));
  double bineq[m_ineq];
  const int m_eq = 0;
  double *Aeq = NULL;
  double *beq = NULL;
  const double rhobeg = 1.0;
  const double rhoend = 1e-6;
  const double ftarget = -INFINITY;
  const int iprint = 0;
  const int maxfun = 500*n;
  for (int i = 0; i < n; ++ i)
  {
    for (int j = 0; j < m_ineq; ++ j)
      Aineq[j*n+i] = random_gen(-1.0, 1.0);
    x[i] = x_global[i];
    xl[i] = -1.0;
    xu[i] = 1.0;
  }
  if(strcmp(algo, "bobyqa") == 0)
  {
    rc = prima_bobyqa(&fun_inner, n, x, &f, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "lincoa") == 0)
  {
    rc = prima_lincoa(&fun_inner, n, x, &f, &cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "newuoa") == 0)
  {
    rc = prima_newuoa(&fun_inner, n, x, &f, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "uobyqa") == 0)
  {
    rc = prima_uobyqa(&fun_inner, n, x, &f, &nf, rhobeg, rhoend, ftarget, maxfun, iprint);
  }
  (void)rc;
  *f_global = f;
}


int main(int argc, char * argv[])
{
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

  double x[n];
  double xl[n];
  double xu[n];
  double f = 0.0;
  double cstrv = 0.0;
  double nlconstr[m_nlcon];
  const int m_ineq = 5;
  double *Aineq = malloc(n*m_ineq*sizeof(double));
  double bineq[m_ineq];
  const int m_eq = 0;
  double *Aeq = NULL;
  double *beq = NULL;
  const double rhobeg = 1.0;
  const double rhoend = 1e-6;
  const double ftarget = -INFINITY;
  const int iprint = PRIMA_MSG_RHO;
  const int maxfun = 500*n;

  for (int j = 0; j < m_ineq; ++ j)
    bineq[j] = random_gen(-1.0, 1.0);
  for (int j = 0; j < m_nlcon; ++ j)
    nlconstr[j] = 0.0;

  for (int i = 0; i < n; ++ i)
  {
    for (int j = 0; j < m_ineq; ++ j)
      Aineq[j*n+i] = random_gen(-1.0, 1.0);
    x[i] = random_gen(-1.0, 1.0);
    xl[i] = -1.0;
    xu[i] = 1.0;
  }
  int nf = 0;
  if(strcmp(algo, "bobyqa") == 0)
  {
    rc = prima_bobyqa(&fun_outer, n, x, &f, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "cobyla") == 0)
  {
    rc = 0;
  }
  else if(strcmp(algo, "lincoa") == 0)
  {
    rc = prima_lincoa(&fun_outer, n, x, &f, &cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "newuoa") == 0)
  {
    rc = prima_newuoa(&fun_outer, n, x, &f, &nf, rhobeg, rhoend, ftarget, maxfun, 2*n+1, iprint);
  }
  else if(strcmp(algo, "uobyqa") == 0)
  {
    rc = prima_uobyqa(&fun_outer, n, x, &f, &nf, rhobeg, rhoend, ftarget, maxfun, iprint);
  }
  else
  {
    printf("incorrect algo\n");
    return 1;
  }
  const char *msg = prima_get_rc_string(rc);

  printf("f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", f, cstrv, nlconstr[0], rc, msg, nf);
  return 0;
}
