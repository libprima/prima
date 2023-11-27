// A stress test on excessively large problems.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define N_MAX 2000
#define M_INEQ_MAX 1000
#define M_NLCON 200

int n = 0;
int m_ineq = 0;
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
  for (int i = 0; i < MIN(M_NLCON, n-1); ++ i)
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

  double x0[N_MAX];
  double xl[N_MAX];
  double xu[N_MAX];
  prima_problem problem;
  prima_init_problem(&problem, N_MAX);
  problem.x0 = x0;
  problem.calcfc = &fun_con;
  problem.calfun = &fun;
  prima_options options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_RHO;
  options.maxfun = 500*N_MAX;
  double *Aineq = malloc(N_MAX*M_INEQ_MAX*sizeof(double));
  double bineq[M_INEQ_MAX];
  problem.Aineq = Aineq;
  problem.bineq = bineq;
  problem.xl = xl;
  problem.xu = xu;
  for (int j = 0; j < M_INEQ_MAX; ++ j)
    bineq[j] = random_gen(-1.0, 1.0);

  for (int i = 0; i < N_MAX; ++ i)
  {
    for (int j = 0; j < m_ineq; ++ j)
      Aineq[j*N_MAX+i] = random_gen(-1.0, 1.0);
    x0[i] = random_gen(-1.0, 1.0);
    xl[i] = -1.0;
    xu[i] = 1.0;
  }
  int algorithm = 0;
  prima_result result;
  if(strcmp(algo, "bobyqa") == 0)
  {
    algorithm = PRIMA_BOBYQA;
    problem.n = 1600;
  }
  else if(strcmp(algo, "cobyla") == 0)
  {
    algorithm = PRIMA_COBYLA;
    problem.n = 800;
    problem.m_nlcon = M_NLCON;
    problem.m_ineq = 600;
  }
  else if(strcmp(algo, "lincoa") == 0)
  {
    algorithm = PRIMA_LINCOA;
    problem.n = 1000;
    problem.m_ineq = 1000;
  }
  else if(strcmp(algo, "newuoa") == 0)
  {
    algorithm = PRIMA_NEWUOA;
    problem.n = 1600;
  }
  else if(strcmp(algo, "uobyqa") == 0)
  {
    algorithm = PRIMA_UOBYQA;
    problem.n = 100;
  }
  else
  {
    printf("incorrect algo\n");
    return 1;
  }
  rc = prima_minimize(algorithm, &problem, &options, &result);
  printf("f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", result.f, result.cstrv, result.nlconstr ? result.nlconstr[0] : 0.0, rc, result.message, result.nf);
  prima_free_problem(&problem);
  prima_free_result(&result);
  return 0;
}
