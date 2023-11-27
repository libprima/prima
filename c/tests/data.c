// Test for data callback argument.

#include "prima/prima.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define M_NLCON 1
const int n = 2;
int debug = 0;
static int int_data = 0xff;
void * data_ref = &int_data;

static void fun(const double x[], double *f, const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;

  static int count = 0;
  if (debug)
  {
    ++ count;
    printf("count=%d\n", count);
  }

  // check data is ok
  if (data != data_ref)
  {
    printf("invalid data\n");
    *f = NAN;
  }
}

static void fun_con(const double x[], double *f, double constr[], const void *data)
{
  const double x1 = x[0];
  const double x2 = x[1];
  *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;
  constr[0] = x1*x1 + x2*x2 - 13;// ||x||^2<=13

  static int count = 0;
  if (debug)
  {
    ++ count;
    printf("count=%d\n", count);
  }

  // check data is ok
  if (data != data_ref)
  {
    printf("invalid data\n");
    *f = NAN;
  }
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

  double x0[] = {0, 0};
  double xl[] = {-6.0, -6.0};
  double xu[] = {6.0, 6.0};
  prima_problem problem;
  prima_init_problem(&problem, n);
  problem.calcfc = &fun_con;
  problem.calfun = &fun;
  problem.x0 = x0;
  prima_options options;
  prima_init_options(&options);
  options.iprint = PRIMA_MSG_RHO;
  options.maxfun = 500*n;
  options.data = data_ref;
  problem.m_ineq = 3;
  double Aineq[3*2] = {1.0, 0.0,
                       0.0, 1.0,
                       1.0, 1.0};
  double bineq[3] = {4.0,
                     3.0,
                     10.0};
  problem.Aineq = Aineq;
  problem.bineq = bineq;
  problem.xl = xl;
  problem.xu = xu;
  prima_result result;
  int algorithm = 0;
  if(strcmp(algo, "bobyqa") == 0)
  {
    algorithm = PRIMA_BOBYQA;
  }
  else if(strcmp(algo, "cobyla") == 0)
  {
    algorithm = PRIMA_COBYLA;
    problem.m_nlcon = M_NLCON;
  }
  else if(strcmp(algo, "lincoa") == 0)
  {
    algorithm = PRIMA_LINCOA;
  }
  else if(strcmp(algo, "newuoa") == 0)
  {
    algorithm = PRIMA_NEWUOA;
  }
  else if(strcmp(algo, "uobyqa") == 0)
  {
    algorithm = PRIMA_UOBYQA;
  }
  else
  {
    printf("incorrect algo\n");
    return 1;
  }
  int rc = prima_minimize(algorithm, &problem, &options, &result);
  printf("f*=%g cstrv=%g nlconstr=%g rc=%d msg='%s' evals=%d\n", result.f, result.cstrv, result.nlconstr ? result.nlconstr[0] : 0.0, rc, result.message, result.nf);
  prima_free_problem(&problem);
  prima_free_result(&result);
  return (fabs(result.x[0]-3)>2e-2 || fabs(result.x[1]-2)>2e-2);
}
