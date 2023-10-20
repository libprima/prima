
/* Dedicated to the late Professor M. J. D. Powell FRS (1936--2015). */

#include "prima/prima.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXFUN_DIM_DFT 500

int prima_init_options(prima_options *opt)
{
  if (opt)
  {
    memset(opt, 0, sizeof(prima_options));
    opt->maxfun = -1;// interpreted as MAXFUN_DIM_DFT*n
    opt->rhobeg = 1.0;
    opt->rhoend = 1e-6;
    opt->iprint = PRIMA_MSG_NONE;
    opt->ftarget = -INFINITY;
    opt->npt = -1;// interpreted as 2*n+1
    opt->f0 = NAN;
    return 0;
  }
  else
    return PRIMA_NULL_OPTIONS;
}

/* implemented in fortran (*_c.f90) */
int cobyla_c(const int m_nlcon, const prima_objcon calcfc, const void *data, const int n, double x[], double *f, double *cstrv, double nlconstr[],
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
             double f0, const double nlconstr0[],
             int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint, int *info);
int bobyqa_c(prima_obj calfun, const void *data, const int n, double x[], double *f, const double xl[], const double xu[],
             int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, int *info);
int newuoa_c(prima_obj calfun, const void *data, const int n, double x[], double *f,
             int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, int *info);
int uobyqa_c(prima_obj calfun, const void *data, const int n, double x[], double *f,
             int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint, int *info);
int lincoa_c(prima_obj calfun, const void *data, const int n, double x[], double *f,
             double *cstrv,
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
             int *nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, int *info);

int prima_check_options(prima_options *opt, int n, int alloc_bounds)
{
  if (!opt)
    return PRIMA_NULL_OPTIONS;
  if (alloc_bounds && !opt->xl)
  {
    opt->xl = (double*)malloc(n * sizeof(double));
    if (!opt->xl)
      return PRIMA_MEMORY_ALLOCATION_FAILS;
    opt->_allocated_xl = 1;
    for (int i = 0; i< n; ++ i)
      opt->xl[i] = -INFINITY;
  }
  if (alloc_bounds && !opt->xu)
  {
    opt->xu = (double*)malloc(n * sizeof(double));
    if (!opt->xu)
      return PRIMA_MEMORY_ALLOCATION_FAILS;
    opt->_allocated_xu = 1;
    for (int i = 0; i< n; ++ i)
      opt->xu[i] = INFINITY;
  }
  if (opt->maxfun < 0)
    opt->maxfun = MAXFUN_DIM_DFT*n;
  if (opt->npt < 0)
    opt->npt = 2*n+1;
  return 0;
}

int prima_free_options(prima_options *opt)
{
  if (opt)
  {
    if (opt->_allocated_xl)
    {
      free(opt->xl);
      opt->xl = NULL;
      opt->_allocated_xl = 0;
    }
    if (opt->_allocated_xu)
    {
      free(opt->xu);
      opt->xu = NULL;
      opt->_allocated_xu = 0;
    }
    if (opt->_allocated_nlconstr0)
    {
      free(opt->nlconstr0);
      opt->nlconstr0 = NULL;
      opt->_allocated_nlconstr0 = 0;
    }
    return 0;
  }
  else
    return PRIMA_NULL_OPTIONS;
}

int prima_init_result(prima_results *results)
{
  if (results)
  {
    memset(results, 0, sizeof(prima_results));
    results->f = 0.0;
    results->cstrv = 0.0;
    return 0;
  }
  else
    return PRIMA_NULL_RESULT;
}

int prima_free_results(prima_results *results)
{
  if (results)
  {
    if (results->nlconstr)
    {
      free(results->nlconstr);
      results->nlconstr = NULL;
      results->_m_nlcon = 0;
    }
    return 0;
  }
  else
    return PRIMA_NULL_RESULT;
}

/* these functions just call the fortran compatibility layer and return the status code */
int prima_cobyla(const prima_objcon calcfc, const int n, double x[], prima_options *opt, prima_results *results)
{
  int info = prima_check_options(opt, n, 1);
  if (info == 0)
    info = prima_init_result(results);
  if (info == 0)
  {
    // reuse or (re)allocate nlconstr
    if (results->_m_nlcon != opt->m_nlcon)
    {
      if (results->nlconstr)
        free(results->nlconstr);
      results->nlconstr = (double*)calloc(opt->m_nlcon, sizeof(double));
      if (!(results->nlconstr))
        info = PRIMA_MEMORY_ALLOCATION_FAILS;
      results->_m_nlcon = opt->m_nlcon;
    }
  }
  if (info == 0)
  {
    // evaluate f0, nlconstr0 if either one is not provided
    if (opt->f0 == NAN || !opt->nlconstr0)
    {
      if (!opt->nlconstr0)
      {
        opt->nlconstr0 = (double*)calloc(n, sizeof(double));
        if (!opt->nlconstr0)
          return PRIMA_MEMORY_ALLOCATION_FAILS;
        opt->_allocated_nlconstr0 = 1;
      }
      calcfc(x, &(opt->f0), opt->nlconstr0, opt->data);
    }
  }
  if (info == 0)
  {
    cobyla_c(opt->m_nlcon, calcfc, opt->data, n, x, &(results->f), &(results->cstrv), results->nlconstr,
            opt->m_ineq, opt->Aineq, opt->bineq, opt->m_eq, opt->Aeq, opt->beq,
            opt->xl, opt->xu, opt->f0, opt->nlconstr0, &(results->nf), opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->iprint, &info);
  }
  return info;
}

int prima_bobyqa(const prima_obj calfun, const int n, double x[], prima_options *opt, prima_results *results)
{
  int info = prima_check_options(opt, n, 1);
  if (info == 0)
    info = prima_init_result(results);
  if (info == 0)
  {
    bobyqa_c(calfun, opt->data, n, x, &(results->f), opt->xl, opt->xu, &(results->nf), opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->npt, opt->iprint, &info);
  }
  return info;
}

int prima_newuoa(const prima_obj calfun, const int n, double x[], prima_options *opt, prima_results *results)
{
  int info = prima_check_options(opt, n, 0);
  if (info == 0)
    info = prima_init_result(results);
  if (info == 0)
  {
    newuoa_c(calfun, opt->data, n, x, &(results->f), &(results->nf), opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->npt, opt->iprint, &info);
  }
  return info;
}

int prima_uobyqa(const prima_obj calfun, const int n, double x[], prima_options *opt, prima_results *results)
{
  int info = prima_check_options(opt, n, 0);
  if (info == 0)
    info = prima_init_result(results);
  if (info == 0)
  {
    uobyqa_c(calfun, opt->data, n, x, &(results->f), &(results->nf), opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->iprint, &info);
  }
  return info;
}

int prima_lincoa(const prima_obj calfun, const int n, double x[], prima_options *opt, prima_results *results)
{
  int info = prima_check_options(opt, n, 1);
  if (info == 0)
    info = prima_init_result(results);
  if (info == 0)
  {
    lincoa_c(calfun, opt->data, n, x, &(results->f), &(results->cstrv),
            opt->m_ineq, opt->Aineq, opt->bineq, opt->m_eq, opt->Aeq, opt->beq,
            opt->xl, opt->xu, &(results->nf), opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->npt, opt->iprint, &info);
  }
  return info;
}

const char *prima_get_rc_string(const int rc)
{
  switch (rc)
  {
    case PRIMA_SMALL_TR_RADIUS:
      return "Trust region radius reaches its lower bound";
    case PRIMA_FTARGET_ACHIEVED:
      return "The target function value is reached";
    case PRIMA_TRSUBP_FAILED:
      return "A trust region step failed to reduce the model";
    case PRIMA_MAXFUN_REACHED:
      return "Maximum number of function evaluations reached";
    case PRIMA_MAXTR_REACHED:
      return "Maximum number of trust region iterations reached";
    case PRIMA_NAN_INF_X:
      return "The input X contains NaN of Inf";
    case PRIMA_NAN_INF_F:
      return "The objective or constraint functions return NaN or +Inf";
    case PRIMA_NAN_INF_MODEL:
      return "NaN or Inf occurs in the model";
    case PRIMA_NO_SPACE_BETWEEN_BOUNDS:
      return "No space between bounds";
    case PRIMA_DAMAGING_ROUNDING:
      return "Rounding errors are becoming damaging";
    case PRIMA_ZERO_LINEAR_CONSTRAINT:
      return "One of the linear constraints has a zero gradient";
    case PRIMA_INVALID_INPUT:
      return "Invalid input";
    case PRIMA_ASSERTION_FAILS:
      return "Assertion fails";
    case PRIMA_VALIDATION_FAILS:
      return "Validation fails";
    case PRIMA_MEMORY_ALLOCATION_FAILS:
      return "Memory allocation failed";
    case PRIMA_NULL_OPTIONS:
      return "NULL options";
    case PRIMA_NULL_RESULT:
      return "NULL result";
    default:
      return "Invalid return code";
  }
}
