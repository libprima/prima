
/* Dedicated to the late Professor M. J. D. Powell FRS (1936--2015). */

#include "prima/prima.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXFUN_DIM_DFT 500

int prima_init_options(prima_options *options)
{
  if (options)
  {
    memset(options, 0, sizeof(prima_options));
    options->maxfun = -1;// interpreted as MAXFUN_DIM_DFT*n
    options->rhobeg = 1.0;
    options->rhoend = 1e-6;
    options->iprint = PRIMA_MSG_NONE;
    options->ftarget = -INFINITY;
    options->npt = -1;// interpreted as 2*n+1
    return 0;
  }
  else
    return PRIMA_NULL_OPTIONS;
}

int prima_init_problem(prima_problem *problem, int n)
{
  if (problem)
  {
    memset(problem, 0, sizeof(prima_problem));
    problem->n = n;
    problem->f0 = NAN;
    return 0;
  }
  else
    return PRIMA_NULL_PROBLEM;
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

int prima_check_problem(prima_problem *problem, prima_options *options, int alloc_bounds, int use_constr)
{
  if (!problem)
    return PRIMA_NULL_PROBLEM;
  if (!options)
    return PRIMA_NULL_OPTIONS;
  if (!problem->x0)
    return PRIMA_NULL_X0;
  if ((use_constr && !problem->calcfc) || (!use_constr && !problem->calfun))
    return PRIMA_NULL_FUNCTION;
  if (alloc_bounds && !problem->xl)
  {
    problem->xl = (double*)malloc(problem->n * sizeof(double));
    if (!problem->xl)
      return PRIMA_MEMORY_ALLOCATION_FAILS;
    problem->_allocated_xl = 1;
    for (int i = 0; i< problem->n; ++ i)
      problem->xl[i] = -INFINITY;
  }
  if (alloc_bounds && !problem->xu)
  {
    problem->xu = (double*)malloc(problem->n * sizeof(double));
    if (!problem->xu)
      return PRIMA_MEMORY_ALLOCATION_FAILS;
    problem->_allocated_xu = 1;
    for (int i = 0; i< problem->n; ++ i)
      problem->xu[i] = INFINITY;
  }
  if (options->maxfun < 0)
    options->maxfun = MAXFUN_DIM_DFT*problem->n;
  if (options->npt < 0)
    options->npt = 2*problem->n+1;
  return 0;
}

int prima_free_problem(prima_problem *problem)
{
  if (problem)
  {
    if (problem->_allocated_xl)
    {
      free(problem->xl);
      problem->xl = NULL;
      problem->_allocated_xl = 0;
    }
    if (problem->_allocated_xu)
    {
      free(problem->xu);
      problem->xu = NULL;
      problem->_allocated_xu = 0;
    }
    if (problem->_allocated_nlconstr0)
    {
      free(problem->nlconstr0);
      problem->nlconstr0 = NULL;
      problem->_allocated_nlconstr0 = 0;
    }
    return 0;
  }
  else
    return PRIMA_NULL_PROBLEM;
}

int prima_init_result(prima_result *result, prima_problem *problem)
{
  if (result)
  {
    memset(result, 0, sizeof(prima_result));
    result->f = 0.0;
    result->cstrv = 0.0;
    if (!problem)
      return PRIMA_NULL_PROBLEM;
    if (!problem->x0)
      return PRIMA_NULL_X0;
    result->x = (double*)malloc(problem->n * sizeof(double));
    if (!result->x)
      return PRIMA_MEMORY_ALLOCATION_FAILS;
    memcpy(result->x, problem->x0, problem->n * sizeof(double));
    return 0;
  }
  else
    return PRIMA_NULL_RESULT;
}

int prima_free_result(prima_result *result)
{
  if (result)
  {
    if (result->nlconstr)
    {
      free(result->nlconstr);
      result->nlconstr = NULL;
      result->_m_nlcon = 0;
    }
    return 0;
  }
  else
    return PRIMA_NULL_RESULT;
}

/* these functions just call the fortran compatibility layer and return the status code */
int prima_cobyla(prima_problem *problem, prima_options *options, prima_result *result)
{
  int info = prima_check_problem(problem, options, 1, 1);
  if (info == 0)
    info = prima_init_result(result, problem);
  if (info == 0)
  {
    // reuse or (re)allocate nlconstr
    if (result->_m_nlcon != problem->m_nlcon)
    {
      if (result->nlconstr)
        free(result->nlconstr);
      result->nlconstr = (double*)calloc(problem->m_nlcon, sizeof(double));
      if (!(result->nlconstr))
        info = PRIMA_MEMORY_ALLOCATION_FAILS;
      result->_m_nlcon = problem->m_nlcon;
    }
  }
  if (info == 0)
  {
    // evaluate f0, nlconstr0 if either one is not provided
    if (problem->f0 == NAN || !problem->nlconstr0)
    {
      if (!problem->nlconstr0)
      {
        problem->nlconstr0 = (double*)calloc(problem->n, sizeof(double));
        if (!problem->nlconstr0)
          return PRIMA_MEMORY_ALLOCATION_FAILS;
        problem->_allocated_nlconstr0 = 1;
      }
      problem->calcfc(result->x, &(problem->f0), problem->nlconstr0, options->data);
    }
  }
  if (info == 0)
  {
    cobyla_c(problem->m_nlcon, problem->calcfc, options->data, problem->n, result->x, &(result->f), &(result->cstrv), result->nlconstr,
            problem->m_ineq, problem->Aineq, problem->bineq, problem->m_eq, problem->Aeq, problem->beq,
            problem->xl, problem->xu, problem->f0, problem->nlconstr0, &(result->nf), options->rhobeg, options->rhoend, options->ftarget, options->maxfun, options->iprint, &info);
    result->status = info;
    result->message = prima_get_rc_string(info);
  }
  return info;
}

int prima_bobyqa(prima_problem *problem, prima_options *options, prima_result *result)
{
  int info = prima_check_problem(problem, options, 1, 0);
  if (info == 0)
    info = prima_init_result(result, problem);
  if (info == 0)
  {
    bobyqa_c(problem->calfun, options->data, problem->n, result->x, &(result->f), problem->xl, problem->xu, &(result->nf), options->rhobeg, options->rhoend, options->ftarget, options->maxfun, options->npt, options->iprint, &info);
    result->status = info;
    result->message = prima_get_rc_string(info);
  }
  return info;
}

int prima_newuoa(prima_problem *problem, prima_options *options, prima_result *result)
{
  int info = prima_check_problem(problem, options, 0, 0);
  if (info == 0)
    info = prima_init_result(result, problem);
  if (info == 0)
  {
    newuoa_c(problem->calfun, options->data, problem->n, result->x, &(result->f), &(result->nf), options->rhobeg, options->rhoend, options->ftarget, options->maxfun, options->npt, options->iprint, &info);
    result->status = info;
    result->message = prima_get_rc_string(info);
  }
  return info;
}

int prima_uobyqa(prima_problem *problem, prima_options *options, prima_result *result)
{
  int info = prima_check_problem(problem, options, 0, 0);
  if (info == 0)
    info = prima_init_result(result, problem);
  if (info == 0)
  {
    uobyqa_c(problem->calfun, options->data, problem->n, result->x, &(result->f), &(result->nf), options->rhobeg, options->rhoend, options->ftarget, options->maxfun, options->iprint, &info);
    result->status = info;
    result->message = prima_get_rc_string(info);
  }
  return info;
}

int prima_lincoa(prima_problem *problem, prima_options *options, prima_result *result)
{
  int info = prima_check_problem(problem, options, 1, 0);
  if (info == 0)
    info = prima_init_result(result, problem);
  if (info == 0)
  {
    lincoa_c(problem->calfun, options->data, problem->n, result->x, &(result->f), &(result->cstrv),
            problem->m_ineq, problem->Aineq, problem->bineq, problem->m_eq, problem->Aeq, problem->beq,
            problem->xl, problem->xu, &(result->nf), options->rhobeg, options->rhoend, options->ftarget, options->maxfun, options->npt, options->iprint, &info);
    result->status = info;
    result->message = prima_get_rc_string(info);
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
    case PRIMA_NULL_PROBLEM:
      return "NULL problem";
    case PRIMA_NULL_X0:
      return "NULL x0";
    case PRIMA_NULL_RESULT:
      return "NULL result";
    case PRIMA_NULL_FUNCTION:
      return "NULL function";
    default:
      return "Invalid return code";
  }
}
