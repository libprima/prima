
/* Dedicated to the late Professor M. J. D. Powell FRS (1936--2015). */

#include "prima/prima.h"
#include "stddef.h"
#include "math.h"

int prima_init_options(prima_options * opt)
{
  if (opt != NULL)
  {
    opt->maxfun = 100;
    opt->rhobeg = 1.0;
    opt->rhoend = 1e-6;
    opt->iprint = PRIMA_MSG_NONE;
    opt->ftarget = -INFINITY;
    opt->npt = -1;// interpreted as 2*n+1
    opt->data = NULL;
    opt->m_ineq = 0;
    opt->Aineq = NULL;
    opt->bineq = NULL;
    opt->m_eq = 0;
    opt->Aeq = NULL;
    opt->beq = NULL;
    opt->m_nlcon = 0;
    return 0;
  }
  return 1;
}

/* implemented in fortran (*_c.f90) */
int cobyla_c(const int m_nlcon, const prima_objcon calcfc, const void *data, const int n, double x[], double *f, double *cstrv, double nlconstr[],
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
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

/* these functions just call the fortran compatibility layer and return the status code */
int prima_cobyla(const prima_objcon calcfc, const int n, double x[], double *f,
                 double *cstrv, double nlconstr[],
                 const double xl[], const double xu[],
                 int *nf, const prima_options * opt)
{
  if (!opt)
    return PRIMA_NULL_OPTIONS;
  int info = 0;
  cobyla_c(opt->m_nlcon, calcfc, opt->data, n, x, f, cstrv, nlconstr,
           opt->m_ineq, opt->Aineq, opt->bineq, opt->m_eq, opt->Aeq, opt->beq,
           xl, xu, nf, opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->iprint, &info);
  return info;
}

int prima_bobyqa(const prima_obj calfun, const int n, double x[], double *f, const double xl[], const double xu[],
                 int *nf, const prima_options * opt)
{
  if (!opt)
    return PRIMA_NULL_OPTIONS;
  int info = 0;
  bobyqa_c(calfun, opt->data, n, x, f, xl, xu, nf, opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->npt < 0 ? 2*n+1 : opt->npt, opt->iprint, &info);
  return info;
}

int prima_newuoa(const prima_obj calfun, const int n, double x[], double *f,
                 int *nf, const prima_options * opt)
{
  if (!opt)
    return PRIMA_NULL_OPTIONS;
  int info = 0;
  newuoa_c(calfun, opt->data, n, x, f, nf, opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->npt < 0 ? 2*n+1 : opt->npt, opt->iprint, &info);
  return info;
}

int prima_uobyqa(const prima_obj calfun, const int n, double x[], double *f,
                 int *nf, const prima_options * opt)
{
  if (!opt)
    return PRIMA_NULL_OPTIONS;
  int info = 0;
  uobyqa_c(calfun, opt->data, n, x, f, nf, opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->iprint, &info);
  return info;
}

int prima_lincoa(const prima_obj calfun, const int n, double x[], double *f, double *cstrv,
                 const double xl[], const double xu[],
                 int *nf, const prima_options * opt)
{
  if (!opt)
    return PRIMA_NULL_OPTIONS;
  int info = 0;
  lincoa_c(calfun, opt->data, n, x, f, cstrv,
           opt->m_ineq, opt->Aineq, opt->bineq, opt->m_eq, opt->Aeq, opt->beq,
           xl, xu, nf, opt->rhobeg, opt->rhoend, opt->ftarget, opt->maxfun, opt->npt < 0 ? 2*n+1 : opt->npt, opt->iprint, &info);
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
    default:
      return "Invalid return code";
  }
}
