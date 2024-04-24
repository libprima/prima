#ifndef PRIMA_INTERNAL_H
#define PRIMA_INTERNAL_H

#include "prima/prima.h"

unsigned int get_random_seed(void);

// Function to check whether the problem matches the algorithm
prima_rc_t prima_check_problem(const prima_problem_t problem, const prima_algorithm_t algorithm);

// Function to initialize the result
prima_rc_t prima_init_result(prima_result_t *const result, const prima_problem_t problem);

// Functions implemented in Fortran (*_c.f90)
int cobyla_c(const int m_nlcon, const prima_objcon_t calcfc, const void *data, const int n, double x[], double *const f, double *const cstrv, double nlconstr[],
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
             const double f0, const double nlconstr0[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint, const double ctol,
             const prima_callback_t callback, int *const info);

int bobyqa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f, const double xl[], const double xu[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, const prima_callback_t callback, int *const info);

int newuoa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f,
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, const prima_callback_t callback, int *const info);

int uobyqa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f,
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int iprint, const prima_callback_t callback, int *const info);

int lincoa_c(prima_obj_t calfun, const void *data, const int n, double x[], double *const f,
             double *const cstrv, const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[], const double xl[], const double xu[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget, const int maxfun, const int npt, const int iprint, const double ctol,
             const prima_callback_t callback, int *const info);

#endif
