// A test for the `data` argument in the objective/constraint callback function.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Make PRIMA available
#include "prima/prima.h"

#define M_NLCON 1

const int n = 2;
int debug = 0;
static int int_data = 0xff;
void * data_ref = &int_data;

// Objective function for unconstrained, bound constrained, and linearly-constrained problems
static void fun(const double x[], double *const f, const void *data)
{
    const double x1 = x[0];
    const double x2 = x[1];
    *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2)+0.1*(x1+x2)-10;

    static int nf = 0;
    if (debug) {
        ++ nf;
        printf("Number of function evaluations = %d\n", nf);
    }

    // Check whether data is OK
    if (data != data_ref) {
        printf("Invalid data!\n");
        *f = NAN;
    }
}

// Objective & constraint function for nonlinearly-constrained problems
static void fun_con(const double x[], double *const f, double constr[], const void *data)
{
    const double x1 = x[0];
    const double x2 = x[1];
    *f = 5*(x1-3)*(x1-3)+7*(x2-2)*(x2-2) + 0.1*(x1+x2) - 10;
    constr[0] = x1*x1 + x2*x2 - 13;  // ||x||^2<=13

    static int nf = 0;
    if (debug) {
        ++ nf;
        printf("Number of function evaluations = %d\n", nf);
    }

    // Check whether data is OK
    if (data != data_ref) {
        printf("Invalid data!\n");
        *f = NAN;
    }
}

// Main function
int main(int argc, char * argv[])
{
    char *algo = "uobyqa";
    prima_algorithm_t algorithm = PRIMA_UOBYQA;
    if (argc > 1)
        algo = argv[1];
    printf("Algorithm = %s\n", algo);

    if (argc > 2)
        debug = (strcmp(argv[2], "debug") == 0);
    printf("Debug = %d\n", debug);

    // Set up the options
    prima_options_t options;
    prima_init_options(&options);
    options.iprint = PRIMA_MSG_RHO;
    options.maxfun = 500*n;
    options.data = data_ref;

    // Data for the problem
    double x0[] = {0.0,
        0.0};
    double xl[] = {-6.0,
        -6.0};
    double xu[] = {6.0,
        6.0};
    double Aineq[3*2] = {1.0, 0.0,
        0.0, 1.0,
        1.0, 1.0};
    double bineq[3] = {4.0,
        3.0,
        10.0};

    // Define the algorithm and the problem according to `algo`
    prima_problem_t problem;
    prima_init_problem(&problem, n);
    problem.x0 = x0;
    if(strcmp(algo, "uobyqa") == 0) {
        algorithm = PRIMA_UOBYQA;
        problem.calfun = &fun;
    }
    else if(strcmp(algo, "newuoa") == 0) {
        algorithm = PRIMA_NEWUOA;
        problem.calfun = &fun;
    }
    else if(strcmp(algo, "bobyqa") == 0) {
        algorithm = PRIMA_BOBYQA;
        problem.calfun = &fun;
        problem.xl = xl;
        problem.xu = xu;
    }
    else if(strcmp(algo, "lincoa") == 0) {
        algorithm = PRIMA_LINCOA;
        problem.calfun = &fun;
        problem.xl = xl;
        problem.xu = xu;
        problem.m_ineq = 3;
        problem.Aineq = Aineq;
        problem.bineq = bineq;
    }
    else if(strcmp(algo, "cobyla") == 0) {
        algorithm = PRIMA_COBYLA;
        problem.m_nlcon = M_NLCON;
        problem.calcfc = &fun_con;
        problem.xl = xl;
        problem.xu = xu;
        problem.m_ineq = 3;
        problem.Aineq = Aineq;
        problem.bineq = bineq;
    }
    else {
        printf("Invalid algorithm %s!\n", algo);
        return 1;
    }

    // Call the solver
    prima_result_t result;
    const prima_rc_t rc = prima_minimize(algorithm, problem, options, &result);

    // Print the result
    printf("f* = %g, cstrv = %g, nlconstr = {%g}, rc = %d, msg = '%s', evals = %d\n", result.f, result.cstrv, result.nlconstr ? result.nlconstr[0] : 0.0, rc, result.message, result.nf);

    // Check the result
    int success = (fabs(result.x[0] - 3) > 2e-2 || fabs(result.x[1] - 2) > 2e-2);

    // Free the result
    prima_free_result(&result);

    return success;
}
