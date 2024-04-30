// A stress test on excessively large problems

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Make PRIMA available
#include "prima/prima.h"

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

// Objective function for unconstrained, bound constrained, and linearly-constrained problems
static void fun(const double x[], double *const f, const void *data)
{
    // Objective: Rosenbrock function
    *f = 0.0;
    for (int i = 0; i < n-1; ++ i)
        *f += (x[i] - 1.0) * (x[i] - 1.0) + alpha * (x[i+1] - x[i]*x[i]) * (x[i+1] - x[i]*x[i]);

    static int nf = 0;
    if (debug) {
        ++ nf;
        printf("Number of function evaluations = %d\n", nf);
    }
    (void)data;
}

// Objective & constraint function for nonlinearly-constrained problems
static void fun_con(const double x[], double *const f, double constr[], const void *data)
{
    // Objective: Rosenbrock function
    *f = 0.0;
    for (int i = 0; i < n-1; ++ i)
        *f += (x[i] - 1.0) * (x[i] - 1.0) + alpha * (x[i+1] - x[i]*x[i]) * (x[i+1] - x[i]*x[i]);

    // Constraint: x_{i+1} <= x_i^2
    for (int i = 0; i < MIN(M_NLCON, n-1); ++ i)
        constr[i] = x[i+1] - x[i] * x[i];

    static int nf = 0;
    if (debug) {
        ++ nf;
        printf("Number of function evaluations = %d\n", nf);
    }
    (void)data;
}

// A function generating a seed that alters weekly
unsigned int get_random_seed(void)
{
    // Set the random seed to year/week
    char buf[10] = {0};
    time_t t = time(NULL);
    struct tm timeinfo;
    localtime_s(&timeinfo, &t);
    int rc = strftime(buf, 10, "%y%W", &timeinfo);
    if (!rc)
        return 42;
    else
        return atoi(buf);
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

    unsigned int seed = get_random_seed();
    printf("Random seed = %u\n", seed);
    srand(seed);

    // Set up the options
    prima_options_t options;
    prima_init_options(&options);
    options.iprint = PRIMA_MSG_RHO;
    options.maxfun = 500*N_MAX;

    // Data for the problem
    double x0[N_MAX];
    double xl[N_MAX];
    double xu[N_MAX];
    for (int i = 0; i < N_MAX; ++ i) {
        x0[i] = random_gen(-1.0, 1.0);
        xl[i] = random_gen(-2.0, -1.0);
        xu[i] = random_gen(1.0, 2.0);
    }
    double *Aineq = malloc(N_MAX*M_INEQ_MAX*sizeof(double));
    for (int i = 0; i < N_MAX; ++ i) {
        for (int j = 0; j < m_ineq; ++ j)
            Aineq[j*N_MAX+i] = random_gen(-1.0, 1.0);
    }
    double bineq[M_INEQ_MAX];
    for (int j = 0; j < M_INEQ_MAX; ++ j)
        bineq[j] = random_gen(-1.0, 1.0);

    // Define the algorithm and the problem according to `algo`
    prima_problem_t problem;
    prima_init_problem(&problem, N_MAX);
    problem.x0 = x0;
    if(strcmp(algo, "uobyqa") == 0) {
        algorithm = PRIMA_UOBYQA;
        problem.n = 100;
        problem.calfun = &fun;
    }
    else if(strcmp(algo, "newuoa") == 0) {
        algorithm = PRIMA_NEWUOA;
        problem.n = 1600;
        problem.calfun = &fun;
    }
    else if(strcmp(algo, "bobyqa") == 0) {
        algorithm = PRIMA_BOBYQA;
        problem.n = 1600;
        problem.calfun = &fun;
        problem.xl = xl;
        problem.xu = xu;
    }
    else if(strcmp(algo, "lincoa") == 0) {
        algorithm = PRIMA_LINCOA;
        problem.n = 1000;
        problem.m_ineq = 1000;
        problem.calfun = &fun;
        problem.xl = xl;
        problem.xu = xu;
        problem.Aineq = Aineq;
        problem.bineq = bineq;
    }
    else if(strcmp(algo, "cobyla") == 0) {
        algorithm = PRIMA_COBYLA;
        problem.n = 800;
        problem.m_nlcon = M_NLCON;
        problem.m_ineq = 600;
        problem.calcfc = &fun_con;
        problem.xl = xl;
        problem.xu = xu;
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

    // Free the result
    prima_free_result(&result);

    return 0;
}
