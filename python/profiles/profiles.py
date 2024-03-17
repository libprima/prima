import numpy as np
import sys
import os
from excludelist import excludelist
from optiprofiler import set_cutest_problem_options, find_cutest_problems, run_benchmark
import argparse
from time import time


def pdfo_uobyqa(fun, x0):
    from pdfo import pdfo

    res = pdfo(fun, x0, method='uobyqa')
    return res.x


def prima_uobyqa(fun, x0):
    from prima import minimize
    
    res = minimize(fun, x0, method='uobyqa')
    return res.x


def pdfo_newuoa(fun, x0):
    from pdfo import pdfo

    res = pdfo(fun, x0, method='newuoa')
    return res.x


def prima_newuoa(fun, x0):
    from prima import minimize

    res = minimize(fun, x0, method='newuoa')
    return res.x


def pdfo_bobyqa(fun, x0, lb, ub):
    from pdfo import pdfo
    from scipy.optimize import Bounds

    bounds = Bounds(lb, ub)
    res = pdfo(fun, x0, method='bobyqa', bounds=bounds)
    return res.x


def prima_bobyqa(fun, x0, lb, ub):
    from prima import minimize
    from scipy.optimize import Bounds

    bounds = Bounds(lb, ub)
    res = minimize(fun, x0, method='bobyqa', bounds=bounds)
    return res.x


def pdfo_lincoa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint

    # Until we include the elimination of linear equality constraints in PRIMA,
    # we disable it in PDFO to make the comparison fair.
    options = {'eliminate_lin_eq': False}
    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    res = pdfo(fun, x0, method='lincoa', bounds=bounds, constraints=constraints, options=options)
    return res.x


def prima_lincoa(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq):
    from prima import minimize, Bounds, LinearConstraint

    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    res = minimize(fun, x0, method='lincoa', bounds=bounds, constraints=constraints)
    return res.x


def pdfo_cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq):
    from pdfo import pdfo
    from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint

    # Until we include the elimination of linear equality constraints in PRIMA,
    # we disable it in PDFO to make the comparison fair.
    options = {'eliminate_lin_eq': False}
    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ub_x0 = c_ub(x0)
    if c_ub_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    res = pdfo(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options=options)
    return res.x


def prima_cobyla(fun, x0, lb, ub, a_ub, b_ub, a_eq, b_eq, c_ub, c_eq):
    from prima import minimize, Bounds, LinearConstraint, NonlinearConstraint

    bounds = Bounds(lb, ub)
    constraints = []
    if b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    c_ub_x0 = c_ub(x0)
    if c_ub_x0.size > 0:
        constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    c_eq_x0 = c_eq(x0)
    if c_eq_x0.size > 0:
        constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    res = minimize(fun, x0, method='cobyla', bounds=bounds, constraints=constraints)
    return res.x


def get_problems(description):
    cutest_problem_names = find_cutest_problems(description)
    return list(filter(lambda x: x not in excludelist(description), cutest_problem_names))


if __name__ == '__main__':
    # If we run this script from a directory other than the one that contains it, pycutest's call to importlib will fail,
    # unless we insert the current working directory into the path.
    sys.path.insert(0, os.getcwd())
    
    parser = argparse.ArgumentParser(description='Generate performance profiles comparing PRIMA to PDFO. '
                                     'By default no profiles are run. To run one or more profiles use the flags '
                                     'specified below. Multiple flags may be specified to run multiple profiles, '
                                     'but please note they will be run sequentially.')
    parser.add_argument('-j', '--n_jobs', type=int, default=None, help='Number of jobs to run in parallel')
    parser.add_argument('-n', '--newuoa', action='store_true', help='Run the NEWUOA benchmark')
    parser.add_argument('-u', '--uobyqa', action='store_true', help='Run the UOBYQA benchmark')
    parser.add_argument('-b', '--bobyqa', action='store_true', help='Run the BOBYQA benchmark')
    parser.add_argument('-l', '--lincoa', action='store_true', help='Run the LINCOA benchmark')
    parser.add_argument('-c', '--cobyla', action='store_true', help='Run the COBYLA benchmark')
    args = parser.parse_args()
    
    os.environ['PYCUTEST_CACHE'] = os.getcwd()
    
    if args.newuoa:
        start = time()
        print("Running profiles for NEWUOA")
        set_cutest_problem_options(n_min=1, n_max=200, m_min = 0, m_max=0)
        cutest_problem_names = get_problems('unconstrained')
        run_benchmark([pdfo_newuoa, prima_newuoa], ['PDFO-NEWUOA', 'PRIMA-NEWUOA'], cutest_problem_names, benchmark_id='newuoa', n_jobs=args.n_jobs)
        print(f'Completed NEWUOA profile in {time() - start:.2f} seconds')
    
    if args.uobyqa:
        start = time()
        print("Running profiles for UOBYQA")
        set_cutest_problem_options(n_min=1, n_max=100, m_min = 0, m_max=0)
        cutest_problem_names = get_problems('unconstrained')
        run_benchmark([pdfo_uobyqa, prima_uobyqa], ['PDFO-UOBYQA', 'PRIMA-UOBYQA'], cutest_problem_names, benchmark_id='uobyqa', n_jobs=args.n_jobs)
        print(f'Completed UOBYQA profile in {time() - start:.2f} seconds')
    
    if args.bobyqa:
        start = time()
        print("Running profiles for BOBYQA")
        set_cutest_problem_options(n_min=1, n_max=200, m_min = 0, m_max=0)
        cutest_problem_names = get_problems('bound')
        run_benchmark([pdfo_bobyqa, prima_bobyqa], ['PDFO-BOBYQA', 'PRIMA-BOBYQA'], cutest_problem_names, benchmark_id='bobyqa', n_jobs=args.n_jobs)
        print(f'Completed BOBYQA profile in {time() - start:.2f} seconds')
    
    if args.lincoa:
        start = time()
        print("Running profiles for LINCOA")
        set_cutest_problem_options(n_min=1, n_max=200, m_min=1, m_max=20_000)
        cutest_problem_names = get_problems('adjacency linear')
        run_benchmark([pdfo_lincoa, prima_lincoa], ['PDFO-LINCOA', 'PRIMA-LINCOA'], cutest_problem_names, benchmark_id='lincoa', n_jobs=args.n_jobs)
        print(f'Completed LINCOA profile in {time() - start:.2f} seconds')
    
    if args.cobyla:
        start = time()
        print("Running profiles for COBYLA")
        set_cutest_problem_options(n_min=1, n_max=100, m_min=1, m_max=20_000)
        cutest_problem_names = get_problems('quadratic other')
        run_benchmark([pdfo_cobyla, prima_cobyla], ['PDFO-COBYLA', 'PRIMA-COBYLA'], cutest_problem_names, benchmark_id='cobyla', n_jobs=args.n_jobs)
        print(f'Completed COBYLA profile in {time() - start:.2f} seconds')

