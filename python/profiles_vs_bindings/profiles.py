import numpy as np
import sys
import os
from excludelist import excludelist
from optiprofiler import benchmark
import argparse
from time import time
import prima.backends.pyprima.common.linalg
prima.backends.pyprima.common.linalg.USE_NAIVE_MATH = True
from prima import minimize, Bounds, LinearConstraint, NonlinearConstraint
from prima.backends.pybindings import __cmake_build_type__
assert __cmake_build_type__ == 'Debug', "Fortran must be in debug mode to compare against the Python implementation"



nondefault_options = lambda n, f0: {
    'ftarget' : f0 - 314, # if this is doable, 3.14 otherwise
    'maxfev' : 271*n,
    'npt' : int(min(3.14*n, n**1.23)),
    'rhobeg' : 2.71828,
    'rhoend' : 3.14159*1.0e-4,
}

def get_pyprima_options(n, f0):
    options = {'backend': 'Python'}
    if os.environ.get('NONDEFAULT_PYPRIMA') == 'True':
        additional_options = nondefault_options(n, f0)
        # Change the option name
        additional_options['maxfun'] = additional_options.pop('maxfev')
        options |= additional_options
    return options


def get_bindings_options(n, f0):
    options = {'backend': 'Fortran'}
    if os.environ.get('NONDEFAULT_BINDINGS') == 'True':
        additional_options = nondefault_options(n, f0)
        options |= additional_options
    return options


def pyprima_cobyla(fun, x0, lb=None, ub=None, a_ub=None, b_ub=None, a_eq=None, b_eq=None, c_ub=None, c_eq=None):
    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub) if lb is not None and ub is not None else None
    constraints = []
    if b_ub is not None and b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq is not None and b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    if c_ub is not None:
        c_ub_x0 = c_ub(x0)
        if c_ub_x0.size > 0:
            constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    if c_eq is not None:
        c_eq_x0 = c_eq(x0)
        if c_eq_x0.size > 0:
            constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    options = get_pyprima_options(len(x0), f0)
    if 'npt' in options:
        del options['npt']
    result = minimize(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options=options)
    return result.x


def bindings_cobyla(fun, x0, lb=None, ub=None, a_ub=None, b_ub=None, a_eq=None, b_eq=None, c_ub=None, c_eq=None):
    f0 = fun.__self__._fun(x0)
    bounds = Bounds(lb, ub) if lb is not None and ub is not None else None
    constraints = []
    if b_ub is not None and b_ub.size > 0:
        constraints.append(LinearConstraint(a_ub, -np.inf, b_ub))
    if b_eq is not None and b_eq.size > 0:
        constraints.append(LinearConstraint(a_eq, b_eq, b_eq))
    if c_ub is not None:
        c_ub_x0 = c_ub(x0)
        if c_ub_x0.size > 0:
            constraints.append(NonlinearConstraint(c_ub, -np.inf, np.zeros_like(c_ub_x0)))
    if c_eq is not None:
        c_eq_x0 = c_eq(x0)
        if c_eq_x0.size > 0:
            constraints.append(NonlinearConstraint(c_eq, np.zeros_like(c_eq_x0), np.zeros_like(c_eq_x0)))
    res = minimize(fun, x0, method='cobyla', bounds=bounds, constraints=constraints, options=get_bindings_options(len(x0), f0))
    return res.x


if __name__ == '__main__':
    # If we run this script from a directory other than the one that contains it, pycutest's call to importlib will fail,
    # unless we insert the current working directory into the path.
    sys.path.insert(0, os.getcwd())
    os.environ['PYCUTEST_CACHE'] = os.getcwd()
    
    parser = argparse.ArgumentParser(description='Generate performance profiles comparing PyPRIMA to PRIMA Python (bindings).')
    parser.add_argument('-j', '--n_jobs', type=int, default=1, help='Number of jobs to run in parallel')
    parser.add_argument('--default_only', action='store_true', help='Run only the default options for both PyPRIMA and PRIMA')
    args = parser.parse_args()

    

    def run_three_benchmarks(pyprima_fun, bindings_fun, algorithm, cutest_problem_names, default_only, n_jobs):
        '''
        Proper validation of both default and nondefault options requires 3 runs: Both default, both nondefault, and
        one default one nondefault. The first two should look identical, and so the third run confirms that our
        experiment setup is valid (i.e. it rules out a scenario where even though options are provided, both algorithms
        end up using default options anyway).
        '''
        # Sharing state with multiprocessing is hard when we can't control the function signature,
        # so we resort to using the environment to pass options.
        algorithm = algorithm.lower()
        ALGORITHM = algorithm.upper()
        project_x0 = algorithm == 'lincoa'
        os.environ['NONDEFAULT_PYPRIMA'] = "False"
        os.environ['NONDEFAULT_BINDINGS'] = "False"
        benchmark([pyprima_fun, bindings_fun], plibs='pycutest', solver_names=[f'PyPRIMA-{ALGORITHM}', f'Bindings-{ALGORITHM}'], problem_names=cutest_problem_names, benchmark_id=f'{algorithm}_default_options', n_jobs=n_jobs, project_x0=project_x0)
        if not default_only:
            os.environ['NONDEFAULT_PYPRIMA'] = "True"
            os.environ['NONDEFAULT_BINDINGS'] = "True"
            benchmark([pyprima_fun, bindings_fun], plibs='pycutest', solver_names=[f'PyPRIMA-{ALGORITHM}', f'Bindings-{ALGORITHM}'], problem_names=cutest_problem_names, benchmark_id=f'{algorithm}_nondefault_options', n_jobs=n_jobs, project_x0=project_x0)
            os.environ['NONDEFAULT_PYPRIMA'] = "True"
            os.environ['NONDEFAULT_BINDINGS'] = "False"
            benchmark([pyprima_fun, bindings_fun], plibs='pycutest', solver_names=[f'PyPRIMA-{ALGORITHM}', f'Bindings-{ALGORITHM}'], problem_names=cutest_problem_names, benchmark_id=f'{algorithm}_different_options', n_jobs=n_jobs, project_x0=project_x0)

    start = time()
    print("Running profiles for COBYLA")
    with open('cobyla.txt') as f:
        cutest_problem_names = f.read().splitlines()
    cutest_problem_names = list(filter(lambda x: x not in excludelist('cobyla'), cutest_problem_names))
    run_three_benchmarks(pyprima_cobyla, bindings_cobyla, 'cobyla', cutest_problem_names, args.default_only, args.n_jobs)
    print(f'Completed COBYLA profile in {time() - start:.2f} seconds')
