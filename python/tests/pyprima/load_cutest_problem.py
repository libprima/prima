import numpy as np
from contextlib import redirect_stdout, redirect_stderr
import os


def load_cutest_problem(problem_name):
    """
    Load a CUTEst problem.

    Parameters
    ----------
    problem_name : str
        Name of the CUTEst problem to load.

    Returns
    -------
    fun : callable
        Objective function
    x0 : numpy ndarray
        Starting point of the problem
    lb : numpy ndarray of same shape as x0.
        Lower bounds. If unbounded, will contain -inf
    ub : numpy ndarray of same shape as x0.
        Upper bounds. If unbounded, will contain +inf
    constraints : dict or None
        Dictionary of constraints with the following keys:
          - 'a_ub' - numpy ndarray
          - 'b_ub' - numpy ndarray
          - 'a_eq' - numpy ndarray
          - 'b_eq' - numpy ndarray
          - 'm_nonlinear_ub' - int
          - 'c_ub' - callable
          - 'm_nonlinear_eq' - int
          - 'c_eq' - callable

        All of these keys will be present if a dictionary is returned.
        b_ub.size and b_eq.size can used to check if the relevant constraints
        are present. Similarly m_nonlinear_ub and m_nonlinear_eq can be used
        to check of c_ub and c_eq should be used, respectively.

        If None is returned there are no constraints.
          


    Raises
    ------
    Exception
        If there was a failure in loading the cutest problem
    
    """
    import pycutest
    
    def _build_linear_ub(cutest_problem):
        """
        Build the linear inequality constraints from a CUTEst problem.
        """
        idx_ub = cutest_problem.is_linear_cons & ~cutest_problem.is_eq_cons
        idx_ub_cl = cutest_problem.cl[idx_ub] > -1e20
        idx_ub_cu = cutest_problem.cu[idx_ub] < 1e20
        a_ub = []
        b_ub = []
        for i, index in enumerate(np.flatnonzero(idx_ub)):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            if idx_ub_cl[i]:
                a_ub.append(-g_val)
                b_ub.append(c_val - cutest_problem.cl[index])
            if idx_ub_cu[i]:
                a_ub.append(g_val)
                b_ub.append(cutest_problem.cu[index] - c_val)
        return np.reshape(a_ub, (-1, cutest_problem.n)), np.array(b_ub)
    
    def _build_linear_eq(cutest_problem):
        """
        Build the linear equality constraints from a CUTEst problem.
        """
        idx_eq = cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
        a_eq = []
        b_eq = []
        for index in np.flatnonzero(idx_eq):
            c_val, g_val = cutest_problem.cons(np.zeros(cutest_problem.n), index, True)
            a_eq.append(g_val)
            b_eq.append(0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]) - c_val)
        return np.reshape(a_eq, (-1, cutest_problem.n)), np.array(b_eq)
    
    def _build_c_ub(cutest_problem):
        """
        Build the nonlinear inequality constraints of a CUTEst problem.
        """
        def c_ub(x):
            idx_ub = ~(cutest_problem.is_linear_cons | cutest_problem.is_eq_cons)
            idx_ub_cl = cutest_problem.cl[idx_ub] > -1e20
            idx_ub_cu = cutest_problem.cu[idx_ub] < 1e20
            c = []
            for i, index in enumerate(np.flatnonzero(idx_ub)):
                c_val = cutest_problem.cons(x, index)
                if idx_ub_cl[i]:
                    c.append(cutest_problem.cl[index] - c_val)
                if idx_ub_cu[i]:
                    c.append(c_val - cutest_problem.cu[index])
            return np.array(c)
        return c_ub
    
    def _build_c_eq(cutest_problem):
        """
        Build the nonlinear equality constraints of a CUTEst problem.
        """
        def c_eq(x):
            idx_eq = ~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons
            c = []
            for index in np.flatnonzero(idx_eq):
                c_val = cutest_problem.cons(x, index)
                c.append(0.5 * (cutest_problem.cl[index] + cutest_problem.cu[index]) - c_val)
            return np.array(c)
        return c_eq
    
    # Preprocess the problem name.
    if not isinstance(problem_name, str):
        raise TypeError('The argument problem_name must be a string.')

    # Attempt to load the CUTEst problem.
    cutest_problem = None
    print(f'Loading CUTEst problem {problem_name}.')
    try:
        with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                cutest_problem = pycutest.import_problem(problem_name)
    except Exception as err:
        print(f'Failed to load CUTEst problem {problem_name}: {err}')

    # If the problem is not successfully loaded, raise an exception.
    if cutest_problem is None:
        raise Exception(f'Failed to load CUTEst problem {problem_name}.')
    
    # The problem is successfully loaded.
    # Cutest uses 1e20 to signify infinity/unbounded
    CUTEST_INF = 1e20
    lb = np.array(cutest_problem.bl)
    lb[lb <= -CUTEST_INF] = -np.inf
    ub = np.array(cutest_problem.bu)
    ub[ub >= CUTEST_INF] = np.inf
    if cutest_problem.m > 0:
        constraints = {'c_ub': _build_c_ub(cutest_problem), 'c_eq': _build_c_eq(cutest_problem)}
        constraints['a_ub'], constraints['b_ub'] = _build_linear_ub(cutest_problem)
        constraints['a_eq'], constraints['b_eq'] = _build_linear_eq(cutest_problem)
        idx_nonlinear_ub = ~(cutest_problem.is_linear_cons | cutest_problem.is_eq_cons)
        constraints['m_nonlinear_ub'] = np.count_nonzero(cutest_problem.cl[idx_nonlinear_ub] > -1e20) + np.count_nonzero(cutest_problem.cu[idx_nonlinear_ub] < 1e20)
        constraints['m_nonlinear_eq'] = np.count_nonzero(~cutest_problem.is_linear_cons & cutest_problem.is_eq_cons)
    else:
        constraints = None
    print(f'CUTEst problem {cutest_problem.name} (n = {cutest_problem.n}) successfully loaded.')
    return cutest_problem.obj, cutest_problem.x0, cutest_problem.bl, cutest_problem.bu, constraints