# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

import platform
import sys
import warnings
from inspect import stack

import numpy as np

python_version = sys.version_info.major
if python_version >= 3:
    from inspect import signature

scalar_types = (int, float, np.generic)
eps = np.finfo(np.float64).eps
solver_list = ['uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla']
invoker_list = solver_list[:]
invoker_list.append('pdfo')


class OptimizeResult(dict):
    """Result structure of the DFO algorithms.

    Attributes
    ----------
    x: ndarray, shape (n,)
        The (approximate) solution array.
    success: bool
       Flag indicating whether the optimizer exited successfully.
    status: int
        Flag characterizing the exit condition:
          0  The lower bound for the trust region radius is reached
          1  The target function value is achieved
          2  A trust region step failed to reduce the quadratic model
          3  The objective function has been evaluated `maxfev` times
 4, 7, 8, 9  Rounding errors become severe in the Fortran code
         13  All variables are fixed by the constraints
         -1  NaN occurs in `x`
         -2  The objective/constraint function returns NaN or nearly infinite values (only in the classical mode)
         -3  NaN occurs in the models
         -4  Constraints are infeasible

        exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran code but cannot be returned by pdfo or its
        solvers.
    message: str
        Message related to the exit condition flag. If `options['quiet']` is set to True, this message will not be
        printed.
    fun: float
        Returns the computed objective function value at the solution `x`.
    nfev: int
        Number of function evaluations.
    constrviolation: float
        Constraint violation at the solution `x`. It is set to 0 if the problem is unconstrained.
    fhist: ndarray, shape (nfev,)
        History of the objective function evaluations performed during the computation. Its size is `nfev` and is
        ordered by iteration. Its minimum among the feasible point should be fun.
    chist: ndarray, shape (nfev,)
        History of the constraint violations computed during the computation. If the problem is unconstrained, `chist`
        is set to None.
    method: str
        The name of the method that was used to solve the problem.
    constr_modified: bool
        An indicator specifying if the constraints have been modified during the algorithm (LINCOA may modify the
        constraints if it cannot find a feasible starting point).
    warnings: list
        A recording of every warning raised during the computation.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    # get an element of the optimization result structure
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError('The following attribute does not exist: {}.'.format(name))

    # set an element of the optimization result structure
    __setattr__ = dict.__setitem__

    # delete an element of the optimization result structure
    __delattr__ = dict.__delitem__

    # display the optimization result structure
    def __repr__(self):
        if self.keys():
            maxlength = max(map(len, self.keys())) + 1
            return '\n'.join([k.rjust(maxlength) + ': ' + repr(self[k]) for k in sorted(self.keys())])
        else:
            return self.__class__.__name__ + '()'


class Bounds:
    """Bound structure.

    The Numpy types of the arrays are converted to fit the Fortran backend and their sizes are checked.

    Attributes
    ----------
    lb: ndarray, shape (n,)
        The lower-bound vector of bound constraints.
    ub: ndarray, shape (n,)
        The upper-bound vector of bound constraints.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    def __init__(self, lb=None, ub=None):
        if isinstance(lb, scalar_types):
            lb = [lb]
        if isinstance(ub, scalar_types):
            ub = [ub]
        self.lb = np.asarray(lb if lb is not None else [], dtype=np.float64)
        self.ub = np.asarray(ub if ub is not None else [], dtype=np.float64)

        # reshape the flat matrices
        if len(self.lb.shape) > 1 and np.prod(self.lb.shape) == self.lb.size:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) > 1 and np.prod(self.ub.shape) == self.ub.size:
            self.ub = self.ub.reshape(self.ub.size)

        # check the length of the attributes
        if len(self.lb.shape) != 1 or len(self.ub.shape) != 1 or self.lb.size != self.ub.size:
            raise AttributeError('The sizes of the bounds are inconsistent; checks the shapes of the arrays.')


class LinearConstraint:
    """Linear constraint structure.

    The Numpy types of the arrays are converted to fit the Fortran backend and their sizes are checked.

    Attributes
    ----------
    A: ndarray, shape (m,n)
        The coefficient matrix of linear constraints.
    lb: ndarray, shape (n,)
        The lower-bound vector of linear constraints.
    ub: ndarray, shape (n,)
        The upper-bound vector of linear constraints.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    def __init__(self, a=None, lb=None, ub=None):
        if isinstance(a, scalar_types):
            a = [[a]]
        if isinstance(lb, scalar_types):
            lb = [lb]
        if isinstance(ub, scalar_types):
            ub = [ub]
        self.A = np.asarray(a if a is not None else [[]], dtype=np.float64, order='F')
        self.lb = np.asarray(lb if lb is not None else [], dtype=np.float64)
        self.ub = np.asarray(ub if ub is not None else [], dtype=np.float64)

        # reshape the flat matrices
        if len(self.lb.shape) > 1 and np.prod(self.lb.shape) == self.lb.size:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) > 1 and np.prod(self.ub.shape) == self.ub.size:
            self.ub = self.ub.reshape(self.ub.size)
        if len(self.lb.shape) == 0:
            self.lb = self.lb.reshape(lb.size)
        if len(self.ub.shape) == 0:
            self.ub = self.ub.reshape(ub.size)
        if len(self.A.shape) in [0, 1]:
            self.A = self.A.reshape((1, self.A.size))

        # check the length of the attributes
        if not (len(self.lb.shape) == 1 and len(self.A.shape) == 2 and len(self.ub.shape) == 1 and
                self.lb.size in [0, self.A.shape[0]] and self.ub.size in [0, self.A.shape[0]]) or \
                (self.lb.size == 0 and self.ub.size == 0 and self.A.size > 0):
            raise AttributeError('The sizes of linear constraints are inconsistent; check the shapes of the arrays.')


class NonlinearConstraint:
    """Nonlinear constraint structure.

    The Numpy types of the arrays are converted to fit the Fortran backend.

    Attributes
    ----------
    fun: callable
        The constraint function, which accepts a vector `x` at input and returns a vector of shape (n,).
    lb: ndarray, shape (n,)
        The lower-bound vector of the nonlinear constraints.
    ub: ndarray, shape (n,)
        The upper-bound vector of the nonlinear constraints.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    def __init__(self, fun, lb=None, ub=None):
        if not callable(fun):
            raise ValueError('The constraint function should be a callable object.')

        def float_fun(x):
            fx = np.asarray(fun(x))

            # scalars are converted to one-dimensional array
            if isinstance(fx, (int, float, np.generic)):
                fx = [fx]

            if not hasattr(fx, '__len__'):
                raise ValueError('The output of the constraint function has a wrong type.')

            fx = np.float64(fx)

            # when this function is called, the lower and upper bounds are well defined
            if fx.size > 0 and (len(self.lb.shape) > 1 or len(self.ub.shape) > 1 or
                                (self.lb.size == 0 and self.ub.size == 0) or self.lb.size not in [0, fx.size] or
                                self.ub.size not in [0, fx.size]):
                raise AttributeError(
                    'The size of the vector returned by the constraint function is inconsistent with the constraint'
                    'bounds; check the shapes of the arrays.')

            return fx

        if isinstance(lb, scalar_types):
            lb = [lb]
        if isinstance(ub, scalar_types):
            ub = [ub]

        self.fun = float_fun
        self.lb = np.asarray(lb if lb is not None else [], dtype=np.float64)
        self.ub = np.asarray(ub if ub is not None else [], dtype=np.float64)

        # reshape the flat matrices
        if len(self.lb.shape) > 1 and np.prod(self.lb.shape) == self.lb.size:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) > 1 and np.prod(self.ub.shape) == self.ub.size:
            self.ub = self.ub.reshape(self.ub.size)
        if len(self.lb.shape) == 0:
            self.lb = self.lb.reshape(lb.size)
        if len(self.ub.shape) == 0:
            self.ub = self.ub.reshape(ub.size)

        if len(self.lb.shape) != 1 or len(self.ub.shape) != 1 or \
                (self.lb.size != self.ub.size and self.lb.size > 1 and self.ub.size > 1):
            warnings.warn(
                'The sizes of the constraint bounds are inconsistent; check the shapes of the arrays.', Warning)


def prepdfo(fun, x0, args=(), method=None, bounds=None, constraints=(), options=None):
    """Pre-processing of the arguments.

    Parameters
    ----------
    fun: callable
        The objective function, which accepts a vector `x` at input and returns a scalar.
    x0: ndarray, shape (n,)
        The initial guess. The size of `x0` should be consistent with the objective function.
    args: tuple, optional
        The extra-arguments to pass to the objective function. For example,

            ``pdfo(fun, x0, args, options)``

        is equivalent to

            ``pdfo(lambda x: fun(x, args), x0, options=options)``

    method: str, optional
        The name of the Powell method that will be used.
    bounds: either ndarray of tuple with shape(n,2), or Bounds, optional
        Bound constraints of the problem. The bounds can be specified in two different ways:
            1. Instance of `Bounds` class.
            2. Sequence of (lb, ub) pairs for each element in `x`. To specify that `x[i]` is unbounded below, set
               `bounds[i, 0]` to -np.inf; set `bounds[i, 1]` to np.inf if `x[i]` is unbounded above.
    constraints: dict, LinearConstraint, NonlinearConstraint or list of them, optional
        Constraints of the problem. It can be one of the three cases below.
            1. A dictionary with fields:
                type: str
                    Constraint type: 'eq' for equality constraints and 'ineq' for inequality constraints.
                fun: callable
                    The constraint function.
            2. Instances of LinearConstraint or NonlinearConstraint.
            3. A list, each of whose elements can be a dictionary described in 1, an instance of LinearConstraint or an
               instance of NonlinearConstraint.
    options: dict, optional
        The options passed to the solver. It is a structure that contains optionally:
            rhobeg: float, optional
                Initial value of the trust region radius, which should be a positive scalar. `options['rhobeg']` should
                be typically set roughly to one tenth of the greatest expected change to a variable. By default, it is
                1.
            rhoend: float, optional
                Final value of the trust region radius, which should be a positive scalar. `options['rhoend']` should
                indicate typically the accuracy required in the final values of the variables. Moreover,
                `options['rhoend']` should be no more than `options['rhobeg']` and is by default 1e-6.
            maxfev: int, optional
                Upper bound of the number of calls of the objective function `fun`. Its value must be not less than
                `options['npt']`+1. By default, it is 500*n.
            npt: int, optional
                Number of interpolation points of each model used in Powell's Fortran code.
            ftarget: float, optional
                Target value of the objective function. If a feasible iterate achieves an objective function value lower
                or equal to `options['ftarget']`, the algorithm stops immediately. By default, it is -np.inf.
            scale: bool, optional
                Flag indicating whether to scale the problem. If it is True, the variables will be scaled according to
                the bounds constraints if any. By default, it is False.
            quiet: bool, optional
                Flag of quietness of the interface. If it is set to True, the output message will not be printed. This
                flag does not interfere with the warning and error printing.
            classical: bool, optional
                Flag indicating whether to call the classical Powell code or not. By default, it is False.
            debug: bool, optional
                Debugging flag. By default, it is False.
            chkfunval: bool, optional
                Flag used when debugging. If both `options['debug']` and `options['chkfunval']` are True, an extra
                function evaluation would be performed to check whether the returned objective function value is
                consistent with the returned `x`. By default, it is False.

    Returns
    -------
    The preprocessed `fun`, `x0`, `method`, `bounds`, `options`, and
    prob_info: dict
        An internal dictionary containing the problem information.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    fun_name = stack()[0][3]  # name of the current function
    list_warnings = []

    if len(stack()) < 3 or stack()[1][3].lower() not in invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(invoker_list)))
    invoker = stack()[1][3].lower()

    # saving of the raw data in prob_info before preprocessing
    prob_info = dict()
    prob_info['raw_data'] = \
        {'objective': fun, 'x0': x0, 'args': args, 'bounds': bounds, 'constraints': constraints, 'options': options}

    if len(stack()) >= 4 and stack()[2][3].lower() == 'pdfo':
        # the invoker is a solver called by pdfo, then prepdfo should have been called in pdfo
        prob_info['infeasible'] = False
        prob_info['nofreex'] = False
        return fun, x0, bounds, constraints, options, method, prob_info

    if fun is None:
        def fun(x_loc, *args_loc):
            return np.float64(0)

        warn_message = '{}: there is no objective function.'.format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    if args is not None and not hasattr(args, '__len__'):
        args = [args]

    def fun_c(x):
        try:
            fun_x = fun(x) if hasattr(args, '__len__') and len(args) == 0 else fun(x, *args)
        except TypeError:
            raise TypeError('{}: the number of parameters is inconsistent with `args`.'.format(invoker))

        if hasattr(fun_x, '__len__') and len(fun_x) == 1:
            fun_x = fun_x[0]
        elif (hasattr(fun_x, '__len__') or not isinstance(fun_x, scalar_types)) and fun_x is not None:
            raise ValueError('{}: the objective function should return a scalar.'.format(invoker))

        # if the objective function returns None, then we are solving a feasibility problem. The returned value of fun_c
        # will be interpreted as np.nan by the following statement
        return np.float64(fun_x)

    prob_info['raw_data']['objective'] = fun_c

    if isinstance(x0, scalar_types):
        x0 = [x0]
    if not hasattr(x0, '__len__'):
        raise ValueError('{}: the initial guess should be a scalar or a vector.'.format(invoker))
    try:
        x0_c = np.asarray(x0, dtype=np.float64)

        if x0_c.size == 0:
            raise ValueError('{}: the initial guess should not be empty.'.format(invoker))

        # reshape the initial guess
        if len(x0_c.shape) > 1 and np.prod(x0_c.shape) == x0_c.size:
            x0_c = x0_c.reshape(x0_c.size)
    except ValueError:
        raise ValueError('{}: the initial guess should contain only scalars.'.format(invoker))
    if len(x0_c.shape) > 1:
        raise ValueError('{}: the initial guess should be a scalar or a vector.'.format(invoker))

    # to be clear, the length of x0 is denoted lenx0 instead of n
    lenx0 = x0_c.size

    if method is None and invoker != 'pdfo':
        method = invoker
    if method is not None and not isinstance(method, str):
        raise ValueError('{}: the method name should be a string.'.format(invoker))

    lb, ub, infeasible, fixed_indices, fixed_values = _bounds_validation(invoker, bounds, lenx0)
    prob_info['raw_data']['bounds'] = (lb, ub)
    prob_info['infeasible_bound'] = infeasible
    prob_info['fixedx'] = fixed_indices
    prob_info['fixedx_value'] = fixed_values

    constraints_c, infeasible, trivial = _constraints_validation(invoker, constraints, lenx0)
    prob_info['raw_data']['constraints'] = constraints_c
    prob_info['infeasible_linear'] = infeasible
    prob_info['trivial_linear'] = trivial

    # after preprocessing the linear/bound constraints, the problem may turn out infeasible, or x may turn out fixed by
    # the bounds
    prob_info['infeasible'] = any(np.r_[prob_info['infeasible_bound'], prob_info['infeasible_linear']])
    prob_info['nofreex'] = all(prob_info['fixedx'])
    if prob_info['nofreex']:
        prob_info['constrv_fixedx'] = _constr_violation(invoker, prob_info['fixedx_value'], lb, ub, constraints_c)

    # reduce the problem if some variables are fixed by the bound constraints
    prob_info['raw_dim'] = lenx0
    prob_info['raw_type'] = _problem_type(lb, ub, constraints_c)
    prob_info['reduced'] = False

    if any(fixed_indices) and not (prob_info['nofreex'] or prob_info['infeasible']):
        fun_c, x0_c, lb, ub, constraints_c = _reduce_problem(fun_c, x0_c, lb, ub, constraints_c, fixed_indices)
        lenx0 = x0_c.size
        prob_info['reduced'] = True

    # problem dimension after reduction
    prob_info['refined_dim'] = lenx0

    # problem type after reduction
    prob_info['refined_type'] = _problem_type(lb, ub, constraints_c)

    # can the invoker handle the given problem? This should be done after the problem type has bee 'refined'
    if not _prob_solv_match(prob_info['refined_type'], invoker.lower()):
        if invoker.lower() == 'pdfo':
            raise SystemError(
                '{}: UNEXPECTED ERROR: problem and solver do not match; it should not happen when the invoker is pdfo '
                'or the problem is not a structure.'.format(fun_name))
        else:
            raise ValueError(
                '{}: a {} problem received; {} cannot solve '
                'it.'.format(fun_name, prob_info['refined_type'].replace('-', ' '), invoker))

    # validate and preprocess options, adopt default options if needed. This should be done after reducing the problem,
    # because BOBYQA requires rhobeg <= min(ub-lb)/2
    options_c, method = _options_validation(invoker, options, method, lenx0, lb, ub, list_warnings)

    # revise x0 for bound and linearly constrained problems. This is necessary for LINCOA, which accepts only a feasible
    # x0. Should we do this even if there are nonlinear constraints? For now, we do not, because doing so may
    # dramatically increase the infeasibility of x0 with respect to the nonlinear constraints
    if prob_info['refined_type'] in ['bound-constrained', 'linearly-constrained'] and not prob_info['nofreex'] and \
            not prob_info['infeasible']:
        x0_old = x0_c
        # another possibility for bound-constrained problems: xind = (x0 < lb) | (x0 > ub);
        # x0(xind) = (lb(xind) + ub(xind))/2;
        result = _project(x0_c, lb, ub, constraints_c)
        x0_c = result.x
        if np.linalg.norm(x0_old - x0_c) > eps * max(1.0, np.linalg.norm(x0_old)):
            warn_message = '{}: x0 is revised to satisfy the constraints.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    # scale the problem if necessary and if intended, x_before_scaling = scaling_factor.*x_after_scaling + shift
    # This should be done after revising x0, which can affect the shift
    prob_info['scaled'] = False
    if options_c['scale'] and not prob_info['nofreex'] and not prob_info['infeasible']:
        fun_c, x0_c, lb, ub, constraints_c, scaling_factor, shift, substantially_scaled = \
            _scale_problem(fun_c, x0_c, lb, ub, constraints_c, list_warnings)

        # scale and shift the problem so that
        #   1. for the variables that have both lower bound and upper bound, the bounds become [-1, 1]
        #   2. the other variables will be shifted so that the corresponding component of x0 becomes 0
        prob_info['scaled'] = True
        prob_info['scaling_factor'] = scaling_factor
        prob_info['shift'] = shift

        # if the problem is substantially scaled, then rhobeg and rhoend may need to be revised
        if substantially_scaled:
            options_c['rhobeg'] = np.float64(1.0)

    # select a solver if invoker == 'pdfo'
    if invoker.lower() == 'pdfo':
        if prob_info['refined_type'] == 'bound-constrained':
            # lb and ub will be used for defining rhobeg if bobyqa is selected
            prob_info['lb'] = lb
            prob_info['ub'] = ub

        method = _solver_selection(invoker, method, options_c, prob_info, list_warnings)

    prob_info['warnings'] = list_warnings

    # the refined data can be useful when debugging. It will be used in _postpdfo even if the debug mode is not enabled
    prob_info['refined_data'] = \
        {'objective': fun_c, 'x0': x0_c, 'lb': lb, 'ub': ub, 'constraints': constraints_c, 'options': options_c}

    if not options_c['debug']:
        # do not carry the raw data with us unless in debug mode
        del prob_info['raw_data']

    if prob_info['refined_type'] == 'bound-constrained' and invoker.lower() == 'pdfo':
        # if the invoker is pdfo, lb and ub may be used for defining rhobeg in case bobyqa is later selected as the
        # solver
        prob_info['lb'] = lb
        prob_info['ub'] = ub

    return fun_c, x0_c, {'lb': lb, 'ub': ub}, constraints_c, options_c, method, prob_info


def _bounds_validation(invoker, bounds, lenx0):
    """Validation and pre-processing of the bounds.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    bounds: ndarray of tuple, shape(n,2) or Bounds
        The same as in prepdfo.
    lenx0: integer
        The size of the problem.

    Returns
    -------
    The preprocessed `lb`, `ub`, and
    infeasible: ndarray, shape (n,)
        A boolean array indicating the infeasible constraints.
    fixed_indices: ndarray, shape (m,), m <= n
        An integer array of the i such that `lb[i] == ub[i]`.
    fixed_values: ndarray, shape (m,), m <= n
        A float array consisting of the `lb[i]` for the i in `fixed_indices`.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    if bounds is not None and (not hasattr(bounds, '__len__') or len(bounds) != 0):
        try:
            from scipy.optimize import Bounds as ScipyBounds

            if isinstance(bounds, ScipyBounds):
                bounds = Bounds(lb=bounds.lb, ub=bounds.ub)
        except ImportError:
            pass

        if not (isinstance(bounds, Bounds) or
                (hasattr(bounds, '__len__') and
                 all(hasattr(bound, '__len__') and len(bound) == 2 and
                     (isinstance(bound[0], scalar_types) or bound[0] is None) and
                     (isinstance(bound[1], scalar_types) or bound[1] is None) for bound in bounds))):
            raise ValueError(
                '{}: the bounds should be an instance of the `Bounds` class or a sequence of scalar '
                '2-tuples.'.format(invoker))

        try:
            if hasattr(bounds, '__len__'):
                bounds_c = np.asarray(bounds, order='F')
                lb = np.asarray(bounds_c[:, 0], dtype=np.float64, order='F')
                ub = np.asarray(bounds_c[:, 1], dtype=np.float64, order='F')
            else:
                lb = np.asarray(bounds.lb, dtype=np.float64, order='F')
                ub = np.asarray(bounds.ub, dtype=np.float64, order='F')
        except ValueError:
            raise ValueError('{}: the bound elements should be scalars.'.format(invoker))

        if lb.size != lenx0 or ub.size != lenx0:
            raise ValueError('{}: the bound size should be equal to the number of variables.'.format(invoker))
    else:
        lb = np.full(lenx0, -np.inf, dtype=np.float64)
        ub = np.full(lenx0, np.inf, dtype=np.float64)

    # NaN bounds are set to infinite values not to take them into account
    lb[np.isnan(lb)] = -np.inf
    ub[np.isnan(ub)] = np.inf

    infeasible = (lb > ub)
    fixed_indices = (np.abs(lb - ub) < 2 * eps)
    fixed_values = (lb[fixed_indices] + ub[fixed_indices]) / 2

    return lb, ub, infeasible, fixed_indices, fixed_values


def _constraints_validation(invoker, constraints, lenx0):
    """Validation and pre-processing of the constraints.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    constraints: dict, LinearConstraint, NonlinearConstraint or list of them
        The same as in prepdfo.
    lenx0: integer
        The size of the problem.

    Returns
    -------
    The preprocessed `constraints`, and
    infeasible: ndarray, shape (n,)
        A boolean array indicating the infeasible constraints.
    trivial: ndarray, shape (n,)
        A boolean array indicating the trivial constraints.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    if constraints is not None:
        try:
            from scipy.optimize import LinearConstraint as ScipyLinearConstraint
            from scipy.optimize import NonlinearConstraint as ScipyNonlinearConstraint
            linear_constraint_types = (LinearConstraint, ScipyLinearConstraint)
            nonlinear_constraint_types = (NonlinearConstraint, ScipyNonlinearConstraint)
        except ImportError:
            ScipyLinearConstraint = type(None)
            ScipyNonlinearConstraint = type(None)
            linear_constraint_types = LinearConstraint
            nonlinear_constraint_types = NonlinearConstraint

        # convert the constraints as a list
        if isinstance(constraints, dict) or not (hasattr(constraints, '__len__')):
            constraints_c = [constraints]
        else:
            constraints_c = constraints

        # create the linear/nonlinear sub-lists
        list_linear = []
        list_nonlinear = []
        for constraint in constraints_c:
            if isinstance(constraint, linear_constraint_types):
                A = np.asarray(constraint.A, dtype=np.float64, order='F')
                lb = np.asarray(constraint.lb, dtype=np.float64)
                ub = np.asarray(constraint.ub, dtype=np.float64)
                A[np.isnan(A)] = 0  # not to use those variables as constraints
                lb[np.isnan(lb)] = -np.inf
                ub[np.isnan(ub)] = np.inf
                list_linear.append(LinearConstraint(a=A, lb=lb, ub=ub))
            elif isinstance(constraint, nonlinear_constraint_types) or \
                    (isinstance(constraint, dict) and {'type', 'fun'} <= set(constraint.keys()) and
                     isinstance(constraint['type'], str) and constraint['type'] in ['eq', 'ineq'] and
                     callable(constraint['fun']) and
                     (python_version == 2 or len(signature(constraint['fun']).parameters) > 0)):
                if isinstance(constraint, nonlinear_constraint_types):
                    fun = constraint.fun
                    lb = np.asarray(constraint.lb, dtype=np.float64)
                    ub = np.asarray(constraint.ub, dtype=np.float64)
                    lb[np.isnan(lb)] = -np.inf
                    ub[np.isnan(ub)] = np.inf
                    list_nonlinear.append(NonlinearConstraint(fun=fun, lb=lb, ub=ub))
                else:
                    list_nonlinear.append(constraint)
            else:
                # the constraint is neither linear nor nonlinear
                raise ValueError(
                    "{}: the constraints should be instances of the `LinearConstraint` or `NonlinearConstraint` "
                    "classes, or a dictionary with field 'type' and 'fun'.".format(invoker))

        # linear constraints build
        a_linear = np.asarray([[]], dtype=np.float64, order='F').reshape(0, lenx0)
        lb_linear = np.asarray([], dtype=np.float64)
        ub_linear = np.asarray([], dtype=np.float64)

        for linear_constraint in list_linear:
            # the type of linear_constraint is surely LinearConstraint
            a_local = linear_constraint.A
            if a_local.size == 0:
                a_local = a_local.reshape(0, lenx0)
            if a_local.shape[1] != lenx0:
                raise ValueError(
                    '{}: the number of columns in A is inconsistent with the number of variables'.format(invoker))

            if linear_constraint.lb.size != 0:
                lb_local = linear_constraint.lb
            else:
                lb_local = np.full(a_local.shape[0], -np.inf)

            if linear_constraint.ub.size != 0:
                ub_local = linear_constraint.ub
            else:
                ub_local = np.full(a_local.shape[0], np.inf)

            a_linear = np.concatenate((a_linear, a_local), axis=0)
            lb_linear = np.r_[lb_linear, lb_local]
            ub_linear = np.r_[ub_linear, ub_local]

        # removing of the abnormal constraints and checking infeasibility
        if lb_linear.size == 0 and ub_linear.size == 0:
            infeasible = np.asarray([], dtype=bool)
            trivial = np.asarray([], dtype=bool)
        else:
            row_norm_inf = np.max(np.abs(a_linear), 1)
            zero = (row_norm_inf == 0)
            infeasible_zero = np.logical_or(np.logical_and(zero, lb_linear > 0), np.logical_and(zero, ub_linear < 0))
            trivial_zero = np.logical_or(np.logical_and(zero, lb_linear <= 0), np.logical_and(zero, ub_linear >= 0))
            row_norm_inf[zero] = 1.0
            lb_linear_norm = lb_linear / row_norm_inf
            ub_linear_norm = ub_linear / row_norm_inf
            lb_ub_or = np.logical_or(np.logical_and(np.isinf(lb_linear_norm), lb_linear_norm > 0),
                                     np.logical_and(np.isinf(ub_linear_norm), ub_linear_norm < 0))
            infeasible = np.logical_or(infeasible_zero, lb_ub_or)
            lb_ub_and = np.logical_and(np.logical_and(np.isinf(lb_linear_norm), lb_linear_norm < 0),
                                       np.logical_and(np.isinf(ub_linear_norm), ub_linear_norm > 0))
            trivial = np.logical_or(trivial_zero, lb_ub_and)
            a_linear = a_linear[np.logical_not(trivial), :]
            lb_linear = lb_linear[np.logical_not(trivial)]
            ub_linear = ub_linear[np.logical_not(trivial)]

        # nonlinear constraints build
        try:
            from .gethuge import gethuge
        except ImportError:
            import_error_so('gethuge')
        hugecon = gethuge('con')

        def fun_nonlinear(x):
            fun_x = np.asarray([], dtype=np.float64)

            for nonlinear_constraint in list_nonlinear:
                if isinstance(nonlinear_constraint, nonlinear_constraint_types):
                    constraint_x = nonlinear_constraint.fun(x)
                else:
                    constraint_x = nonlinear_constraint['fun'](x)

                if constraint_x is None:
                    # if the constraint function returned anything, we convert the default None value to NaN, that can
                    # be understood by Fortran
                    constraint_x = [np.nan]
                elif isinstance(constraint_x, scalar_types):
                    constraint_x = [constraint_x]

                if not hasattr(constraint_x, '__len__'):
                    raise ValueError('{}: the constraint function should return a vector or a scalar.'.format(invoker))

                constraint_x = np.asarray(constraint_x, dtype=np.float64)

                # use extreme barrier to cope with the 'hidden constraints'
                constraint_x[np.logical_or(np.isnan(constraint_x), constraint_x > hugecon)] = hugecon

                # this part is NOT extreme barrier. We replace extremely negative values of cineq (which leads to no
                # constraint violation) by -hugecon. Otherwise, NaN or Inf may occur in the interpolation models
                constraint_x[constraint_x < -hugecon] = -hugecon

                if len(constraint_x.shape) != 1:
                    raise ValueError('{}: the constraint function should return a vector or a scalar.'.format(invoker))

                lenm = constraint_x.size
                if isinstance(nonlinear_constraint, nonlinear_constraint_types):
                    if nonlinear_constraint.lb.size not in [0, lenm] or \
                            nonlinear_constraint.ub.size not in [0, lenm] or \
                            (nonlinear_constraint.lb.size == 0 and nonlinear_constraint.ub.size == 0):
                        raise ValueError(
                            '{}: the size of the vector returned by the constraint function is inconsistent with the '
                            'constraint bounds; check the shapes of the arrays.'.format(invoker))
                    if nonlinear_constraint.lb.size > 0 and nonlinear_constraint.ub.size > 0:
                        constraint_x = np.concatenate((nonlinear_constraint.lb - constraint_x,
                                                       constraint_x - nonlinear_constraint.ub))
                    elif nonlinear_constraint.lb.size > 0:
                        constraint_x = nonlinear_constraint.lb - constraint_x
                    else:
                        # necessarily, nonlinear_constraint.ub is well-defined
                        constraint_x = constraint_x - nonlinear_constraint.ub
                elif nonlinear_constraint['type'] == 'eq':
                    constraint_x = np.concatenate((constraint_x, -constraint_x))

                fun_x = np.concatenate((fun_x, constraint_x))

            return fun_x

        linear_constraints = None if len(list_linear) == 0 else LinearConstraint(a_linear, lb=lb_linear, ub=ub_linear)
        nonlinear_constraints = None if len(list_nonlinear) == 0 else {'type': 'ineq', 'fun': fun_nonlinear}
    else:
        linear_constraints = None
        nonlinear_constraints = None
        infeasible = np.asarray([], dtype=bool)
        trivial = np.asarray([], dtype=bool)

    return {'linear': linear_constraints, 'nonlinear': nonlinear_constraints}, infeasible, trivial


def _options_validation(invoker, options, method, lenx0, lb, ub, list_warnings):
    """Validation and pre-processing of the options.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    options: dict
        The same as in prepdfo.
    method: str
        The same as in prepdfo.
    lenx0: integer
        The size of the initial guess.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    list_warnings: list
        The same as in prepdfo.

    Returns
    -------
    The preprocessed `options` and `method`.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # NOTE: We treat field names case-sensitively

    if options is not None and not isinstance(options, dict):
        raise ValueError('{}: the options should be defined as a dictionary.'.format(invoker))

    # possible solvers
    options = dict() if options is None else options.copy()
    fun_name = stack()[0][3]  # name of the current function

    if invoker not in invoker_list:
        raise SystemError('{}: {} serves only {}'.format(fun_name, fun_name, ', '.join(invoker_list)))

    # which fields are specified?
    if options is not None:
        option_fields = options.keys()
    else:
        option_fields = []

    maxfev = 500 * lenx0
    rhobeg = 1  # the default rhobeg and rhoend will be revised for BOBYQA
    rhoend = 1e-6
    ftarget = -np.inf
    classical = False  # call the classical Powell code?
    scale = False  # scale the problem according to bounds or not
    quiet = True
    debugflag = False
    chkfunval = False

    # what is the solver? We need this information to decide which fields are 'known' (e.g., expected), and also to set
    # npt and rhobeg, rhoend
    if invoker == 'pdfo':
        if method is not None and method.lower() not in invoker_list:
            warn_message = '{}: unknown solver specified; {} will select one automatically.'.format(invoker, invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
            method = None
    else:  # invoker is in {'uobyqa', ..., 'cobyla'}
        if method is not None and method.lower() != invoker:
            warn_message = '{}: a solver different from {} is specified; it is ignored.'.format(invoker, invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        method = invoker

    # check unknown fields according to solver
    if method is not None and method.lower() in ['bobyqa', 'lincoa', 'newuoa']:
        known_field = ['npt', 'maxfev', 'rhobeg', 'rhoend', 'ftarget', 'classical', 'scale', 'quiet', 'debug',
                       'chkfunval', 'solver']
    else:
        known_field = ['maxfev', 'rhobeg', 'rhoend', 'ftarget', 'classical', 'scale', 'quiet', 'debug', 'chkfunval']
    unknown_field = list(set(option_fields).difference(set(known_field)))
    if len(unknown_field) > 0:
        warn_message = '{}: unknown option(s): {}; they are ignored.'.format(invoker, ', '.join(unknown_field))
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    # set default npt according to solver
    # If method == '', then invoker must be pdfo, and a solver will be selected later; when the solver is chosen, a
    # valid npt will be defined. So we do not need to consider the case with method == '' here. Note we have to take
    # maxfev into consideration when selecting the solver, because npt is at most maxfev-1!
    if method is None:
        npt = np.nan
    else:
        if method.lower() in ['bobyqa', 'lincoa', 'newuoa']:
            npt = 2 * lenx0 + 1
        elif method.lower() == 'uobyqa':
            npt = (lenx0 + 1) * (lenx0 + 2) / 2
        else:
            # the method is necessarily COBYLA
            npt = lenx0 + 1

    # revise default rhobeg and rhoend according to solver
    if method is not None and method.lower() == 'bobyqa':
        rhobeg_bobyqa = min(rhobeg, np.min(ub - lb) / 2)
        if ('scale' in option_fields and isinstance(options['scale'], bool) and not options['scale']) or \
                (not scale and not ('scale' in option_fields and isinstance(options['scale'], bool))):
            # if we are going to scale the problem later, then we keep the default value for rhoend; otherwise, we scale
            # it as follows
            rhoend = (rhoend / rhobeg) * rhobeg_bobyqa
        rhobeg = rhobeg_bobyqa

    # validate the user-specified options; adopt the default values if needed
    validated = False
    if 'npt' in option_fields and method is not None and method in ['bobyqa', 'lincoa', 'newuoa']:
        # only newuoa, bobyqa and lincoa accept an npt option
        if not isinstance(options['npt'], scalar_types) or options['npt'] < lenx0 + 2 or \
                options['npt'] > (lenx0 + 1) * (lenx0 + 2) / 2:
            warn_message = \
                '{}: invalid npt. for {}, it should be an integer and n+2 <= npt <= (n+1)*(n+2)/2; it is set to ' \
                '2n+1.'.format(invoker, method)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['npt'] has not got a valid value yet
        # for uobyqa and cobyla or empty solver, we adopt the 'default npt' defined above, although it will NOT be used
        # by the solver
        options['npt'] = npt
    # if prepdfo is called by pdfo, and if pdfo has been called without precising the method, npt will be set to np.nan.
    # If we do not check whether this value is np.nan before casting it into np.int32, it will lead to an error
    if not np.isnan(options['npt']):
        options['npt'] = np.int32(options['npt'])

    # validate options['maxfev']
    validated = False
    if 'maxfev' in option_fields:
        if not isinstance(options['maxfev'], scalar_types) or options['maxfev'] <= 0:
            warn_message = \
                '{}: invalid maxfev; it should be a positive integer; it is set to {}.'.format(invoker, maxfev)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif method is None and options['maxfev'] <= lenx0 + 1:
            options['maxfev'] = lenx0 + 2  # the smallest possible value
            validated = True
            warn_message = \
                '{}: invalid maxfev; it should be a positive integer at least n+2; it is set to n+2.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif method is not None and options['maxfev'] <= options['npt']:
            options['maxfev'] = options['npt'] + 1
            validated = True
            if method.lower() in ['bobyqa', 'lincoa', 'newuoa']:
                warn_message = \
                    '{}: invalid maxfev; {} requires maxfev > npt; it is set to npt+1.'.format(invoker, method)
            elif method.lower() == 'uobyqa':
                warn_message = \
                    '{}: invalid maxfev; {} requires maxfev > (n+1)*(n+2)/2; it is set to ' \
                    '(n+1)*(n+2)/2+1.'.format(invoker, method)
                options['maxfev'] = (lenx0 + 1) * (lenx0 + 2) / 2 + 1
            else:
                # the method is necessarily COBYLA
                warn_message = '{}: invalid maxfev; {} requires maxfev > n+1; it is set to n+2.'.format(invoker, method)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['maxfev'] has not got a valid value yet
        options['maxfev'] = max(maxfev, npt+1)
    options['maxfev'] = np.int32(options['maxfev'])

    # validate options['rhobeg']
    validated = False
    if 'rhobeg' in option_fields:
        if not isinstance(options['rhobeg'], scalar_types) or options['rhobeg'] <= 0:
            warn_message = \
                '{}: invalid rhobeg; it should be a positive number; it is set to ' \
                'max({}, rhoend).'.format(invoker, rhobeg)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif method is not None and method.lower() == 'bobyqa' and options['rhobeg'] > np.min(ub - lb) / 2:
            warn_message = \
                '{}: invalid rhobeg; {} requires rhobeg <= min(ub-lb)/2; it is set to ' \
                'min(ub-lb)/2.'.format(invoker, method)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
            options['rhobeg'] = np.min(ub - lb) / 2
            validated = True
        else:
            validated = True

    if not validated:  # options['rhobeg'] has not got a valid value yet
        if 'rhoend' in option_fields and isinstance(options['rhoend'], scalar_types):
            options['rhobeg'] = max(rhobeg, options['rhoend'])
        else:
            options['rhobeg'] = rhobeg
    options['rhobeg'] = np.float64(max(options['rhobeg'], eps))

    # validate options['rhoend']
    validated = False
    if 'rhoend' in option_fields:
        if not isinstance(options['rhoend'], scalar_types) or options['rhoend'] > options['rhobeg'] or \
                0 >= options['rhoend']:
            warn_message = \
                '{}: invalid rhoend; we should have rhobeg >= rhoend > 0; it is set to ' \
                '{}*rhobeg.'.format(invoker, rhoend / rhobeg)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['rhoend'] has not got a valid value yet
        if rhobeg > 0:
            options['rhoend'] = (rhoend / rhobeg) * options['rhobeg']
        else:
            options['rhoend'] = 0
    options['rhoend'] = np.float64(max(options['rhoend'], eps))

    # validate options['ftarget']
    validated = False
    if 'ftarget' in option_fields:
        if not isinstance(options['ftarget'], scalar_types):
            warn_message = '{}: invalid ftarget; it should be real number; it is set to {}.'.format(invoker, ftarget)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['ftarget'] has not got a valid value yet
        options['ftarget'] = ftarget
    options['ftarget'] = np.float64(options['ftarget'])

    # validate options['classical']
    validated = False
    if 'classical' in option_fields:
        if not isinstance(options['classical'], bool):
            warn_message = \
                '{}: invalid scale flag; it should be True or False; it is set to {}'.format(invoker, classical)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['classical'] has not got a valid value yet
        options['classical'] = scale
    options['classical'] = np.bool(options['classical'])
    if options['classical']:
        warn_message = \
            "{}: in classical mode, which is recommended only for research purpose; set options['classical']=False " \
            "to disable classical mode.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    # validate options['scale']
    validated = False
    if 'scale' in option_fields:
        if not isinstance(options['scale'], bool):
            warn_message = '{}: invalid scale flag; it should be True or False; it is set to {}'.format(invoker, scale)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['scale'] has not got a valid value yet
        options['scale'] = scale
    options['scale'] = np.bool(options['scale'])

    # validate options['quiet']
    validated = False
    if 'quiet' in option_fields:
        if not isinstance(options['quiet'], bool):
            warn_message = '{}: invalid quiet flag; it should be True or False; it is set to {}'.format(invoker, quiet)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['quiet'] has not got a valid value yet
        options['quiet'] = quiet
    options['quiet'] = np.bool(options['quiet'])

    # validate options['debug']
    validated = False
    if 'debug' in option_fields:
        if not isinstance(options['debug'], bool):
            warn_message = \
                '{}: invalid debug flag; it should be True or False; it is set to {}'.format(invoker, debugflag)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['debug'] has not got a valid value yet
        options['debug'] = debugflag
    options['debug'] = np.bool(options['debug'])
    if options['debug']:
        warn_message = "{}: in debug mode; set options['debug']=False to disable debug.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)
        if options['quiet']:
            options['quiet'] = False
            warn_message = "options['quiet'] is set to False because options['debug'] = True.".format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    # validate options['chkfunval']
    validated = False
    if 'chkfunval' in option_fields:
        if not isinstance(options['chkfunval'], bool):
            warn_message = \
                '{}: invalid chkfunval flag; it should be True or False; it is set to {}'.format(invoker, chkfunval)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif not options['debug']:
            warn_message = \
                '{}: chkfunval = True but debug = False; chkfunval is set to false; set both flags to true to check ' \
                'function values.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['chkfunval'] has not got a valid value yet
        options['chkfunval'] = chkfunval
    options['chkfunval'] = np.bool(options['chkfunval'])
    if options['chkfunval']:
        warn_message = \
            "{}: checking whether fx = fun(x) and possibly conval = con(x) at exit, which will cost an extra " \
            "function/constraint evaluation; set options['chkfunval'] = false to disable the check.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    return options, method


def _constr_violation(invoker, x, lb, ub, constraints):
    """Constraint violation of the iterate.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    x: ndarray, shape (n,)
        The current iterate.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: NonlinearConstraint
                The nonlinear constraints of the problem.

    Returns
    -------
    constr_violation: np.float64
        The constraint violation scalar.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    if not (hasattr(x, '__len__') and hasattr(lb, '__len__') and hasattr(ub, '__len__')):
        raise TypeError('{}: UNEXPECTED ERROR: the variable vector and the bounds should be vectors.'.format(invoker))

    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())):
        raise ValueError('{}: UNEXPECTED ERROR: the constraint should be defined as internal type.'.format(invoker))

    x = np.asarray(x, dtype=np.float64)
    lb = np.asarray(lb, dtype=np.float64)
    ub = np.asarray(ub, dtype=np.float64)
    if x.size != lb.size or x.size != ub.size:
        raise ValueError(
            '{}: UNEXPECTED ERROR: the sizes of the variable vector and the bounds are inconsistent.'.format(invoker))

    constr_violation = max(0, np.max(np.r_[lb - x, x - ub] / max(1, np.max(np.abs(np.r_[lb, ub])))))

    if constraints['linear'] is not None:
        a, b = _linear_constraints_constr(constraints['linear'])

        constr_violation = max(constr_violation, (np.dot(a, x) - b) / max(1, np.max(np.abs(b))))

    if constraints['nonlinear'] is not None:
        nonlinear = constraints['nonlinear']

        if not isinstance(nonlinear, dict) or not ({'type', 'fun'} <= set(nonlinear.keys())) or \
                nonlinear['type'] != 'ineq' or not callable(nonlinear['fun']) or \
                (python_version >= 3 and len(signature(nonlinear['fun']).parameters) == 0):
            raise ValueError('{}: UNEXPECTED ERROR: the nonlinear constraints are ill-defined.'.format(invoker))

        nlc = np.asarray(nonlinear['fun'](x), dtype=np.float64)
        constr_violation = max(constr_violation, np.max(nlc))

    return np.float64(constr_violation)


def _problem_type(lb, ub, constraints):
    """The type of the problem considering the constraints.

    Parameters
    ----------
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: NonlinearConstraint
                The nonlinear constraints of the problem.

    Returns
    -------
    constr_type: str
        The type of the constraints.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['prepdfo']
    if len(stack()) < 3 or stack()[1][3].lower() not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    if not (hasattr(lb, '__len__') and hasattr(ub, '__len__')):
        raise TypeError('{}: UNEXPECTED ERROR: the bounds should be vectors.'.format(invoker))

    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())):
        raise ValueError('{}: UNEXPECTED ERROR: the constraint should be defined as internal type.'.format(invoker))

    lb = np.asarray(lb, dtype=np.float64)
    ub = np.asarray(ub, dtype=np.float64)
    if lb.size != ub.size:
        raise ValueError('{}: UNEXPECTED ERROR: the bounds are inconsistent'.format(invoker))

    if constraints['nonlinear'] is not None:
        return 'nonlinearly-constrained'
    elif constraints['linear'] is not None and (constraints['linear'].lb.size > 0 or constraints['linear'].ub.size > 0):
        return 'linearly-constrained'
    elif (len(lb) > 0 and np.max(lb) > -np.inf) or (len(ub) > 0 and np.min(ub) < np.inf):
        return 'bound-constrained'
    else:
        return 'unconstrained'


def _reduce_problem(fun, x0, lb, ub, constraints, fixedx):
    """Reduction of the problem.

    Parameters
    ----------
    fun: callable
        The same as in prepdfo.
    x0: ndarray, shape (n,)
        The same as in prepdfo.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: NonlinearConstraint
                The nonlinear constraints of the problem.
    fixedx: ndarray, shape (n,)
        A boolean array relating the fixed variables.

    Returns
    -------
    The preprocessed `fun`, `x0`, `lb`, `ub`, `constraints`.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['postpdfo']
    if len(stack()) < 3 or stack()[1][3].lower() not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    if not (hasattr(x0, '__len__') and hasattr(lb, '__len__') and
            hasattr(ub, '__len__') and hasattr(fixedx, '__len__')):
        raise TypeError('{}: UNEXPECTED ERROR: the variable vector and the bounds should be vectors.'.format(invoker))

    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())):
        raise ValueError('{}: UNEXPECTED ERROR: the constraint should be defined as internal type.'.format(invoker))

    x0 = np.asarray(x0, dtype=np.float64)
    lb = np.asarray(lb, dtype=np.float64)
    ub = np.asarray(ub, dtype=np.float64)
    fixedx = np.asarray(fixedx, dtype=np.float64)
    if x0.size != lb.size or x0.size != ub.size or x0.size != fixedx.size:
        raise ValueError('{}: UNEXPECTED ERROR: the variable vector and the bounds are inconsistent'.format(invoker))

    # since this function is private and internal, the format of each parameter
    # has already been verified
    freex = np.logical_not(fixedx)
    fixedx_value = (lb[fixedx] + ub[fixedx]) / 2

    def fun_c(freex_value):
        return fun(_fullx(freex_value, fixedx_value, freex, fixedx))

    x0_c = x0[freex]
    lb_c = lb[freex]
    ub_c = ub[freex]

    constraints_c = {'linear': None, 'nonlinear': None}
    if constraints['linear'] is not None:
        a = constraints['linear'].A[:, freex]
        a_fixedx = np.dot(constraints['linear'].A[:, fixedx], fixedx_value)
        lb_lin = constraints['linear'].lb - a_fixedx
        ub_lin = constraints['linear'].ub - a_fixedx

        constraints_c['linear'] = LinearConstraint(a, lb=lb_lin, ub=ub_lin)

    if constraints['nonlinear'] is not None:
        def fun_constraints(freex_value):
            return constraints['nonlinear']['fun'](_fullx(freex_value, fixedx_value, freex, fixedx))

        constraints_c['nonlinear'] = {'type': 'ineq', 'fun': fun_constraints}

    return fun_c, x0_c, lb_c, ub_c, constraints_c


def _linear_constraints_constr(linear_constraint):
    """Construction of the matrices A and b such that A.x <= b.

    Parameters
    ----------
    linear_constraint: LinearConstraint
        The linear constraints.

    Returns
    -------
    a: ndarray
        The coefficient matrix of the inequality.
    b: ndarray
        The right-hand side vector of the inequality.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    fun_name = stack()[0][3]  # name of the current function

    if not isinstance(linear_constraint, LinearConstraint) or len(linear_constraint.A.shape) != 2 or \
            len(linear_constraint.lb.shape) != 1 or len(linear_constraint.ub.shape) != 1 or \
            linear_constraint.A.shape[0] != linear_constraint.lb.size or \
            linear_constraint.lb.size != linear_constraint.ub.size:
        raise ValueError(
            '{}: UNEXPECTED ERROR: the sizes of the coefficient matrix and the lower/upper-bound vectors are '
            'inconsistent.'.format(fun_name))

    a = linear_constraint.A
    lb = linear_constraint.lb
    ub = linear_constraint.ub
    is_lb = not np.logical_and(np.isinf(lb), lb < 0).all()
    is_ub = not np.logical_and(np.isinf(ub), ub > 0).all()
    if is_lb and is_ub:
        b = np.r_[ub, -lb]
        a = np.concatenate((a, -a), axis=0)
    elif is_lb:
        b = -lb
        a = -a
    else:
        b = ub

    return a, b


def _fullx(freex_value, fixedx_value, freex, fixedx):
    """Recover the full variable array from the free and fixed variable values.

    Parameters
    ----------
    freex_value: ndarray, shape (n,)
        A scalar array containing the free variable values.
    fixedx_value: ndarray, shape (m,)
        A scalar array containing the fixed variable values.
    freex: ndarray, shape (n,)
        A boolean array containing the free variable indices.
    fixedx: ndarray, shape (n,)
        A boolean array containing the fixed variable indices.

    Returns
    -------
    x: ndarray, shape (n + m,)
        The full variable array.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['prepdfo', 'postpdfo']
    if len(stack()) < 3 or stack()[1][3].lower() not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    if not (hasattr(freex_value, '__len__') and hasattr(fixedx_value, '__len__') and
            hasattr(freex, '__len__') and hasattr(fixedx, '__len__')):
        raise ValueError('{}: UNEXPECTED ERROR: the variable arrays have wrong types.'.format(invoker))

    freex_value = np.asarray(freex_value, dtype=np.float64)
    fixedx_value = np.asarray(fixedx_value, dtype=np.float64)
    freex = np.asarray(freex, dtype=np.float64)
    fixedx = np.asarray(fixedx, dtype=np.float64)
    if freex.size != fixedx.size or freex_value.size + fixedx_value.size != freex.size:
        raise ValueError('{}: UNEXPECTED ERROR: the variable vector lengths are inconsistent'.format(invoker))

    x = np.empty(freex_value.size + fixedx_value.size, dtype=np.float64)
    x[freex] = freex_value
    x[fixedx] = fixedx_value

    return x


def _prob_solv_match(problem_type, solver):
    """Check whether the problem type and the solver match.

    Parameters
    ----------
    problem_type: str
        The type of the problem.
    solver: str
        The name of the solver.

    Returns
    -------
    match: bool
        A flag indicating the matching of the problem type and the solver.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['_solver_selection', 'prepdfo']
    if len(stack()) < 3 or stack()[1][3].lower() not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    if solver not in ['bobyqa', 'cobyla', 'lincoa', 'newuoa', 'uobyqa', 'pdfo']:
        raise SystemError('{}: UNEXPECTED ERROR: {} is not a known solver.'.format(fun_name, solver))

    if not isinstance(problem_type, str):
        raise ValueError('{}: UNEXPECTED ERROR: the problem type should be a string.'.format(invoker))

    match = True

    if problem_type == 'unconstrained':
        # essentially do nothing. DO NOT remove this case. Otherwise, the case would be included in 'else', which is not
        # correct.
        pass
    elif problem_type == 'bound-constrained':
        if solver in ['newuoa', 'uobyqa']:
            match = False
    elif problem_type == 'linearly-constrained':
        if solver in ['bobyqa', 'newuoa', 'uobyqa']:
            match = False
    elif problem_type == 'nonlinearly-constrained':
        if solver in ['bobyqa', 'lincoa', 'newuoa', 'uobyqa']:
            match = False
    else:
        raise SystemError('{}: UNEXPECTED ERROR: {} is not a known problem type.'.format(fun_name, problem_type))

    return match


def _scale_problem(fun, x0, lb, ub, constraints, list_warnings):
    """Scale the problem.

    Parameters
    ----------
    fun: callable
        The same as in prepdfo.
    x0: ndarray, shape (n,)
        The same as in prepdfo.
        otherwise.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with
        fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: NonlinearConstraint
                The nonlinear constraints of the problem.
    list_warnings: list
        The same as in prepdfo.

    Returns
    -------
    The preprocessed `fun`, `x0`, `lb`, `ub`, `constraints`, and
    scaling_factor: ndarray, shape (n,)
        The scaling factor of each variables.
    shift: np.float64
        The shift of the variables: x_before_scaling = scaling_factor*x_after_scaling + shift.
    substantially_scaled: bool
        A flag indicating if the problem is substantially scaled.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['prepdfo']
    if len(stack()) < 3 or stack()[1][3].lower() not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    if not (hasattr(x0, '__len__') and hasattr(lb, '__len__') and hasattr(ub, '__len__')):
        raise TypeError('{}: UNEXPECTED ERROR: the initial guess and the bounds should be vectors.'.format(invoker))

    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())):
        raise ValueError('{}: UNEXPECTED ERROR: the constraint should be defined as internal type.'.format(invoker))

    x0 = np.asarray(x0, dtype=np.float64)
    lb = np.asarray(lb, dtype=np.float64)
    ub = np.asarray(ub, dtype=np.float64)
    if x0.size != lb.size or x0.size != ub.size:
        raise ValueError(
            '{}: UNEXPECTED ERROR: the sizes of the initial guess and the bounds are inconsistent'.format(invoker))

    if not isinstance(list_warnings, list):
        raise ValueError('{}: UNEXPECTED ERROR: the list of warnings is ill-defined'.format(invoker))

    # x_before_scaling = scaling_factor*x_after_scaling + shift

    # question: What about scaling according to the magnitude of x0, lb, ub, x0-lb, ub-x0?
    # This can be useful if lb and ub reflect the nature of the problem well, and x0 is a reasonable approximation to
    # the optimal solution. Otherwise, it may be a bad idea
    substantially_scaled_threshold = 4

    # we consider the problem substantially scaled_threshold if
    # max(scaling_factor)/min(scaling_factor) > substantially_scaled_threshold
    lenx0 = x0.size
    index_lub = np.logical_and(lb > -np.inf, ub < np.inf)
    scaling_factor = np.ones(lenx0, dtype=np.float64)
    shift = np.zeros(lenx0, dtype=np.float64)

    scaling_factor[index_lub] = (ub[index_lub] - lb[index_lub]) / 2
    shift[index_lub] = (ub[index_lub] + lb[index_lub]) / 2

    # shift x0 to 0 unless both lower and upper bounds are present
    shift[np.logical_not(index_lub)] = x0[np.logical_not(index_lub)]

    def fun_c(x):
        return np.float64(fun(scaling_factor * x + shift))

    x0_c = (x0 - shift) / scaling_factor
    lb_c = (lb - shift) / scaling_factor
    ub_c = (ub - shift) / scaling_factor

    constraints_c = {'linear': None, 'nonlinear': None}
    if constraints['linear'] is not None:
        a = constraints['linear'].A
        lb = constraints['linear'].lb
        ub = constraints['linear'].ub
        constraints_c['linear'] = \
            LinearConstraint(np.dot(a, np.diag(scaling_factor)), lb=lb - np.dot(a, shift), ub=ub - np.dot(a, shift))

    if constraints['nonlinear'] is not None:
        constraints_c['nonlinear'] = \
            {'type': 'ineq', 'fun': lambda x: constraints['nonlinear']['fun'](scaling_factor * x + shift)}

    if any(scaling_factor != 1):
        warn_message = \
            "{}: problem scaled according to bound constraints; do this only if the bounds reflect the scaling of " \
            "variables; if not, set options['scale']=false to disable scaling.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    substantially_scaled = False

    # if (max([scaling_factor; 1./scaling_factor]) > substantially_scaled_threshold)
    if np.max(scaling_factor) / np.min(scaling_factor) > substantially_scaled_threshold:
        substantially_scaled = True
        # this will affect the setting of rhobeg and rhoend: If x is substantially scaled, then rhobeg = 1,
        # rhoend = previously_defined_rhoend/previously_defined_rhobeg

    if np.min(scaling_factor) < eps:
        raise SystemError(
            '{}: UNEXPECTED ERROR: function scale_problem in pdfo returns a wrong scaling.'.format(invoker))

    return fun_c, x0_c, lb_c, ub_c, constraints_c, scaling_factor, shift, substantially_scaled


def _solver_selection(invoker, method, options, prob_info, list_warnings):
    """Select the solver that corresponds to the problem.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    method: str
        The same as in prepdfo.
    options: dict
        The same as in prepdfo.
    prob_info: dict
        An internal dictionary containing the problem information.
    list_warnings: list
        The same as in prepdfo.

    Returns
    -------
    solver: str
        The name of the solver that should be used.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['pdfo']

    if invoker not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    # validate invoker
    if not isinstance(invoker, str):
        raise ValueError('unknown: UNEXPECTED ERROR: invoker should be a string.')

    # validate invoker
    if method is not None and not isinstance(method, str):
        raise ValueError('{}: UNEXPECTED ERROR: method should be a string.'.format(invoker))

    # validate options
    option_fields = {'maxfev', 'rhobeg', 'rhoend'}
    if options is None or not isinstance(options, dict) or not (option_fields <= set(options.keys())) or \
            not isinstance(options['maxfev'], scalar_types) or not isinstance(options['rhobeg'], scalar_types) or \
            not isinstance(options['rhoend'], scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: options should be a valid dictionary.'.format(invoker))

    # validate prob_info
    prob_info_fields = {'refined_type', 'refined_dim'}
    if prob_info is None or not isinstance(prob_info, dict) or not (prob_info_fields <= set(prob_info.keys())) or \
            not isinstance(prob_info['refined_type'], str) or not isinstance(prob_info['refined_dim'], scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: prob_info should be a valid dictionary.'.format(invoker))

    # validate list_warnings
    if not hasattr(list_warnings, '__len__'):
        raise ValueError('{}: UNEXPECTED ERROR: list_warnings should be a list.'.format(invoker))

    solver = method
    ptype = prob_info['refined_type']
    n = prob_info['refined_dim']

    if solver is None or not _prob_solv_match(ptype, solver):
        if solver is not None:
            # do not complain if solver is None
            warn_message = \
                '{}: {} cannot solve a {} problem; {} will select a solver ' \
                'automatically.'.format(invoker, solver, ptype.replace('-', ' '), invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

        if ptype == 'unconstrained':
            if 2 <= n <= 8 and options['maxfev'] >= (n + 1) * (n + 2) / 2:
                solver = 'uobyqa'  # does not need options['npt']
            elif options['maxfev'] <= n + 2:  # options['maxfev'] == n + 2
                solver = 'cobyla'  # does not need options['npt']
            else:
                # Interestingly, we note in our test that LINCOA outperformed NEWUOA on unconstrained CUTEst problems
                # when the dimension was not large (i.e., <=50) or the precision requirement was not high (i.e.,
                # >= 1e-5). Therefore, it is worthwhile to try LINCOA when an unconstrained problem is given.
                # Nevertheless, for the moment, we set the default solver for unconstrained problems to be newuoa.
                solver = 'newuoa'
                options['npt'] = min(2 * n + 1, options['maxfev'] - 1)
        elif ptype == 'bound-constrained':
            if options['maxfev'] <= n + 2:
                solver = 'cobyla'  # does not need options['npt']
            else:
                solver = 'bobyqa'
                options['npt'] = min(2 * n + 1, options['maxfev'] - 1)
                if {'lb', 'ub'} <= set(prob_info.keys()):
                    rhobeg_bobyqa = min(options['rhobeg'], np.min(prob_info['ub'] - prob_info['lb']) / 2)
                    options['rhoend'] = (options['rhoend'] / options['rhobeg']) * rhobeg_bobyqa
                    options['rhobeg'] = max(rhobeg_bobyqa, eps)
                    options['rhoend'] = max(options['rhoend'], eps)
        elif ptype == 'linearly-constrained':
            if options['maxfev'] <= n + 2:
                solver = 'cobyla'  # does not need options['npt']
            else:
                solver = 'lincoa'
                options['npt'] = min(2 * n + 1, options['maxfev'] - 1)
        elif ptype == 'nonlinearly-constrained':
            solver = 'cobyla'  # does not need options['npt']
        else:
            raise SystemError("{}: UNEXPECTED ERROR: unknown problem type '{}' received.".format(fun_name, ptype))

    if solver not in solver_list or not _prob_solv_match(ptype, solver):
        raise SystemError("{}: UNEXPECTED ERROR: invalid solver '{}' selected.".format(fun_name, solver))

    return solver


def _project(x0, lb, ub, constraints, options=None):
    """Projection of the initial guess onto the feasible set.

    Parameters
    ----------
    x0: ndarray, shape (n,)
        The same as in prepdfo.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with
        fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: dict
                The nonlinear constraints of the problem. When ``_project`` is called, the nonlinear constraints are
                None.
    options: dict, optional

    Returns
    -------
    result: OptimizeResult
        The result of the projection.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['prepdfo']
    if len(stack()) < 3 or stack()[1][3].lower() not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(invoker_list)))
    invoker = stack()[1][3].lower()

    # validate x0
    if isinstance(x0, scalar_types):
        x0_c = [x0]
    elif hasattr(x0, '__len__'):
        x0_c = x0
    else:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should be a vector.'.format(invoker))
    try:
        x0_c = np.asarray(x0_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should contain only scalars.'.format(invoker))
    if len(x0_c.shape) != 1:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should be a vector.'.format(invoker))
    lenx0 = x0_c.size

    # validate lb
    if isinstance(lb, scalar_types):
        lb_c = [lb]
    elif hasattr(lb, '__len__'):
        lb_c = lb
    else:
        raise ValueError('{}: UNEXPECTED ERROR: lb should be a vector.'.format(invoker))
    try:
        lb_c = np.asarray(lb_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: lb should contain only scalars.'.format(invoker))
    if len(lb_c.shape) != 1 or lb.size != lenx0:
        raise ValueError('{}: UNEXPECTED ERROR: the size of lb is inconsistent with x0.'.format(invoker))

    # validate ub
    if isinstance(ub, scalar_types):
        ub_c = [ub]
    elif hasattr(ub, '__len__'):
        ub_c = ub
    else:
        raise ValueError('{}: UNEXPECTED ERROR: ub should be a vector.'.format(invoker))
    try:
        ub_c = np.asarray(ub_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: ub should contain only scalars.'.format(invoker))
    if len(ub_c.shape) != 1 or ub.size != lenx0:
        raise ValueError('{}: UNEXPECTED ERROR: the size of ub is inconsistent with x0.'.format(invoker))

    # validate constraints
    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())) or \
            not (isinstance(constraints['linear'], LinearConstraint) or constraints['linear'] is None):
        # the nonlinear constraints will not be taken into account in this function and are, therefore, not validated
        raise ValueError('{}: UNEXPECTED ERROR: The constraints are ill-defined.'.format(invoker))

    # validate options
    if options is not None and not isinstance(options, dict):
        raise ValueError('{}: UNEXPECTED ERROR: The options should be a dictionary.'.format(invoker))

    # projection onto the feasible set
    max_con = 1e20  # decide whether an inequality constraint can be ignored

    if constraints['linear'] is None:
        # direct projection onto the bound constraints
        x_proj = np.min((np.max((x0_c, lb_c), axis=0), ub_c), axis=0)
        return OptimizeResult(x=x_proj)
    elif np.equal(constraints['linear'].lb, constraints['linear'].ub).all() and np.max(lb_c) <= -max_con and \
            np.min(ub_c) >= max_con:
        # the linear constraints are all equality constraints
        try:
            from scipy.linalg import lstsq

            a = constraints['linear'].A
            b = constraints['linear'].lb
            xi, _, _, _ = lstsq(a, b - np.dot(a, x0_c))
            x_proj = np.min((np.max((x0_c + xi, lb_c), axis=0), ub_c), axis=0)

            return OptimizeResult(x=x_proj)
        except ImportError:
            # we can try to project the initial guess onto the feasible set by solving the associated optimization
            # problem.
            pass

    if constraints['linear'] is not None:
        try:
            from scipy.optimize import minimize
            from scipy.optimize import Bounds as ScipyBounds
            from scipy.optimize import LinearConstraint as ScipyLinearConstraint

            linear = constraints['linear']

            # to be more efficient, SciPy asks to separate the equality and the inequality constraints into two
            # different LinearConstraint structures
            pc_args_ineq, pc_args_eq = dict(), dict()
            pc_args_ineq['A'], pc_args_eq['A'] = np.asarray([[]], order='F'), np.asarray([[]], order='F')
            pc_args_ineq['A'] = pc_args_ineq['A'].reshape(0, linear.A.shape[1])
            pc_args_eq['A'] = pc_args_eq['A'].reshape(0, linear.A.shape[1])
            pc_args_ineq['lb'], pc_args_eq['lb'] = np.asarray([]), np.asarray([])
            pc_args_ineq['ub'], pc_args_eq['ub'] = np.asarray([]), np.asarray([])

            for i in range(linear.lb.size):
                if linear.lb[i] != linear.ub[i]:
                    pc_args_ineq['A'] = np.concatenate((pc_args_ineq['A'], linear.A[i:i+1, :]), axis=0)
                    pc_args_ineq['lb'] = np.r_[pc_args_ineq['lb'], linear.lb[i]]
                    pc_args_ineq['ub'] = np.r_[pc_args_ineq['ub'], linear.ub[i]]
                else:
                    pc_args_eq['A'] = np.concatenate((pc_args_eq['A'], linear.A[i:i+1, :]), axis=0)
                    pc_args_eq['lb'] = np.r_[pc_args_eq['lb'], linear.lb[i]]
                    pc_args_eq['ub'] = np.r_[pc_args_eq['ub'], linear.ub[i]]

            if pc_args_ineq['lb'].size > 0 and pc_args_eq['lb'].size > 0:
                project_constraints = [ScipyLinearConstraint(**pc_args_ineq), ScipyLinearConstraint(**pc_args_eq)]
            elif pc_args_ineq['lb'].size > 0:
                project_constraints = ScipyLinearConstraint(**pc_args_ineq)
            else:
                project_constraints = ScipyLinearConstraint(**pc_args_eq)

            ax_ineq = np.dot(pc_args_ineq['A'], x0_c)
            ax_eq = np.dot(pc_args_eq['A'], x0_c)
            if np.greater(ax_ineq, pc_args_ineq['ub']).any() or np.greater(pc_args_ineq['lb'], ax_ineq).any() or \
                    np.not_equal(ax_eq, pc_args_eq['lb']).any() or \
                    np.greater(x0_c, ub_c).any() or np.greater(lb_c, x0_c).any():
                return minimize(lambda x: np.dot(x - x0_c, x - x0_c) / 2, x0_c, jac=lambda x: (x - x0_c),
                                bounds=ScipyBounds(lb_c, ub_c), constraints=project_constraints)
            else:
                # do not perform any projection if the initial guess is feasible
                return OptimizeResult(x=x0_c)

        except ImportError:
            return OptimizeResult(x=x0_c)

    return OptimizeResult(x=x0_c)


def _augmented_linear_constraint(n, bounds, constraints):
    """Concatenate bound and linear constraints into one constraint.

    Parameters
    ----------
    n: int
        The size of the problem.
    bounds: dict
        The bounds of the problem, defined as a dictionary with fields:
            lb: ndarray, shape (n,)
                The lower bounds of the problem.
            ub: ndarray, shape (n,)
                The upper bounds of the problem.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with
        fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: NonlinearConstraint
                The nonlinear constraints of the problem.

    Returns
    -------
    a_aug: ndarray, shape (m,n)
        The coefficient matrix of the augmented linear constraints.
    b_aug: ndarray, shape (m,)
        The right-hand side vector of the augmented linear constraints.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    fun_name = stack()[0][3]  # name of the current function
    if len(stack()) < 3 or stack()[1][3].lower() not in invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(invoker_list)))
    invoker = stack()[1][3].lower()

    if not isinstance(n, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: the size of the problem should be a scalar.'.format(invoker))
    try:
        n = np.int32(n)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: the size of the problem should be an integer.'.format(invoker))

    if not isinstance(bounds, dict) or not ({'lb', 'ub'} <= set(bounds.keys())) or \
            not hasattr(bounds['lb'], '__len__') or len(bounds['lb']) != n or not hasattr(bounds['ub'], '__len__') or \
            len(bounds['ub']) != n:
        raise ValueError('{}: UNEXPECTED ERROR: the bounds are ill-defined.'.format(invoker))

    if not isinstance(constraints, dict) or not ({'linear'} <= set(constraints.keys())) or \
            not (constraints['linear'] is None or isinstance(constraints['linear'], LinearConstraint)):
        raise ValueError('{}: UNEXPECTED ERROR: the constraints are ill-defined.'.format(invoker))

    idmatrix = np.eye(n)
    lb, ub = bounds['lb'], bounds['ub']
    lb_kept_indices = np.logical_not(np.logical_and(np.isinf(lb), lb < 0))
    ub_kept_indices = np.logical_not(np.logical_and(np.isinf(ub), ub > 0))
    alb = idmatrix[lb_kept_indices, :]
    aub = idmatrix[ub_kept_indices, :]

    # reshape the empty matrices to avoid concatenate exception
    if aub.size == 0:
        aub = aub.reshape(0, n)
    if alb.size == 0:
        alb = alb.reshape(0, n)

    # remove infinite bounds
    lb, ub = lb[lb_kept_indices], ub[ub_kept_indices]

    # construction of the augmented matrices
    if constraints['linear'] is None:
        aineq = np.array([[]], dtype=np.float64)
        bineq = np.array([], dtype=np.float64)
    else:
        aineq, bineq = _linear_constraints_constr(constraints['linear'])
    if aineq.size == 0:
        aineq = aineq.reshape(0, n)

    a_aug = np.concatenate((aineq, -alb, aub), axis=0)
    b_aug = np.concatenate((bineq, -lb, ub), axis=0)
    if not (a_aug.size == 0 and b_aug.size == 0) and \
            not (len(a_aug.shape) == 2 and a_aug.shape[0] == b_aug.size and a_aug.shape[1] == n):
        raise SystemError('{}: UNEXPECTED ERROR: invalid augmented linear constraints.'.format(invoker))

    return a_aug, b_aug


def postpdfo(x, fx, exitflag, output, method, nf, fhist, options, prob_info, constrviolation=0, chist=None):
    """Post-processing of the arguments.

    Parameters
    ----------
    x: ndarray, shape (n,)
        The (approximate) solution array.
    fx: np.float64
        The value of the objective function at `x`.
    exitflag: int
        The flag indicating the exit condition of the solver.
    output: dict
        A dictionary containing all the fields that should be added to the optimizer result.
    method: str, optional
        The name of the Powell method that was used.
    nf: int
        The number of function evaluations.
    fhist: ndarray, shape (m,)
        The history of every objective function evaluations.
    options: dict, optional
        The same as in pdfo. It has been preprocessed by prepdfo. `options['quiet']` will be used.
    prob_info: dict
        An internal dictionary containing the problem information.
    constrviolation: np.float64, optional
        The constraint violation at `x`.
    chist: ndarray, shape (m,), optional
        The history of constraint violations.

    Returns
    -------
    result: OptimizeResult
        The results of the solver, represented as an instance of ``OptimizeResult``.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    # With extreme barrier (implemented when options['classical']=False), all the function values that are NaN or larger
    # than hugefun are replaced by hugefun; all the constraint values that are NaN or larger than hugecon are replaced
    # by hugecon. hugefun and hugecon are defined in pdfoconst.F, and can be obtained by gethuge.
    try:
        from .gethuge import gethuge
    except ImportError:
        import_error_so('gethuge')
    hugefun = gethuge('fun')
    hugecon = gethuge('con')

    # possible solvers
    fun_name = stack()[0][3]  # name of the current function

    if len(stack()) < 3 or stack()[1][3].lower() not in invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(invoker_list)))
    invoker = stack()[1][3].lower()

    # validate x
    if not hasattr(x, '__len__') and \
            not isinstance(x, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: x should be a scalar or a vector.'.format(invoker))
    try:
        x_c = np.asarray(x, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: x should contain only scalars.'.format(invoker))
    if len(x_c.shape) > 1:
        raise ValueError('{}: UNEXPECTED ERROR: x should be a vector.'.format(invoker))

    # validate fx
    if hasattr(fx, '__len__') and len(fx) == 1:
        fx_c = np.float64(fx[0])
    else:
        fx_c = np.float64(fx)
    if not isinstance(fx_c, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: fx should be a scalar.'.format(invoker))

    # validate exitflag
    if not isinstance(exitflag, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: exitflag should be a scalar.'.format(invoker))
    exitflag_c = np.int32(exitflag)
    if exitflag_c != exitflag:
        raise ValueError('{}: UNEXPECTED ERROR: exitflag should not be a floating number.'.format(invoker))

    # validate output
    if output is None or not isinstance(output, dict):
        raise ValueError('{}: UNEXPECTED ERROR: output should be a valid dictionary.'.format(invoker))

    # validate method
    if method is None or not isinstance(method, str):
        raise ValueError('{}: UNEXPECTED ERROR: method should be a string.'.format(invoker))

    # validate nf
    if not isinstance(nf, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: nf should be a scalar.'.format(invoker))
    nf_c = np.int32(nf)
    if nf_c != nf:
        raise ValueError('{}: UNEXPECTED ERROR: nf should not be a floating number.'.format(invoker))

    # validate fhist
    if not hasattr(fhist, '__len__') and not isinstance(fhist, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: fhist should be a scalar of a vector.'.format(invoker))
    try:
        fhist_c = np.asarray(fhist[:nf], dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: fhist should contain nf scalars.'.format(invoker))
    if len(fhist_c.shape) != 1:
        raise ValueError('{}: UNEXPECTED ERROR: fhist should be a vector.'.format(invoker))

    # validate constrviolation
    if not np.isnan(constrviolation) and not isinstance(constrviolation, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: constrviolation should be a scalar.'.format(invoker))
    if np.isnan(constrviolation):
        constrviolation_c = constrviolation
    else:
        constrviolation_c = np.float64(constrviolation)

    # validate chist
    if not (chist is None and (method in ['pdfo', 'newuoa', 'uobyqa'] or nf == 0)) and \
            not hasattr(chist, '__len__') and not isinstance(chist, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: chist should be a scalar or a vector.'.format(invoker))
    if chist is None:
        chist_c = chist
    else:
        try:
            chist_c = np.asarray(chist[:nf], dtype=np.float64)
        except ValueError:
            raise ValueError('{}: UNEXPECTED ERROR: chist should contain nf scalars.'.format(invoker))
    if chist_c is not None and len(chist_c.shape) != 1:
        raise ValueError('{}: UNEXPECTED ERROR: chist should be a vector.'.format(invoker))

    # if the invoker is a solver called by pdfo, then let pdfo do the
    # postprocessing
    output['x'] = x_c
    output['fun'] = fx_c
    output['status'] = exitflag_c
    output['success'] = exitflag_c in [0, 1]
    if len(stack()) >= 4 and stack()[2][3].lower() == 'pdfo':
        output['nfev'] = nf_c
        output['constrviolation'] = constrviolation_c
        output['fhist'] = fhist_c
        output['chist'] = chist_c

        return OptimizeResult(**output)

    # if the solver is not called by pdfo (can be pdfo directly), perform the post-processing
    # validate options
    option_fields = {'quiet', 'debug', 'classical', 'chkfunval'}
    if options is None or not isinstance(options, dict) or not (option_fields <= set(options.keys())) or \
            not isinstance(options['quiet'], (bool, np.bool)) or not isinstance(options['debug'], (bool, np.bool)) or \
            not isinstance(options['classical'], (bool, np.bool)) or \
            not isinstance(options['chkfunval'], (bool, np.bool)):
        raise ValueError('{}: UNEXPECTED ERROR: options should be a valid dictionary.'.format(invoker))

    # Manage the extreme barriers
    if not options['classical']:
        if (fhist_c > hugefun).any() or np.isnan(fhist_c).any():
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} returns an fhist with NaN or values larger than hugefun={}; this is '
                'impossible with extreme barrier.'.format(invoker, method, hugefun))
        elif fhist_c.size > 0 and np.max(fhist_c) == hugecon:
            warn_message = \
                '{}: extreme barrier is invoked; function values that are NaN or larger than hugefun={} are replaced ' \
                'by hugefun.'.format(invoker, hugefun)
            warnings.warn(warn_message, Warning)
            output['warnings'].append(warn_message)

        if method == 'cobyla' and chist_c is not None and hasattr(chist_c, '__len__') and \
                ((chist_c > hugecon).any() or np.isnan(chist_c).any()):
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} returns an chist with NaN or values larger than hugecon={}; this is '
                'impossible with extreme barrier.'.format(invoker, method, hugecon))
        elif chist_c is not None and chist_c.size > 0 and np.max(chist_c) == hugecon:
            warn_message = '{}: extreme barrier is invoked; function values that are NaN or larger than hugecon={} ' \
                           'are replaced by hugecon.'.format(invoker, hugecon)
            warnings.warn(warn_message, Warning)
            output['warnings'].append(warn_message)

    # validate prob_info
    prob_info_fields = \
        {'infeasible', 'warnings', 'scaled', 'reduced', 'fixedx', 'fixedx_value', 'refined_type', 'raw_type',
         'infeasible_linear', 'infeasible_bound'}
    if prob_info is None or not isinstance(prob_info, dict) or not (prob_info_fields <= set(prob_info.keys())) or \
            not isinstance(prob_info['infeasible'], (bool, np.bool)) or \
            not hasattr(prob_info['warnings'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, str), prob_info['warnings'])) or \
            not isinstance(prob_info['scaled'], (bool, np.bool)) or \
            not isinstance(prob_info['reduced'], (bool, np.bool)) or not hasattr(prob_info['fixedx'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, (bool, np.bool_)), prob_info['fixedx'])) or \
            not hasattr(prob_info['fixedx_value'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, scalar_types), prob_info['fixedx_value'])) or \
            not hasattr(prob_info['infeasible_linear'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, (bool, np.bool_)), prob_info['infeasible_linear'])) or \
            not hasattr(prob_info['infeasible_bound'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, (bool, np.bool_)), prob_info['infeasible_bound'])) or \
            not isinstance(prob_info['refined_type'], str) or not isinstance(prob_info['raw_type'], str):
        raise ValueError('{}: UNEXPECTED ERROR: prob_info should be a valid dictionary.'.format(invoker))

    prob_info_keys = prob_info.keys()
    if prob_info['scaled']:
        prob_info_fields_scaled = {'scaling_factor', 'shift'}
        if not (prob_info_fields_scaled <= set(prob_info_keys)) or \
                not hasattr(prob_info['scaling_factor'], '__len__') or \
                not all(map(lambda pi: isinstance(pi, scalar_types), prob_info['scaling_factor'])) or \
                not hasattr(prob_info['shift'], '__len__') or \
                not all(map(lambda pi: isinstance(pi, scalar_types), prob_info['shift'])):
            raise ValueError(
                '{}: UNEXPECTED ERROR: prob_info should contain scaling factors if the problem has been '
                'scaled.'.format(invoker))
    prob_info_c = prob_info.copy()

    # validate the value of the inputs
    if nf_c < 0 or (nf_c == 0 and not prob_info['infeasible'] and exitflag_c > 0):
        raise ValueError('{}: UNEXPECTED ERROR: {} returns a nf <= 0 unexpectedly'.format(invoker, method))

    # the problem was (possibly) scaled, scale it back
    # The scaling affects constrviolation when there are bound constraint. Hence constrviolation has to be recalculated
    # so that it equals the constraint violation of the returned x with respect to the original problem.  Ideally, chist
    # should also be recalculated. However, it is impossible because we do not save the history of x. Therefore, when
    # prob_info['scaled'] == True, chist is not the history of constraint violation of the original problem but the
    # scaled one. It it not consistent with constrviolation. Without saving of history of x, we cannot do better.
    # Before recalculating constrviolation, save the one returned by the solver, because it will be used in debug mode
    # when checking whether fx is consistent with fhist and chist. See the definition of fhistf for details.
    constrv_returned = constrviolation_c
    if prob_info_c['scaled']:
        # First calculate the residuals of the linear constraints. This must be calculated before x is scaled back.
        # Otherwise, we would have to scale also the linear constraints back to get the correct residuals
        linear = prob_info_c['refined_data']['constraints']['linear']
        if linear is not None:
            Ax = np.dot(linear.A, x_c)
            r = np.r_[linear.lb - Ax, Ax - linear.ub]
            b = np.r_[-linear.lb, linear.ub]
        else:
            r, b = np.asarray([np.nan]), np.asarray([np.nan])

        # scale x back
        x_c = prob_info_c['scaling_factor'] * x_c + prob_info_c['shift']
        output['x'] = x_c

        # scale the bounds back
        lb = prob_info_c['scaling_factor'] * prob_info_c['refined_data']['lb']
        lb += prob_info_c['shift']
        ub = prob_info_c['scaling_factor'] * prob_info_c['refined_data']['ub']
        ub += prob_info_c['shift']

        # We only need to calculate constrviolation for lincoa and cobyla, because uobyqa and newuoa do not handle
        # constrained problems, while bobyqa is a feasible method and should return constrviolation = 0 regardless of
        # the scaling unless something goes wrong.
        if method == 'lincoa':
            # LINCOA returns a relative constraint violation
            conv_n = np.concatenate((r, lb - x_c, x_c - ub))
            conv_n = np.nanmax((np.zeros_like(conv_n), conv_n), axis=0)
            conv_d = np.abs(np.concatenate((b, lb, ub)))
            conv_d = np.nanmax((np.ones_like(conv_d), conv_d), axis=0)
            constrviolation_c = np.max(conv_n / conv_d)
        else:
            nlc = np.asarray([-np.inf], dtype=np.float64)
            if 'nlc' in output.keys():
                nlc = np.asarray(output['nlc'], dtype=np.float64)
            conv = np.concatenate((r, lb - x_c, x_c - ub, nlc))
            if np.isnan(conv).all():
                constrviolation_c = np.nan
            else:
                constrviolation_c = np.nanmax(np.append(conv, 0))

    # the problem was (possibly) reduced, get the full x
    if prob_info_c['reduced']:
        x_c = _fullx(x_c, prob_info_c['fixedx_value'], np.logical_not(prob_info_c['fixedx']), prob_info_c['fixedx'])

    # set output.{nf, constrviolation, fhist, chist, method}
    output['nfev'] = nf_c
    output['constrviolation'] = constrviolation_c
    output['fhist'] = fhist_c
    output['chist'] = chist_c
    output['method'] = method

    # revise constrviolation and chist according to problem type
    max_c = 0 if chist_c is None or chist_c.size == 0 else np.nanmax(chist_c)
    if prob_info_c['refined_type'] == 'unconstrained' and (constrviolation_c > 0 or max_c > 0):
        raise ValueError(
            '{}: UNEXPECTED ERROR: {} returns positive constrviolations for an unconstrained '
            'problem.'.format(invoker, method))

    if prob_info_c['raw_type'] == 'unconstrained':
        if 'constrviolation' in output.keys():
            del output['constrviolation']
        if 'chist' in output.keys():
            del output['chist']
    elif nf_c > 0 and prob_info_c['refined_type'] == 'unconstrained' and prob_info_c['raw_type'] != 'unconstrained':
        output['constrviolation'] = np.float64(0)
        output['chist'] = np.zeros(nf_c)

    # record the returned message
    if exitflag_c == 0:
        output['message'] = \
            'Return from {} because the lower bound for the trust region radius is reached.'.format(method)
    elif exitflag_c == 1:
        output['message'] = 'Return from {} because the target function value is achieved.'.format(method)
    elif exitflag_c == 2:
        output['message'] = \
            'Return from {} because a trust region step has failed to reduce the quadratic model.'.format(method)
    elif exitflag_c == 3:
        output['message'] = \
            'Return from {} because the objective function has been evaluated maxfev times.'.format(method)
    elif exitflag_c == 4:
        output['message'] = 'Return from {} because of much cancellation in a denominator.'.format(method)
    elif exitflag_c == 5:
        output['message'] = 'Return from {} because npt is not in the required interval.'.format(method)
    elif exitflag_c == 6:
        output['message'] = \
            'Return from {} because one of the differences xu(i) - xl(i) is less than 2*rhobeg.'.format(method)
    elif exitflag_c == 7:
        output['message'] = 'Return from {} because rounding errors are becoming damaging.'.format(method)
    elif exitflag_c == 8:
        output['message'] = 'Return from {} because rounding errors prevent reasonable changes to x.'.format(method)
    elif exitflag_c == 9:
        output['message'] = 'Return from {} because the denominator of the updating formula is zero.'.format(method)
    elif exitflag_c == 10:
        output['message'] = 'Return from {} because n should not be less than 2.'.format(method)
    elif exitflag_c == 11:
        output['message'] = 'Return from {} because maxfev is less than npt+1.'.format(method)
    elif exitflag_c == 12:
        output['message'] = 'Return from {} because the gradient of a constraint is zero.'.format(method)
    elif exitflag_c == 13:
        output['message'] = 'Return from {} because all the variables are fixed by the bounds.'.format(method)
    elif exitflag_c == -1:
        output['message'] = 'Return from {} because NaN occurs in x.'.format(method)
    elif exitflag_c == -2:
        if method == 'cobyla':
            output['message'] = \
                'Return from {} because the objective function returns an NaN or nearly infinite value, or the ' \
                'constraints return a NaN.'.format(method)
        else:
            output['message'] = \
                'Return from {} because the objective function returns an NaN or nearly infinite value.'.format(method)
    elif exitflag_c == -3:
        output['message'] = 'Return from {} because NaN occurs in the models.'.format(method)
    elif exitflag_c == -4:
        if np.any(prob_info['infeasible_linear']):
            output['InfeasibleLinear'] = np.where(prob_info['infeasible_linear'])[0]

        if np.any(prob_info['infeasible_bound']):
            output['InfeasibleBound'] = np.where(prob_info['infeasible_bound'])[0]

        output['message'] = 'Return from {} because the constraints are infeasible.'.format(method)
    else:
        raise ValueError('{}: UNEXPECTED ERROR: {} returns an invalid exitflag {}.'.format(invoker, method, exitflag_c))

    if 'warnings' in output.keys():
        warning_list_output = output['warnings']
        del output['warnings']

        if not hasattr(warning_list_output, '__len__'):
            raise SystemError('{}: UNEXPECTED ERROR: the warnings should be defined as a list.'.format(invoker))
    else:
        warning_list_output = []

    if 'warnings' in prob_info_c.keys():
        warning_list = list(prob_info_c['warnings'])
        warning_list.extend(list(warning_list_output))
    else:
        warning_list = []

    if len(warning_list) > 0:
        output['warnings'] = warning_list

    # more careful checks about fx, constrviolation, fhist and chist.
    # We do this only if the coe is in debug mode but not in classical mode. The classical mode cannot pass these checks
    if options['debug'] and not options['classical'] and nf_c > 0:
        if 'raw_data' not in prob_info_keys:
            raise ValueError("{}: UNEXPECTED ERROR: 'raw_data' should be a field of prob_info".format(invoker))

        # check whether fx is 'optimal'
        fhistf = fhist_c
        if method in ['bobyqa', 'lincoa', 'cobyla']:
            fhistf = fhistf[chist_c <= max(constrv_returned, 0)]

        if np.isnan(fhistf).all():
            min_f = np.nan
        else:
            min_f = np.nanmin((fx_c, np.nanmin(fhistf)))

        if fx != min_f and not (np.isnan(fx) and np.isnan(min_f)) and method != 'lincoa' and \
                'constr_modified' in output.keys() and output['constr_modified']:
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} returns an fhist that does not match nf or fx'.format(invoker, method))

        # check whether constrviolation is correct
        cobyla_prec = np.float64(1e-12)
        lincoa_prec = np.float64(1e-14)

        # COBYLA cannot ensure fx=fun(x) or conval=con(x) due to rounding errors. Instead of checking the equality, we
        # check whether the relative error is within cobyla_prec.
        constrviolation = np.float64(0)
        if 'constrviolation' in output.keys():
            constrviolation = output['constrviolation']
        if method == 'bobyqa' and np.nanmax((constrviolation, np.nanmax(chist_c))) > 0:
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} is a feasible solver yet it returns positive '
                'constrviolations.'.format(invoker, method))

        if (method == 'lincoa' and not output['constr_modified']) or method == 'cobyla':
            linear = prob_info_c['raw_data']['constraints']['linear']
            lb, ub = prob_info_c['raw_data']['bounds']

            if linear is not None:
                try:
                    Ax = np.dot(linear.A, x_c)
                    r = np.r_[linear.lb - Ax, Ax - linear.ub]
                    b = np.r_[-linear.lb, linear.ub]
                except ValueError:
                    raise ValueError(
                        '{}: UNEXPECTED ERROR: the linear constraints are no more consistent'.format(invoker))
            else:
                b, r = np.asarray([np.nan]), np.asarray([np.nan])

            if method == 'lincoa':
                # LINCOA returns a relative constraint violation
                conv_n = np.concatenate((r, lb - x_c, x_c - ub))
                conv_n = np.nanmax((np.zeros_like(conv_n), conv_n), axis=0)
                conv_d = np.abs(np.concatenate((b, lb, ub)))
                conv_d = np.nanmax((np.ones_like(conv_d), conv_d), axis=0)
                conv = np.max(conv_n / conv_d)
            else:
                nlc = np.asarray([-np.inf], dtype=np.float64)
                if 'nlc' in output.keys():
                    nlc = np.asarray(output['nlc'], dtype=np.float64)
                conv = np.concatenate(([0], r, lb - x_c, x_c - ub, nlc))
                conv = np.nanmax(conv)
                # COBYLA returns an absolute constraint violation (there is nothing to compare with, because con is a
                # black-box)

            if not (np.isnan(conv) and np.isnan(constrviolation_c)) and \
                    not (np.isinf(conv) and np.isinf(constrviolation_c)) and \
                    not (abs(constrviolation_c - conv) <= lincoa_prec * max(1, abs(conv)) and method == 'lincoa') and \
                    not (abs(constrviolation_c - conv) <= cobyla_prec * max(1, abs(conv)) and method == 'cobyla'):
                raise ValueError(
                    '{}: UNEXPECTED ERROR: {} returns a CONSTRVIOLATION that does not match x.'.format(invoker, method))

            if np.isnan(fx_c):
                cf = chist_c[np.isnan(fhist_c)]
            else:
                cf = chist_c[fhist_c == fx_c]
            if (cf != constrv_returned).all() and not (np.isnan(constrv_returned) and np.isnan(cf).all()):
                raise ValueError(
                    '{}: UNEXPECTED ERROR: {} returns a CONSTRVIOLATION that does not match '
                    'chist.'.format(invoker, method))

        if options['chkfunval']:
            # check whether fx = fun(x)
            fun_x = prob_info_c['raw_data']['objective'](x_c)
            if np.isnan(fun_x) or (fun_x > hugefun):
                fun_x = hugefun
                # due to extreme barrier (implemented when options['classical']=False), all the function values that are
                # NaN or larger than hugefun are replaced by hugefun.

            # it seems that COBYLA can return fx~=fun(x) due to rounding errors. Therefore, we cannot use "fx != fun_x"
            # to check COBYLA
            if not (np.isnan(fx_c) and np.isnan(fun_x)) and \
                    not (fx_c == fun_x or (abs(fun_x - fx_c) <= cobyla_prec * max(1, abs(fx_c))) and
                         method == 'cobyla'):
                raise ValueError(
                    '{}: UNEXPECTED ERROR: {} returns an fx that does not match x.'.format(invoker, method))

            # check whether nlc = nonlinear(x)
            nonlinear = prob_info_c['raw_data']['constraints']['nonlinear']
            if nonlinear is not None:
                if 'nlc' in output.keys():
                    nlc = output['nlc']
                else:
                    nlc = np.array([], dtype=np.float64)

                if not hasattr(nlc, '__len__'):
                    raise ValueError('{}: UNEXPECTED ERROR: nlc should be recorded as a ndarray.'.format(invoker))

                nlcx = nonlinear['fun'](x_c)
                nlcx[np.logical_or(np.isnan(nlcx), nlcx > hugecon)] = hugecon

                # This part is NOT extreme barrier. We replace extremely negative values of cineq (which leads to no
                # constraint violation) by -hugecon. Otherwise, NaN or Inf may occur in the interpolation models
                nlcx[nlcx < -hugecon] = -hugecon

                max_x = 0 if nlcx.size == 0 else np.nanmax(nlcx)

                if nlc.size != nlcx.size or not np.array_equal(np.isnan(nlc), np.isnan(nlcx)) or \
                        (np.isnan(nlc).any() and (np.abs(nlc - nlcx) > cobyla_prec * max(1, max_x)).any()):
                    raise ValueError(
                        '{}: UNEXPECTED ERROR: {} returns a con(x) that does not match x.'.format(invoker, method))

    if method == 'lincoa' and 'constr_modified' in output.keys():
        del output['constr_modified']

    if 'nlc' in output.keys():
        del output['nlc']

    if not options['quiet']:
        print(output['message'], end='\n\n')

    return OptimizeResult(**output)


def import_error_so(missing_file=None):
    """Raise an error, if an import failed.

    Parameters
    ----------
    missing_file: str
        The name of the missing file.

    Returns
    -------
    None

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    invoker = stack()[1][3].lower()
    if len(stack()) < 3 or invoker not in invoker_list:
        raise SystemError('The function `{}` should not be called manually'.format(stack()[0][3]))

    if missing_file is None:
        missing_file = invoker

    system_os = platform.system()
    system_known = True
    if system_os.lower() in ['darwin', 'linux']:
        system_os = 'unix'
    elif system_os.lower() == 'windows':
        system_os = 'win'
    else:
        system_known = False

    if system_known:
        raise ImportError(
            '{} is missing, please execute `setup.py` and add the interface path to the PYTHONPATH variable (see '
            '`README_py_{}.txt`).'.format(missing_file, system_os))
    else:
        raise ImportError(
            '{} is missing, please execute `setup.py` and add the interface path to the PYTHONPATH variable (see '
            'the README file corresponding to your system).'.format(missing_file))
