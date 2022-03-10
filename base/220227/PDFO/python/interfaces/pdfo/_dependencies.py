# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

import sys
import warnings
from inspect import stack

import numpy as np

python_version = sys.version_info.major
if python_version >= 3:
    # From Python 3, we can check the signature of the function and ensure that
    # the objective and constraint function are correctly defined.
    from inspect import signature

# All the accepted scalar types; np.generic correspond to all NumPy types.
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
         14  A linear feasibility problem has been received and solved
         15  A linear feasibility problem has been received but was not solved
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
    constr_value: ndarray or list of ndarrays
        Values of the constraint functions at the returned `x`. It can be one of the two cases below depending on how
        the `constraints` variable is specified at the input.
        1. If `constraints` is a dictionary or an instance of NonlinearConstraint or LinearConstraint, then
           `constr_value` is an ndarray, whose value is `constraints['fun'](x)`, `constraints.fun(x)`, or
           `constraints.A*x`.
        2. If `constraints` is a list of dictionaries or instances of NonlinearConstraint or LinearConstraint, then
           `constr_value` is a list of ndarrays described in 1, each of which is the value of the corresponding
           component in `constraints`.
        If a nonlinear constraint is trivial (i.e., it has -inf as th lower bound and +inf as th upper bound), then its
        is represented by NaN in constr_value, because such constraints are not evaluated during the computation.
        Trivial nonlinear constraints with unknown dimensions (for example, {'type': 'eq', 'fun': None} or
        NonlinearConstraint(fun, None, None)) are represented by empty arrays in constr_value.
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

    def __delattr__(self, key):
        r"""
        Delete the attribute ``key``. This method is called by the built-in
        command ``del obj.key``, where ``obj`` is an instance of this class.
        """
        super(OptimizeResult, self).__delitem__(key)

    def __dir__(self):
        r"""
        Return the list of the names of the attributes in the current scope.
        This method is called when the built-in function ``dir(obj)`` is called,
        where ``obj`` is an instance of this class.
        """
        return list(self.keys())

    def __getattr__(self, name):
        r"""
        Return the value of the attribute ``name``. This method raises an
        `AttributeError` exception if such an attribute does not exist and is
        called when the built-in access ``obj.name`` is called, where ``obj`` is
        an instance of this class.
        """
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    def __repr__(self):
        r"""
        Return a string representation of an instance of this class, which looks
        like a valid Python expression that can be used to recreate an object
        with the same value (given an appropriate environment).
        """
        items = sorted(self.items())
        args = ', '.join(k + '=' + repr(v) for k, v in items)
        return self.__class__.__name__ + '(' + args + ')'

    def __setattr__(self, key, value):
        r"""
        Assign the value ``value`` to the attribute ``key``. This method is
        called when the built-in assignment ``obj.key = value`` is made, where
        ``obj`` is an instance of this class.
        """
        super(OptimizeResult, self).__setitem__(key, value)

    def __str__(self):
        r"""
        Return an informal string representation of an instance of this class,
        which is designed to be nicely printable. This method is called when the
        built-in functions ``format(obj)`` or ``print(obj)`` is called, where
        ``obj`` is an instance of this class.
        """
        if self.keys():
            items = sorted(self.items())
            width = max(map(len, self.keys())) + 1
            return '\n'.join(k.rjust(width) + ': ' + str(v) for k, v in items)
        else:
            return self.__class__.__name__ + '()'


class Bounds:
    """Bound structure.

    Bounds(lb, ub) specifies a bound constraint

    lb <= x <= ub,

    where x is an n-dimensional vector.

    Attributes
    ----------
    lb: ndarray, shape (n,)
        The lower-bound vector of the constraint. Set ``lb`` to ``-numpy.inf``
        to disable these bounds.
    ub: ndarray, shape (n,)
        The upper-bound vector of the constraint. Set ``ub`` to ``numpy.inf``
        to disable these bounds.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    def __init__(self, lb=None, ub=None):
        # The scalars are converted to arrays with one elements, to process a single type.
        if isinstance(lb, scalar_types):
            lb = [lb]
        if isinstance(ub, scalar_types):
            ub = [ub]

        if not ((lb is None or hasattr(lb, '__len__')) and (ub is None or hasattr(ub, '__len__'))):
            raise AttributeError('The bounds lb and ub should be vectors.')

        # Either lb, ub or both can be set to None or [] not to precise the bound.
        if (lb is None or len(lb) == 0) and ub is not None and len(ub) > 0:
            self.lb = np.full(len(ub), -np.inf, dtype=np.float64)
            self.ub = np.asarray(ub, dtype=np.float64)
        elif lb is not None and len(lb) > 0 and (ub is None or len(ub) == 0):
            self.lb = np.asarray(lb, dtype=np.float64)
            self.ub = np.full(len(lb), np.inf, dtype=np.float64)
        else:
            self.lb = np.asarray(lb if lb is not None else [], dtype=np.float64)
            self.ub = np.asarray(ub if ub is not None else [], dtype=np.float64)

        # Reshape the flat matrices.
        if len(self.lb.shape) > 1 and np.prod(self.lb.shape) == self.lb.size:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) > 1 and np.prod(self.ub.shape) == self.ub.size:
            self.ub = self.ub.reshape(self.ub.size)

        # Broadcast the bounds if necessary.
        if self.lb.size == 1:
            self.lb = np.copy(np.broadcast_to(self.lb, self.ub.shape))
        elif self.ub.size == 1:
            self.ub = np.copy(np.broadcast_to(self.ub, self.lb.shape))

        # Check the length of the attributes.
        if len(self.lb.shape) != 1 or len(self.ub.shape) != 1 or self.lb.size != self.ub.size:
            raise AttributeError('The sizes of the bounds are inconsistent; checks the shapes of the arrays.')

    def __repr__(self):
        return '{}({}, {})'.format(type(self).__name__, self.lb, self.ub)


class LinearConstraint:
    """Linear constraint structure.

    LinearConstraint(A, lb, ub) specifies a linear constraint

    lb <= A*x <= ub,

    where x is an n-dimensional vector.

    Attributes
    ----------
    A: ndarray, shape (m,n)
        The coefficient matrix of the constraint.
    lb: ndarray, shape (m,)
        The lower-bound vector of the constraint. Set ``lb`` to ``-numpy.inf``
        to disable these bounds.
    ub: ndarray, shape (m,)
        The upper-bound vector of the constraint. Set ``ub`` to ``numpy.inf``
        to disable these bounds.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    def __init__(self, a=None, lb=None, ub=None):
        # The scalars are converted to arrays with one elements, to process a single type.
        if isinstance(a, scalar_types):
            a = [[a]]
        if isinstance(lb, scalar_types):
            lb = [lb]
        if isinstance(ub, scalar_types):
            ub = [ub]

        if not ((lb is None or hasattr(lb, '__len__')) and (ub is None or hasattr(ub, '__len__'))):
            raise AttributeError('The bounds lb and ub should be vectors.')

        # Either lb, ub or both can be set to None or [] not to precise the bound.
        self.A = np.asarray(a if a is not None else [[]], dtype=np.float64, order='F')
        if (lb is None or len(lb) == 0) and ub is not None and len(ub) > 0:
            self.lb = np.full(len(ub), -np.inf, dtype=np.float64)
            self.ub = np.asarray(ub, dtype=np.float64)
        elif lb is not None and len(lb) > 0 and (ub is None or len(ub) == 0):
            self.lb = np.asarray(lb, dtype=np.float64)
            self.ub = np.full(len(lb), np.inf, dtype=np.float64)
        else:
            self.lb = np.asarray(lb if lb is not None else [], dtype=np.float64)
            self.ub = np.asarray(ub if ub is not None else [], dtype=np.float64)

        # If any NaN is detected, the it should not altered the constraint. The NaN will be most likely generated during
        # the conversion of the type of the arrays to np.float64, since np.float64(None) is NaN.
        self.A[np.isnan(self.A)] = 0  # not to use those variables as constraints
        self.lb[np.isnan(self.lb)] = -np.inf
        self.ub[np.isnan(self.ub)] = np.inf

        # Reshape the flat matrices.
        if len(self.lb.shape) > 1 and np.prod(self.lb.shape) == self.lb.size:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) > 1 and np.prod(self.ub.shape) == self.ub.size:
            self.ub = self.ub.reshape(self.ub.size)
        if len(self.lb.shape) == 0:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) == 0:
            self.ub = self.ub.reshape(self.ub.size)
        if len(self.A.shape) in [0, 1]:
            self.A = self.A.reshape((1, self.A.size))

        # Broadcast the bounds if necessary.
        if self.lb.size == 1:
            self.lb = np.copy(np.broadcast_to(self.lb, self.ub.shape))
        elif self.ub.size == 1:
            self.ub = np.copy(np.broadcast_to(self.ub, self.lb.shape))

        # If no bounds have been provided, infinite values should be considered.
        if self.lb.size == 0 and self.ub.size == 0:
            self.lb = np.full(self.A.shape[0], -np.inf, dtype=np.float64)
            self.ub = np.full(self.A.shape[0], np.inf, dtype=np.float64)

        # Check the length of the attributes.
        if not (len(self.lb.shape) == 1 and len(self.A.shape) == 2 and len(self.ub.shape) == 1 and
                self.lb.size in [0, self.A.shape[0]] and self.ub.size in [0, self.A.shape[0]]) or \
                (self.lb.size == 0 and self.ub.size == 0 and self.A.size > 0):
            raise AttributeError('The sizes of linear constraints are inconsistent; check the shapes of the arrays.')

    def __repr__(self):
        return '{}({}, {}, {})'.format(type(self).__name__, self.A, self.lb, self.ub)


class NonlinearConstraint:
    """Nonlinear constraint structure.

    NonlinearConstraint(fun, lb, ub) specifies a nonlinear constraint

    lb <= fun(x) <= ub.

    Attributes
    ----------
    fun: callable
        The constraint function, which accepts a vector `x` at input and returns a vector of shape (m,).
    lb: ndarray, shape (m,)
        The lower-bound vector of the constraint. Set ``lb`` to ``-numpy.inf``
        to disable these bounds.
    ub: ndarray, shape (m,)
        The upper-bound vector of the constraint. Set ``ub`` to ``numpy.inf``
        to disable these bounds.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """

    def __init__(self, fun, lb=None, ub=None):
        if not callable(fun) and fun is not None and not (hasattr(fun, '__len__') and len(fun) == 0):
            # If fun is defined as None or [], the constraint should not be considered
            raise ValueError('The constraint function should be a callable object.')

        # The output of the nonlinear constraint function should be an array containing floating point numbers.
        def float_fun(x):
            if callable(fun):
                fx = np.asarray(fun(x))
            else:
                # If fun is not defined as callable, this function will never be called, except if any bounds is
                # defined, which should raise an exception.
                fx = []

            # Scalars are converted to one-dimensional array.
            if isinstance(fx, (int, float, np.generic)):
                fx = [fx]

            if not hasattr(fx, '__len__'):
                raise ValueError('The output of the constraint function has a wrong type.')

            fx = np.float64(fx)

            # When this function is called, the lower and upper bounds are well defined.
            if fx.size > 0 and (len(self.lb.shape) > 1 or len(self.ub.shape) > 1 or
                                (self.lb.size == 0 and self.ub.size == 0) or self.lb.size not in [0, fx.size] or
                                self.ub.size not in [0, fx.size]):
                raise AttributeError(
                    'The size of the vector returned by the constraint function is inconsistent with the constraint '
                    'bounds; check the shapes of the arrays.')

            return fx

        # The scalars are converted to arrays with one elements, to process a single type.
        if isinstance(lb, scalar_types):
            lb = [lb]
        if isinstance(ub, scalar_types):
            ub = [ub]

        if not ((lb is None or hasattr(lb, '__len__')) and (ub is None or hasattr(ub, '__len__'))):
            raise AttributeError('The bounds lb and ub should be vectors.')

        # Either lb, ub or both can be set to None or [] not to precise the bound.
        self.fun = float_fun
        if (lb is None or len(lb) == 0) and ub is not None and len(ub) > 0:
            self.lb = np.full(len(ub), -np.inf, dtype=np.float64)
            self.ub = np.asarray(ub, dtype=np.float64)
        elif lb is not None and len(lb) > 0 and (ub is None or len(ub) == 0):
            self.lb = np.asarray(lb, dtype=np.float64)
            self.ub = np.full(len(lb), np.inf, dtype=np.float64)
        else:
            self.lb = np.asarray(lb if lb is not None else [], dtype=np.float64)
            self.ub = np.asarray(ub if ub is not None else [], dtype=np.float64)

        # If any NaN is detected, the it should not altered the constraint. The NaN will be most likely generated during
        # the conversion of the type of the arrays to np.float64, since np.float64(None) is NaN.
        self.lb[np.isnan(self.lb)] = -np.inf
        self.ub[np.isnan(self.ub)] = np.inf

        # Reshape the flat matrices.
        if len(self.lb.shape) > 1 and np.prod(self.lb.shape) == self.lb.size:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) > 1 and np.prod(self.ub.shape) == self.ub.size:
            self.ub = self.ub.reshape(self.ub.size)
        if len(self.lb.shape) == 0:
            self.lb = self.lb.reshape(self.lb.size)
        if len(self.ub.shape) == 0:
            self.ub = self.ub.reshape(self.ub.size)

        # Broadcast the bounds if necessary.
        if self.lb.size == 1:
            self.lb = np.copy(np.broadcast_to(self.lb, self.ub.shape))
        elif self.ub.size == 1:
            self.ub = np.copy(np.broadcast_to(self.ub, self.lb.shape))

        if len(self.lb.shape) != 1 or len(self.ub.shape) != 1 or \
                (self.lb.size != self.ub.size and self.lb.size > 1 and self.ub.size > 1):
            warnings.warn(
                'The sizes of the constraint bounds are inconsistent; check the shapes of the arrays.', Warning)

    def __repr__(self):
        return '{}({}, {}, {})'.format(type(self).__name__, self.fun, self.lb, self.ub)


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
            2. An instance of LinearConstraint or NonlinearConstraint.
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
            eliminate_lin_eq: bool, optional
                Flag indicating whether the linear equality constraints should be eliminated. By default, it is True.
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

    if len(stack()) >= 4 and stack()[2][3].lower() == 'pdfo':
        # The invoker is a solver called by pdfo, then prepdfo should have been called in pdfo. In this case, we set
        # prob_info to an empty dictionary.
        prob_info = dict()
        return fun, x0, bounds, constraints, options, method, prob_info

    # Save problem information in probinfo.
    # At return, probinfo has the following fields:
    # 1. raw_data: problem data before preprocessing/validating, including fun, x0, constraints, lb, ub, options
    # 2. refined_data: problem data after preprocessing/validating, including fun, x0, constraints, lb, ub, options
    # 3. fixedx: a bool vector indicating which variables are fixed by the constraints
    # 4. fixedx_value: the values of the variables fixed by the constraints
    # 5. nofreex: whether all variables are fixed by bound constraints and/or the linear equality constraints
    # 6. infeasible_bound: a bool vector indicating which bound constraints are infeasible
    # 7. infeasible_linear: a bool vector indicating which linear constraints are infeasible (up to naive tests)
    # 8. infeasible_nonlinear: a bool vector indicating which nonlinear constraints are infeasible (up to naive tests)
    # 9. trivial_linear
    # 10. infeasible: whether the problem is infeasible (up to naive tests)
    # 11. scaled: whether the problem is scaled
    # 12. scaling_factor: vector of scaling factors
    # 13. shift: vector of shifts
    # 14. reduced: whether the problem is reduced (due to fixed variables)
    # 15. space_chg: function of the change of space (if applicable)
    # 16. bounds_in_lin_eq: indices of the bounds added in the linear inequalities (if applicable)
    # 17. raw_type: problem type before reduction
    # 18. raw_dim: problem dimension before reduction
    # 19. refined_type: problem type after reduction
    # 20. refined_dim: problem dimension after reduction
    # 21. feasibility_problem: whether the problem is a feasibility problem
    # 22. user_options_fields: the fields in the user-specified options
    # 23. warnings: warnings during the preprocessing/validation
    # 24. constr_metadata: metadata of each constraint, which is needed to build the output constr_value. It contains:
    #         - linear_indices: the indices of the linear constraints in the argument `constraints`.
    #         - nonlinear_indices: the indices of the nonlinear constraints in the argument `constraints`.
    #         - data: a list of metadata for each constraint (length, bounds, dropped indices, ...).
    prob_info = {'raw_data': {'objective': fun, 'x0': x0, 'args': args, 'bounds': bounds, 'constraints': constraints,
                              'options': options}, 'feasibility_problem': False}

    # If fun is None, then we are dealing with a feasibility problem; set fun to a fake objective function that
    # returns a constant.
    if fun is None:
        prob_info['feasibility_problem'] = True

        def fun(x_loc, *args_loc):
            return np.float64(0)

        warn_message = '{}: there is no objective function. A feasibility problem will be solved.'.format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)
    elif not callable(fun):
        raise ValueError('{}: the objective function should be callable.'.format(invoker))

    # The extra-arguments of the objective function should be given as a list or a tuple.
    if args is not None and not hasattr(args, '__len__'):
        args = [args]

    # Get the extrem barrier for the objective function.
    try:
        from .gethuge import gethuge
    except ImportError:
        import_error_so('gethuge')
    hugefun = gethuge('fun')

    # The objective function should return a floating-point number.
    def fun_c(x):
        try:
            fun_x = fun(x) if hasattr(args, '__len__') and len(args) == 0 else fun(x, *args)
        except TypeError:
            raise TypeError('{}: the number of parameters is inconsistent with `args`.'.format(invoker))

        if hasattr(fun_x, '__len__') and len(fun_x) == 1:
            fun_x = fun_x[0]
        elif (hasattr(fun_x, '__len__') or not isinstance(fun_x, scalar_types)) and fun_x is not None:
            raise ValueError('{}: the objective function should return a scalar.'.format(invoker))

        fun_x = np.float64(fun_x)
        if not np.isfinite(fun_x) or fun_x > hugefun:
            fun_x = np.float64(hugefun)

        return fun_x

    # The initial guess should be an unidimensional ndarray.
    if isinstance(x0, scalar_types):
        x0 = [x0]
    if not hasattr(x0, '__len__'):
        raise ValueError('{}: the initial guess should be a scalar or a vector.'.format(invoker))
    try:
        x0_c = np.asarray(x0, dtype=np.float64)

        if x0_c.size == 0:
            raise ValueError('{}: the initial guess should not be empty.'.format(invoker))

        # Reshape the initial guess.
        if len(x0_c.shape) > 1 and np.prod(x0_c.shape) == x0_c.size:
            x0_c = x0_c.reshape(x0_c.size)
    except ValueError:
        raise ValueError('{}: the initial guess should contain only scalars.'.format(invoker))
    if len(x0_c.shape) > 1:
        raise ValueError('{}: the initial guess should be a scalar or a vector.'.format(invoker))

    # To clarify, the length of x0 is denoted lenx0 instead of n.
    lenx0 = x0_c.size

    abnormal_x0 = np.logical_or(np.isnan(x0_c), np.isinf(x0_c))
    if abnormal_x0.any():
        x0_c[abnormal_x0] = 0
        warn_message = '{}: x0 contains NaN or infinite values; they are replaced by 0.'.format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    # If prepdfo is called by a solver (i.e. uobyqa, newuoa, bobyqa, lincoa, or cobyla), the selected method is the
    # solver itself.
    if method is None and invoker != 'pdfo':
        method = invoker
    if method is not None and not isinstance(method, str):
        raise ValueError('{}: the method name should be a string.'.format(invoker))

    # Validate the bounds and define its feasibility.
    lb, ub, infeasible_bound, fixed_indices, fixed_values = _bounds_validation(invoker, bounds, lenx0)
    prob_info['raw_data']['bounds'] = (lb, ub)
    prob_info['infeasible_bound'] = infeasible_bound
    prob_info['fixedx'] = fixed_indices
    prob_info['fixedx_value'] = fixed_values

    # Since fixedx_value may be revised in _constraints_validation, we will record it in prob_info only after that. We
    # save its current value in fixedx_value_save, which will be used when calculating the constriant violation at x0.
    fixed_values_save = fixed_values.copy()

    # Validate the linear and nonlinear constraint, and define its feasibility.
    # 1. The constraints will be reduced if some but not all variables are fixed by the bound constraints. See
    #    _bounds_validation for why we do not reduce the problem when all variables are fixed.
    # 2. The 'trivial constraints' will be excluded (if any).
    # 3. In addition, get the indices of infeasible and trivial constraints (if any) and save the information in
    # prob_info.
    free_indices = np.logical_not(fixed_indices)
    fixed_bounds = fixed_indices.copy()  # components fixed only by the bounds
    prob_info['reduced'] = any(fixed_indices) and any(free_indices)
    constraints_c, raw_constraints_c, infeasible_linear, infeasible_nonlinear, trivial, infeasible, prob_info = \
        _constraints_validation(invoker, constraints, lenx0, fixed_indices, fixed_values, x0_c, prob_info)
    prob_info['raw_data']['constraints'] = raw_constraints_c
    prob_info['infeasible_linear'] = infeasible_linear
    prob_info['infeasible_nonlinear'] = infeasible_nonlinear
    prob_info['infeasible'] = infeasible
    prob_info['trivial_linear'] = trivial
    fixed_values = prob_info['fixedx_value']  # it may be changed if the problem is infeasible

    # Problem type before reduction
    prob_info['raw_type'] = _problem_type(lb, ub, constraints_c)
    prob_info['raw_dim'] = lenx0

    # Reduce fun, x0, lb, and ub if some but not all variables are fixed by the bound constraints. A copy of the
    # original bounds are kept for sake of clarity if some variables are fixed by the bound constraints, and the rest by
    # the linear equality constraints.
    lb_old, ub_old = lb.copy(), ub.copy()
    if prob_info['reduced']:
        x0_c = x0_c[free_indices]
        lb = lb[free_indices]
        ub = ub[free_indices]
        lenx0 = x0_c.size

        def fun_c_reduced(freex_value):
            return fun_c(_fullx(freex_value, fixed_values, free_indices, fixed_indices))
    else:
        fun_c_reduced = fun_c

    # After pre-processing the linear/bound constraints, the problem may turn out infeasible, or x may turn out fixed by
    # the bounds.
    prob_info['nofreex'] = all(prob_info['fixedx'])

    # Validate and preprocess options, adopt default options if needed. This should be done after reducing the problem,
    # because BOBYQA requires rhobeg <= min(ub-lb)/2.
    # user_options_fields is a cell array that contains the names of all the user-defined options (even if the options
    # turns out invalid). It will be needed if the user does not specify a solver or specifies a wrong solver. In such a
    # scenario, we will select the solver later, and the options may have to be revised accordingly. We will raise a
    # warning when revising an option that is in user_options_fields. No warning is needed if we are dealing with an
    # option that is not in user_options_fields.
    options_c, prob_info['user_options_fields'], method = \
        _options_validation(invoker, options, method, lenx0, lb, ub, list_warnings)

    # Scale the problem if necessary and if intended, x_before_scaling = scaling_factor.*x_after_scaling + shift.
    # This should be done after revising x0, which can affect the shift.
    prob_info['scaled'] = False
    prob_info['scaling_factor'] = np.ones(lenx0)
    prob_info['shift'] = np.zeros_like(x0_c)
    if options_c['scale'] and not prob_info['nofreex'] and not prob_info['infeasible']:
        # Scale and shift the problem so that all the bounds become [-1, 1]. It is done only if all variables have both
        # lower and upper bounds.
        fun_c_reduced, x0_c, lb, ub, constraints_c, scaling_factor, shift, _ = \
            _scale_problem(fun_c_reduced, x0_c, lb, ub, constraints_c, list_warnings)

        # Scale and shift the problem so that:
        #   1. for the variables that have both lower bound and upper bound, the bounds become [-1, 1].
        #   2. the other variables will be shifted so that the corresponding component of x0 becomes 0.
        prob_info['scaled'] = True
        prob_info['scaling_factor'] = scaling_factor
        prob_info['shift'] = shift

    # If the problem contains any linear equality constraints, remove them from the constraints and force the iterates
    # to lie on the affine space generate by them as long as possible. The affine space is computed using the QR
    # factorization with pivoting (or the SVD factorization) of the Jacobian of the linear equality constraints.
    if not prob_info['nofreex'] and not prob_info['infeasible'] and options_c['eliminate_lin_eq']:
        space_chg, lb, ub, constraints_c, prob_info = _eliminate_linear_equalities(
            invoker, constraints_c, x0_c, lb, ub, prob_info, list_warnings)

        # The space dimension reduction may engender that NPT is no more in the required interval.
        options_c = _pre_npt_elimination(invoker, lb.size, prob_info['user_options_fields'], options_c, list_warnings)
    else:
        space_chg = None
    prob_info['space_chg'] = space_chg

    if prob_info['nofreex']:
        # If the variables are fixed by both the bounds and the linear equality constraints, we should revert the
        # scaling. Note that only the values fixed by the linear equality constraints are scaled, not the one fixed by
        # the bounds.
        fun_c_reduced = fun_c
        if prob_info['scaled']:
            x0_c = prob_info['scaling_factor'] * x0_c + prob_info['shift']
            lb = prob_info['scaling_factor'] * lb + prob_info['shift']
            ub = prob_info['scaling_factor'] * ub + prob_info['shift']
            fixed_lin = np.logical_not(fixed_bounds)
            prob_info['fixedx_value'][fixed_lin] *= prob_info['scaling_factor']
            prob_info['fixedx_value'][fixed_lin] += prob_info['shift']
            constraints_c = raw_constraints_c
            prob_info['scaled'] = False
            prob_info['scaling_factor'] = np.ones(lenx0)
            prob_info['shift'] = np.zeros_like(x0_c)

    # If an explicit formulation for the linear equality constraints has been found, the bound constraints become linear
    # inequality constraints, and the existing constraints are modified to fit the new searching space. The bounds
    # cannot be kept as is, since they may not be bound constraints for the modified space. Consider for example in a
    # three-dimensional space the box-constraints defined by a cube, and observe that there exists a section of it
    # (i.e., an intersection between the cube and an affine hyperplane of a single linear equality constraint) that
    # gives rise to an hexagonal slice (that is hence not a box constraint on the plane).
    if not prob_info['nofreex'] and not prob_info['infeasible'] and space_chg is not None:
        lenx0 = lb.size
        x0_c = np.zeros(lenx0, dtype=np.float64)  # the intercept contains the information about the original x0

        # The nonlinear constraint function is composed with the affine transformation.
        if constraints_c['nonlinear'] is not None:
            # The construction of the vector of variables including the fixed ones has already been made.
            ctr_fct_ori = constraints_c['nonlinear']['fun']
            constraints_c['nonlinear'] = {
                'type': 'ineq',
                'fun': lambda x_red: ctr_fct_ori(space_chg(x_red))
            }

    # Problem after reduction.
    prob_info['refined_type'] = _problem_type(lb, ub, constraints_c)
    prob_info['refined_dim'] = lenx0

    # Can the invoker handle the given problem? This should be done after the problem type has been 'refined' (for
    # example, NEWUOA can handle a bound-constrained problem if every defined bound is fixed).
    if not _prob_solv_match(prob_info['refined_type'], invoker.lower()):
        if invoker.lower() == 'pdfo':
            raise SystemError(
                '{}: UNEXPECTED ERROR: problem and solver do not match; it should not happen when the invoker is pdfo '
                'or the problem is not a structure.'.format(fun_name))
        else:
            raise ValueError(
                '{}: a {} problem received; {} cannot solve '
                'it.'.format(fun_name, prob_info['refined_type'].replace('-', ' '), invoker))

    # Revise x0 for bound and linearly constrained problems. This is necessary for LINCOA, which accepts only a feasible
    # x0. Should we do this even if there are nonlinear constraints? For now, we do not, because doing so may
    # dramatically increase the infeasibility of x0 with respect to the nonlinear constraints.
    if prob_info['refined_type'] in ['bound-constrained', 'linearly-constrained'] and not prob_info['nofreex'] and \
            not prob_info['infeasible']:
        x0_ori = x0_c.copy()
        result = _project(x0_c, lb, ub, constraints_c)
        x0_c = result.x

        if not prob_info['feasibility_problem'] and \
                np.linalg.norm(x0_ori - x0_c) > eps * max(1.0, np.linalg.norm(x0_ori)):
            warn_message = '{}: x0 is revised to satisfy the constraints.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    if prob_info['infeasible']:
        # The problem turned infeasible when analyzing the linear equality constraints. The original x0 should be
        # restored to include the reduced variables, together with the constraints.
        prob_info['constrv_x0'], prob_info['nlc_x0'] = \
            _constr_violation(invoker, x0_c, lb, ub, constraints_c, prob_info)

        # The constraint violation calculated by _constr_violation does not include the violation of x0 for the bounds
        # corresponding to fixedx; the corresponding values of x0 are in fixedx_value, while the values of the bounds
        # (lb and ub are the same up to eps) are in fixedx_value_save. Thus the violation is
        # abs(fixedx_value-fixedx_value_save).
        max_fixed = np.nanmax(np.abs(fixed_values - fixed_values_save)) if fixed_values.size > 0 else 0
        prob_info['constrv_x0'] = max(prob_info['constrv_x0'], max_fixed)
    elif prob_info['nofreex']:
        prob_info['constrv_fixedx'], prob_info['nlc_fixedx'] = \
            _constr_violation(invoker, prob_info['fixedx_value'], lb_old, ub_old, raw_constraints_c, prob_info)

    # If is possible that both prob_info['reduced'] and prob_info['nofreex'] are True, if the bound constraint fixed
    # some (but not all) constraints and the linear equality constraint fixed the others.
    if space_chg is not None and not prob_info['nofreex']:
        def fun_c_space(freex_value):
            return fun_c_reduced(space_chg(freex_value))
    else:
        fun_c_space = fun_c_reduced  # the variable vector is not reduced

    # Select a solver if the invoker is 'pdfo' and no one is provided.
    if invoker.lower() == 'pdfo':
        # lb and ub will be used for defining rhobeg.
        prob_info['refined_data'] = {'lb': lb, 'ub': ub}

        method = _solver_selection(invoker, method, options_c, prob_info, list_warnings)

    if method.lower() == 'bobyqa' and not prob_info['nofreex'] and not prob_info['infeasible'] and \
            not prob_info['feasibility_problem']:
        # The Fortran code of BOBYQA will revise x0 so that the distance between x0 and the inactive bounds is at least
        # rhobeg. We do it here in order to raise a warning when such a revision occurs. After this, the Fortran code
        # will not revise x0 again. If the options['honour_x0'] = True, then we keep x0 unchanged and revise rhobeg if
        # necessary.
        x0_c, options_c = \
            _pre_rhobeg_x0(invoker, x0_c, lb, ub, prob_info['user_options_fields'], options_c, list_warnings)

    # Store the warnings raised during the pre-processing.
    prob_info['warnings'] = list_warnings

    # The refined data can be useful when debugging. It will be used in postpdfo even if the debug mode is not enabled.
    prob_info['refined_data'] = \
        {'objective': fun_c_space, 'x0': x0_c, 'lb': lb, 'ub': ub, 'constraints': constraints_c, 'options': options_c}

    # When the problem is a linear feasibility problem, PDFO will return the current x0, which has been revised by
    # project. The constraint violation at x0 is needed to set the output. Note that there is no nonlinear constraint in
    # this case.
    if prob_info['feasibility_problem'] and prob_info['refined_type'] != 'nonlinearly-constrained':
        prob_info['constrv_x0'], prob_info['nlc_x0'] = \
            _constr_violation(invoker, x0_c, lb, ub, constraints_c, prob_info)

    if not options_c['debug']:
        # Do not carry the raw data with us unless in debug mode.
        del prob_info['raw_data']

        # The refined data is used only when the problem is scaled. It can also be useful when debugging.
        # Set this field to empty instead of remove it, because postpdfo require that this field exists.
        if not prob_info['scaled']:
            prob_info['refined_data'] = {}

    return fun_c_space, x0_c, {'lb': lb, 'ub': ub}, constraints_c, options_c, method, prob_info


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
                # If the bounds are defined with the SciPy Bounds class, they are converted to the local Bounds class so
                # that some required pre-processing are done.
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

        # Extract lower and upper bounds from the bounds arguments.
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

        # If both lb and ub are empty, the problem should be considered as unbounded.
        if lb.size == 0 and ub.size == 0:
            lb = np.full(lenx0, -np.inf, dtype=np.float64)
            ub = np.full(lenx0, np.inf, dtype=np.float64)

        if lb.size != lenx0 or ub.size != lenx0:
            raise ValueError('{}: the bound size should be equal to the number of variables.'.format(invoker))
    else:
        lb = np.full(lenx0, -np.inf, dtype=np.float64)
        ub = np.full(lenx0, np.inf, dtype=np.float64)

    # NaN bounds are set to infinite values not to take them into account. The conversion is not made in the class
    # Bounds since the user may define the bounds as a list of 2-tuples, in which case the class is never called.
    lb[np.isnan(lb)] = -np.inf
    ub[np.isnan(ub)] = np.inf

    # Check the infeasibility of the bounds.
    infeasible_lb, lb_mod = np.logical_and(lb > 0, np.isinf(lb)), lb.copy()
    infeasible_ub, ub_mod = np.logical_and(ub < 0, np.isinf(ub)), ub.copy()
    lb_mod[infeasible_lb] = -np.inf
    ub_mod[infeasible_ub] = np.inf
    infeasible = np.logical_or(infeasible_lb, np.logical_or(infeasible_ub, lb > ub))
    fixed_indices = (np.abs(lb_mod - ub_mod) < 2 * eps)
    fixed_values = (lb_mod[fixed_indices] + ub_mod[fixed_indices]) / 2

    return lb, ub, infeasible, fixed_indices, fixed_values


def _constraints_validation(invoker, constraints, lenx0, fixed_indices, fixed_values, x0, prob_info):
    """Validation and pre-processing of the constraints.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    constraints: dict, LinearConstraint, NonlinearConstraint or list of them
        The same as in prepdfo.
    lenx0: integer
        The size of the problem.
    fixed_indices: ndarray, shape (n,)
        The boolean vector indicating whether the variable is fixed or not.
    fixed_values: ndarray, shape (n,)
        The values of the fixed variables.
    x0: ndarray, shape (n,)
        The same as in prepdfo.
    prob_info: dict
        The problem information.

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
            linear_constraint_types = LinearConstraint
            nonlinear_constraint_types = NonlinearConstraint

        # convert the constraints as a list
        if isinstance(constraints, dict) or not (hasattr(constraints, '__len__')):
            is_list = False
            constraints_c = [constraints]
        else:
            is_list = True
            constraints_c = constraints

        # create the linear/nonlinear sub-lists
        list_linear = []
        list_nonlinear = []
        for constraint in constraints_c:
            if isinstance(constraint, linear_constraint_types):
                # If the user provided a linear constraint through the LinearConstraint class of SciPy, it is converted
                # to the local LinearConstraint so that some pre-processing are done on the matrices.
                list_linear.append(LinearConstraint(a=constraint.A, lb=constraint.lb, ub=constraint.ub))
            elif isinstance(constraint, nonlinear_constraint_types) or \
                    (isinstance(constraint, dict) and {'type', 'fun'} <= set(constraint.keys()) and
                     isinstance(constraint['type'], str) and constraint['type'] in ['eq', 'ineq'] and
                     (callable(constraint['fun']) or constraint['fun'] is None)):
                if isinstance(constraint, nonlinear_constraint_types):
                    # If the user provided a nonlinear constraint through the NonlinearConstraint class of SciPy, it is
                    # converted to the local NonlinearConstraint so that some pre-processing are done on the vectors and
                    # the constraint function.
                    list_nonlinear.append(NonlinearConstraint(fun=constraint.fun, lb=constraint.lb, ub=constraint.ub))
                else:
                    list_nonlinear.append(constraint)
            else:
                # The constraint is neither linear nor nonlinear.
                raise ValueError(
                    "{}: each constraint should be an instance of the `LinearConstraint` or `NonlinearConstraint` "
                    "class, or a dictionary with fields 'type' and 'fun'.".format(invoker))

        # Create the constraint metadata so that the list of constraint values can be constructed in the post-processing.
        linear_constr_indices = [i for (i, con) in enumerate(constraints_c) if isinstance(con, linear_constraint_types)]
        nonlinear_constr_indices = list(set(range(len(constraints_c))) - set(linear_constr_indices))
        prob_info['constr_meta'] = {
            'is_list': is_list,
            'linear_indices': linear_constr_indices,
            'nonlinear_indices': nonlinear_constr_indices,
            'data': [None] * len(constraints_c)
        }

        # If no bounds are provided with some nonlinear constraints, they should not be considered. Moreover, the bounds
        # in the nonlinear constraints should be feasible.
        reduced_list_nonlinear = []
        infeasible_nonlinear = np.asarray([], dtype=bool)

        # The following registered the state of every nonlinear constraints, i.e., whether it should be considered or
        # not.
        list_nonlinear_bound_types = []
        for nonlinear_constraint, i_meta in zip(list_nonlinear, prob_info['constr_meta']['nonlinear_indices']):
            if isinstance(nonlinear_constraint, nonlinear_constraint_types):
                lb_nonlinear, ub_nonlinear = nonlinear_constraint.lb, nonlinear_constraint.ub
                if (lb_nonlinear.size > 0 and not np.logical_and(np.isinf(lb_nonlinear), lb_nonlinear < 0).all()) or \
                        (ub_nonlinear.size > 0 and not np.logical_and(np.isinf(ub_nonlinear), ub_nonlinear > 0).all()):
                    reduced_list_nonlinear.append(nonlinear_constraint)
                    lbx = np.logical_not(np.logical_and(np.isinf(lb_nonlinear), lb_nonlinear < 0))
                    ubx = np.logical_not(np.logical_and(np.isinf(ub_nonlinear), ub_nonlinear > 0))
                    lb_free = np.all(np.logical_not(lbx))
                    ub_free = np.all(np.logical_not(ubx))
                    list_nonlinear_bound_types.append({
                        'lbx': lbx,
                        'lb_free': lb_free,
                        'ubx': ubx,
                        'ub_free': ub_free
                    })

                    # Set the metadata of the nonlinear constraints that are now accessible.
                    prob_info['constr_meta']['data'][i_meta] = {
                        'trivial': False,
                        'len': max(lb_nonlinear.size, ub_nonlinear.size),
                        'dropped_indices_lb': np.logical_not(lbx),
                        'dropped_indices_ub': np.logical_not(ubx),
                        'lb': lb_nonlinear,
                        'ub': ub_nonlinear
                    }

                    infeasible_lb = np.logical_and(lb_nonlinear > 0, np.isinf(lb_nonlinear))
                    infeasible_ub = np.logical_and(ub_nonlinear < 0, np.isinf(ub_nonlinear))
                    infeasible_nonlinear = \
                        np.r_[infeasible_nonlinear,
                              np.logical_or(infeasible_lb, np.logical_or(infeasible_ub, lb_nonlinear > ub_nonlinear))]
                else:
                    # The nonlinear constraint, defined as a NonlinearConstraint structure, is considered trivial. It
                    # means that the lower bound of each component is -inf and the upper bound of each component is
                    # +inf.
                    prob_info['constr_meta']['data'][i_meta] = {
                        'trivial': True,
                        'len': max(lb_nonlinear.size, ub_nonlinear.size)
                    }
            else:
                # The nonlinear constraint is defined as a dictionary, for which infeasibility cannot be checked
                reduced_list_nonlinear.append(nonlinear_constraint)
                list_nonlinear_bound_types.append(None)

                # Set the metadata of the nonlinear constraints. Since it is defined as a dictionary, we do not have
                # access to any meta-information. However, the structure is constructed with default values, so that it
                # can be modified during the execution of PDFO. In fact, since the dictionary are mutable, the values
                # archived i prob_info['constr_meta']['data'] can be modified on the fly of the execution of PDFO.
                prob_info['constr_meta']['data'][i_meta] = {
                    'trivial': False,
                    'len': -1,
                    'dropped_indices_lb': None,
                    'dropped_indices_ub': None,
                    'lb': None,
                    'ub': None
                }

        list_nonlinear = reduced_list_nonlinear

        # Build the linear constraints.
        a_linear = np.asarray([[]], dtype=np.float64, order='F').reshape(0, lenx0)
        lb_linear = np.asarray([], dtype=np.float64)
        ub_linear = np.asarray([], dtype=np.float64)

        for linear_constraint, i_meta in zip(list_linear, prob_info['constr_meta']['linear_indices']):
            # The type of linear_constraint is necessarily LinearConstraint.
            a_local = linear_constraint.A
            if a_local.size == 0:
                a_local = a_local.reshape(0, lenx0)
            if a_local.shape[1] != lenx0:
                raise ValueError(
                    '{}: the number of columns in A is inconsistent with the number of variables'.format(invoker))

            # If no bounds are provided (either lower or upper), an infinite one is defined.
            if linear_constraint.lb.size != 0:
                lb_local = linear_constraint.lb
            else:
                lb_local = np.full(a_local.shape[0], -np.inf)
            if linear_constraint.ub.size != 0:
                ub_local = linear_constraint.ub
            else:
                ub_local = np.full(a_local.shape[0], np.inf)

            # Set the metadata related to this constraint. Since it a linear constraint, all we need is the complete
            # left-hand side to perform in the post-processing the product Ax.
            prob_info['constr_meta']['data'][i_meta] = {'trivial': False, 'A': a_local}

            # Add the current linear constraint to the global one.
            if a_local.size > 0:
                # If an empty matrix A is provided, it would be reshaped to (0, 1) and hence, the vectors lb_local and
                # ub_local would be set respectively to [-inf] and [+inf]. Thus, it would create an incompatibility if
                # we try to add them to lb_linear and ub_linear.
                a_linear = np.concatenate((a_linear, a_local), axis=0)
                lb_linear = np.r_[lb_linear, lb_local]
                ub_linear = np.r_[ub_linear, ub_local]

        # Remove the abnormal constraints and check infeasibility.
        if prob_info['reduced']:
            free_indices = np.logical_not(fixed_indices)
            a_reduced = a_linear[:, free_indices]
            a_fixed = np.dot(a_linear[:, fixed_indices], fixed_values)
            lb_reduced = lb_linear - a_fixed
            ub_reduced = ub_linear - a_fixed
        else:
            a_reduced = a_linear.copy()
            lb_reduced = lb_linear.copy()
            ub_reduced = ub_linear.copy()

        if a_reduced.size == 0:
            trivial = np.asarray([], dtype=bool)
            infeasible_linear = np.asarray([], dtype=bool)
        else:
            row_norm_inf = np.nanmax(np.abs(a_reduced), 1)
            zero_rows = (row_norm_inf == 0)
            infeasible_zero = np.logical_or(np.logical_and(zero_rows, lb_reduced > 0),
                                            np.logical_and(zero_rows, ub_reduced < 0))
            trivial_zero = np.logical_and(np.logical_and(zero_rows, lb_reduced <= 0),
                                          np.logical_and(zero_rows, ub_reduced >= 0))
            row_norm_inf[zero_rows] = 1.0
            lb_linear_norm = lb_reduced / row_norm_inf
            ub_linear_norm = ub_reduced / row_norm_inf
            lb_ub_and = np.logical_and(np.logical_and(np.isinf(lb_linear_norm), lb_linear_norm < 0),
                                       np.logical_and(np.isinf(ub_linear_norm), ub_linear_norm > 0))
            lb_ub_or = np.logical_or(np.logical_and(np.isinf(lb_linear_norm), lb_linear_norm > 0),
                                     np.logical_and(np.isinf(ub_linear_norm), ub_linear_norm < 0))
            infeasible_linear = np.logical_or(infeasible_zero, np.logical_or(lb_ub_or, lb_reduced > ub_reduced))
            trivial = np.logical_or(trivial_zero, lb_ub_and)

        # Check the infeasibility of the problem
        infeasible = any(np.r_[prob_info['infeasible_bound'], infeasible_linear, infeasible_nonlinear])
        if infeasible:
            fixed_values = x0[fixed_indices]
            prob_info['fixedx_value'] = fixed_values
            a_fixed = np.dot(a_linear[:, fixed_indices], fixed_values)
            lb_reduced = lb_linear - a_fixed
            ub_reduced = ub_linear - a_fixed

        # Reduce the linear constraints according to the fixed variables
        if lb_linear.size != 0 or ub_linear.size != 0:
            a_reduced = a_reduced[np.logical_not(trivial), :]
            lb_reduced = lb_reduced[np.logical_not(trivial)]
            ub_reduced = ub_reduced[np.logical_not(trivial)]

        # Build the nonlinear constraints.
        try:
            from .gethuge import gethuge
        except ImportError:
            # If gethuge cannot be imported, the package is most likely not built.
            import_error_so('gethuge')
        hugecon = gethuge('con')

        # Define the indices of the nonlinear constraints that are not trivial
        non_trivial = [i for i in range(len(constraints_c)) if not prob_info['constr_meta']['data'][i]['trivial']]
        non_linear_non_trivial_indices = [i for i in prob_info['constr_meta']['nonlinear_indices'] if i in non_trivial]

        # Define the global constraint function.
        def fun_nonlinear(x, raw=False):
            fun_x = np.asarray([], dtype=np.float64)

            for nlc_constraint, b_type, i_meta in zip(list_nonlinear, list_nonlinear_bound_types,
                                                      non_linear_non_trivial_indices):
                # Get the value of the constraint function.
                if not raw and prob_info['reduced']:
                    x_full = _fullx(x, fixed_values, free_indices, fixed_indices)
                else:
                    x_full = x
                if isinstance(nlc_constraint, nonlinear_constraint_types):
                    constraint_x = nlc_constraint.fun(x_full)
                elif nlc_constraint['fun'] is not None:
                    constraint_x = nlc_constraint['fun'](x_full)
                else:
                    constraint_x = np.asarray([], dtype=np.float64)

                if constraint_x is None:
                    # If the constraint function returned anything, we convert the default None value to NaN, which can
                    # be understood by Fortran.
                    constraint_x = [np.nan]
                elif isinstance(constraint_x, scalar_types):
                    constraint_x = [constraint_x]

                if not hasattr(constraint_x, '__len__'):
                    raise ValueError('{}: the constraint function should return a vector or a scalar.'.format(invoker))

                constraint_x = np.asarray(constraint_x, dtype=np.float64)

                # Use extreme barrier to cope with the 'hidden constraints'.
                constraint_x[np.logical_or(np.isnan(constraint_x), constraint_x > hugecon)] = hugecon

                # This part is NOT extreme barrier. We replace extremely negative values of cineq (which leads to no
                # constraint violation) by -hugecon. Otherwise, NaN or Inf may occur in the interpolation models.
                constraint_x[constraint_x < -hugecon] = -hugecon

                if len(constraint_x.shape) != 1:
                    raise ValueError('{}: the constraint function should return a vector or a scalar.'.format(invoker))

                # Set the metadata related to this constraint if it has not been done yet.
                if prob_info['constr_meta']['data'][i_meta]['len'] < 0:
                    prob_info['constr_meta']['data'][i_meta]['len'] = constraint_x.size
                    prob_info['constr_meta']['data'][i_meta]['dropped_indices_lb'] = np.full(constraint_x.shape, False)
                    prob_info['constr_meta']['data'][i_meta]['lb'] = np.zeros_like(constraint_x)

                    # if the metadata has not been set, nlc_constraint is necessarily defined as a dictionary.
                    prob_info['constr_meta']['data'][i_meta]['dropped_indices_ub'] = \
                        np.full(constraint_x.shape, not nlc_constraint['type'] == 'eq')
                    if nlc_constraint['type'] == 'eq':
                        prob_info['constr_meta']['data'][i_meta]['ub'] = np.zeros_like(constraint_x)
                    else:
                        prob_info['constr_meta']['data'][i_meta]['ub'] = np.full_like(constraint_x, np.inf)

                lenm = constraint_x.size
                if isinstance(nlc_constraint, nonlinear_constraint_types) and not infeasible:
                    if nlc_constraint.lb.size not in [0, lenm] or \
                            nlc_constraint.ub.size not in [0, lenm] or \
                            (nlc_constraint.lb.size == 0 and nlc_constraint.ub.size == 0):
                        raise ValueError(
                            '{}: the size of the vector returned by the constraint function is inconsistent with the '
                            'constraint bounds; check the shapes of the arrays.'.format(invoker))

                    # Convert the constraints defined as lb <= c(x) <= ub into c_extended(x) <= 0.
                    constraint_x_tmp = np.array([], dtype=np.float64)
                    if not b_type['lb_free']:
                        lbx = b_type['lbx']
                        constraint_x_tmp = np.r_[constraint_x_tmp, nlc_constraint.lb[lbx] - constraint_x[lbx]]
                    if not b_type['ub_free']:
                        ubx = b_type['ubx']
                        constraint_x_tmp = np.r_[constraint_x_tmp, constraint_x[ubx] - nlc_constraint.ub[ubx]]
                    constraint_x = constraint_x_tmp
                elif isinstance(nlc_constraint, dict) and nlc_constraint['type'] == 'eq' and not infeasible:
                    # Necessarily, nlc_constraint is defined as a dictionary, for which all the constraints have to be
                    # considered.
                    constraint_x = np.r_[-constraint_x, constraint_x]
                elif not infeasible:
                    # nlc_constraint is defined as a dictionary, for which all the constraints have to be considered.
                    # Moreover, it consists of an inequality constraint c(x) >= 0, which has to be inversed,
                    constraint_x *= -1
                else:
                    # The problem is infeasible, only the constraint evaluation should be stored.
                    pass

                # Add the current nonlinear constraint evaluation to the general one.
                fun_x = np.r_[fun_x, constraint_x]

            return fun_x

        # Define the global linear and nonlinear constraints.
        linear_constraints = None if len(list_linear) == 0 else LinearConstraint(a_reduced, lb_reduced, ub_reduced)
        nonlinear_constraints = None if len(list_nonlinear) == 0 else {'type': 'ineq', 'fun': fun_nonlinear}
        raw_linear_constraints = None if len(list_linear) == 0 else LinearConstraint(a_linear, lb_linear, ub_linear)
        raw_nonlinear_constraints = None if len(list_nonlinear) == 0 else {'type': 'ineq',
                                                                           'fun': lambda x: fun_nonlinear(x, raw=True)}
    else:
        # No constraints have been provided.
        linear_constraints = None
        nonlinear_constraints = None
        raw_linear_constraints = None
        raw_nonlinear_constraints = None
        infeasible_linear = np.asarray([], dtype=bool)
        infeasible_nonlinear = np.asarray([], dtype=bool)
        trivial = np.asarray([], dtype=bool)
        infeasible = prob_info['infeasible_bound']

    return {'linear': linear_constraints, 'nonlinear': nonlinear_constraints}, \
           {'linear': raw_linear_constraints, 'nonlinear': raw_nonlinear_constraints}, infeasible_linear, \
           infeasible_nonlinear, trivial, infeasible, prob_info


def _eliminate_linear_equalities(invoker, constraints, x0, lb, ub, prob_info, list_warnings):
    """Computation of the affine hyperplane on which the iterates should lie if any. It is done using a QR
    factorization with pivoting of the Jacobian of the linear equality constraints if SciPy is installed, and using its
    SVD otherwise.

    Parameters
    ----------
    invoker: str
        The name of the invoker.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: NonlinearConstraint
                The nonlinear constraints of the problem.
    x0: ndarray, shape (n,)
        The same as in prepdfo.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    prob_info: dict
        The problem information.
    list_warnings: list
        The same as in prepdfo.

    Returns
    -------
    The modified `prob_info`, and
    space_chg: callable
        The function that perform the change of space
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    prob_info: dict
        The problem information.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())):
        raise ValueError('{}: UNEXPECTED ERROR: the constraint should be defined as internal type.'.format(invoker))

    if constraints['linear'] is not None and not prob_info['nofreex'] and not prob_info['infeasible']:
        # Get a boolean array indicating whether the constraints are equality ones or not. We consider that a constraint
        # is an equality one if the relative difference between the lower bound and the upper bound of the constraint is
        # at most 10 times the machine epsilon.
        # Note: it is important to keep the inequality strict! If any bound is infinite, the corresponding constraint
        # would otherwise be considered as equality (because np.inf <= np.inf).
        lin_eq_indices = np.less(
            np.abs(constraints['linear'].ub - constraints['linear'].lb),
            10 * eps * np.maximum(1, np.abs(constraints['linear'].lb) + np.abs(constraints['linear'].ub)))
        if np.any(lin_eq_indices):
            # Some linear equality constraints have been detected, and their indices are stored in lin_eq_indices. We
            # first need to find a point feasible for the linear equality constraints and to store it in intercept.
            # Since x0 has previously been projected onto the bound and linear constraints, it is likely to satisfy, but
            # if not, attempt to project it.
            lin_ineq_indices = np.logical_not(lin_eq_indices)
            a_eq = constraints['linear'].A[lin_eq_indices]
            col_norm = np.sum(np.square(a_eq), axis=0)
            zero_col = np.less(col_norm, 10 * eps * np.sum(col_norm))
            p_zeros = np.arange(x0.size)
            p_zeros = np.r_[p_zeros[np.logical_not(zero_col)], p_zeros[zero_col]]
            p_zeros_inv = np.argsort(p_zeros)
            n_red = x0.size - sum(zero_col)
            a_eq = a_eq[:, np.logical_not(zero_col)]
            b_eq = (constraints['linear'].lb[lin_eq_indices] + constraints['linear'].ub[lin_eq_indices]) / 2
            x0_old = x0.copy()

            try:
                from scipy.linalg import qr

                # Scipy is installed. Compute the QR factorization with pivoting of the Jacobian of the linear equality
                # constraints, and deduced from it the consistency of the system and the reduced form of the linear
                # equality constraints.
                q_eq, r_eq, p_eq = qr(a_eq, pivoting=True)

                # The rank of the Jacobian of the linear equality constraints is determined with a relative error of 10
                # times the machine epsilon.
                rcn = np.cumsum(np.flip(np.sum(np.square(r_eq), axis=1)))
                rank_a_eq = sum(np.greater_equal(rcn, 10 * eps * rcn[-1]))

                # Compute the coefficients of the hyperplane associated with the linear equality constraints.
                qtb = np.dot(q_eq.T, b_eq)
                riv = np.linalg.inv(r_eq[:rank_a_eq, :rank_a_eq])  # r_eq[:rank, :rank] is invertible and assumed small
                comp = np.abs(qtb[rank_a_eq:])
                if any(np.greater(comp, 10 * eps * np.max(np.abs(qtb), initial=1))):
                    # The linear equality constraints are infeasible
                    space_chg = None
                    intercept = None
                    prob_info['infeasible'] = True
                else:
                    # Compute the reduced form.
                    null_basis = -np.dot(riv, r_eq[:rank_a_eq, rank_a_eq:])
                    p_inv = np.argsort(p_eq)  # equivalent to the transposition of the matrix of permutation
                    x0 = x0[p_zeros]
                    x0[:n_red] = x0[p_eq]
                    feasible = np.r_[np.dot(riv, qtb[:rank_a_eq]) + np.dot(null_basis, x0[rank_a_eq:n_red]),
                                     x0[rank_a_eq:]]
                    intercept = np.r_[feasible[p_inv], feasible[n_red:]]
                    intercept = intercept[p_zeros_inv]

                    if rank_a_eq < x0.size:
                        # Compute the reduced form of the bounds.
                        lb = lb[p_zeros]
                        ub = ub[p_zeros]
                        lb[:n_red] = lb[p_eq]
                        ub[:n_red] = ub[p_eq]
                        p_restore = np.arange(x0.size)
                        p_restore[:n_red] = p_restore[p_inv]
                        p_restore = p_restore[p_zeros_inv]

                        # Compute the reduced form of the linear constraints.
                        A_perm = constraints['linear'].A[lin_ineq_indices, :]
                        A_perm = A_perm[:, p_zeros]
                        A_perm[:, :n_red] = A_perm[:, p_eq]
                        constraints['linear'].A = np.c_[
                            np.dot(A_perm[:, :rank_a_eq], null_basis) + A_perm[:, rank_a_eq:n_red], A_perm[:, n_red:]]
                        bound_shift = np.dot(A_perm, feasible)
                        constraints['linear'].lb = constraints['linear'].lb[lin_ineq_indices] - bound_shift
                        constraints['linear'].ub = constraints['linear'].ub[lin_ineq_indices] - bound_shift
                        row_inf = np.logical_and(
                            np.logical_and(np.isinf(lb[:rank_a_eq]), lb[:rank_a_eq] < 0),
                            np.logical_and(np.isinf(ub[:rank_a_eq]), ub[:rank_a_eq] > 0))
                        row_finite = np.logical_not(row_inf)
                        constraints['linear'].A = np.r_[
                            constraints['linear'].A,
                            np.c_[null_basis[row_finite, :], np.zeros((sum(row_finite), x0.size - n_red))]]
                        row_finite_bounds = np.r_[row_finite, np.full(x0.size - rank_a_eq, False)]
                        prob_info['bounds_in_lin_eq'] = np.where(np.logical_and(p_restore < rank_a_eq,
                                                                                row_finite_bounds[p_restore]))[0]
                        constraints['linear'].lb = np.r_[constraints['linear'].lb,
                                                         lb[row_finite_bounds] - feasible[row_finite_bounds]]
                        constraints['linear'].ub = np.r_[constraints['linear'].ub,
                                                         ub[row_finite_bounds] - feasible[row_finite_bounds]]

                        # The reduced bounds consist only of the last rows.
                        lb = lb[rank_a_eq:] - x0[rank_a_eq:]
                        ub = ub[rank_a_eq:] - x0[rank_a_eq:]

                    def space_chg(y):
                        x = np.r_[np.dot(null_basis, y[:n_red-rank_a_eq]), y] + feasible
                        x[:n_red] = x[p_inv]
                        x = x[p_zeros_inv]

                        return x

            except ImportError:
                # SciPy is not installed so that the SVD factorization is used instead.
                u_eq, s_eq, vh_eq = np.linalg.svd(a_eq)
                scn = np.cumsum(np.flip(s_eq))
                rank_a_eq = sum(np.greater_equal(scn, 10 * eps * scn[-1]))
                utb = np.dot(b_eq, u_eq)
                comp = np.abs(utb[rank_a_eq:])
                if any(np.greater(comp, 10 * eps * np.max(np.abs(utb), initial=1))):
                    # The linear equality constraints are infeasible.
                    space_chg = None
                    intercept = None
                    prob_info['infeasible'] = True
                else:
                    # Compute the reduced form.
                    x0 = x0[p_zeros]
                    feasible = np.dot(np.multiply(np.divide(1, s_eq[:rank_a_eq]), utb[:rank_a_eq]),
                                      vh_eq[:rank_a_eq, :])
                    feasible = np.r_[
                        feasible + np.dot(np.dot(vh_eq[rank_a_eq:, :], x0[:n_red]), vh_eq[rank_a_eq:, :]), x0[n_red:]]
                    intercept = feasible[p_zeros_inv]

                    if rank_a_eq < x0.size:
                        # Compute the reduced form of the bounds.
                        lb = lb[p_zeros]
                        ub = ub[p_zeros]

                        # Compute the reduced form of the linear constraints.
                        A_perm = constraints['linear'].A[lin_ineq_indices, :]
                        A_perm = A_perm[:, p_zeros]
                        constraints['linear'].A = np.c_[np.dot(A_perm[:, :n_red], vh_eq[rank_a_eq:, :].T),
                                                        A_perm[:, n_red:]]
                        bound_shift = np.dot(A_perm, feasible)
                        constraints['linear'].lb = constraints['linear'].lb[lin_ineq_indices] - bound_shift
                        constraints['linear'].ub = constraints['linear'].ub[lin_ineq_indices] - bound_shift
                        row_inf = np.logical_and(
                            np.logical_and(np.isinf(lb[:n_red]), lb[:n_red] < 0),
                            np.logical_and(np.isinf(ub[:n_red]), ub[:n_red] > 0))
                        row_finite = np.logical_not(row_inf)
                        constraints['linear'].A = np.r_[
                            constraints['linear'].A,
                            np.c_[vh_eq[rank_a_eq:, row_finite].T, np.zeros((sum(row_finite), x0.size - n_red))]]
                        row_finite_bounds = np.r_[row_finite, np.full(x0.size - n_red, False)]
                        prob_info['bounds_in_lin_eq'] = np.where(np.logical_and(p_zeros_inv < n_red,
                                                                                row_finite_bounds[p_zeros_inv]))[0]
                        constraints['linear'].lb = np.r_[constraints['linear'].lb,
                                                         lb[row_finite_bounds] - feasible[row_finite_bounds]]
                        constraints['linear'].ub = np.r_[constraints['linear'].ub,
                                                         ub[row_finite_bounds] - feasible[row_finite_bounds]]

                        # The reduced bounds consist only of the last rows.
                        lb = np.r_[np.full(n_red - rank_a_eq, -np.inf), lb[n_red:] - x0[n_red:]]
                        ub = np.r_[np.full(n_red - rank_a_eq, np.inf), ub[n_red:] - x0[n_red:]]

                    def space_chg(y):
                        x = np.r_[np.dot(y[:n_red-rank_a_eq], vh_eq[rank_a_eq:, :]), y[n_red-rank_a_eq:]] + feasible
                        x = x[p_zeros_inv]

                        return x

            if space_chg is not None and not prob_info['infeasible'] and rank_a_eq < x0.size:
                # The process may have created rows of zero in the Jacobian matrix of the linear inequality constraints.
                # They should be removed.
                if constraints['linear'].A.size > 0:
                    row_norm_inf = np.nanmax(np.abs(constraints['linear'].A), 1)
                    if np.logical_or(constraints['linear'].lb[row_norm_inf == 0] > 0,
                                     constraints['linear'].ub[row_norm_inf == 0] < 0).any():
                        raise ValueError('{}: the linear inequalities are inconsistent.'.format(invoker))
                    constraints['linear'].A = constraints['linear'].A[row_norm_inf > 0, :]
                    constraints['linear'].lb = constraints['linear'].lb[row_norm_inf > 0]
                    constraints['linear'].ub = constraints['linear'].ub[row_norm_inf > 0]
                    if constraints['linear'].A.size == 0:
                        constraints['linear'] = None
                else:
                    constraints['linear'] = None

            if intercept is not None and np.linalg.norm(x0_old - intercept) > 10 * eps * max(1.0, np.linalg.norm(x0)):
                warn_message = '{}: x0 is revised to satisfy the linear equality constraints.'.format(invoker)
                warnings.warn(warn_message, Warning)
                list_warnings.append(warn_message)

            if not prob_info['infeasible'] and rank_a_eq == x0.size:
                # The coefficient of the vector of variables are fixed by the linear constraints. Set the values of
                # coefficients that are not fixed by the bounds. Note that we choose prob_info['fixedx'].size
                # for the size of the problem and not lenx0, because some values may have been removed earlier
                # because of the bounds.
                space_chg = None
                prob_info['nofreex'] = True
                fixed_values = np.empty(prob_info['fixedx'].size, dtype=np.float64)
                fixed_values[prob_info['fixedx']] = prob_info['fixedx_value']
                fixed_values[np.logical_not(prob_info['fixedx'])] = intercept
                prob_info['fixedx_value'] = fixed_values
                np.ndarray.fill(prob_info['fixedx'], True)
        else:
            # No linear equality constraint has been found.
            space_chg = None
    else:
        # The problem is either infeasible, fixed by the bounds or admits no linear constraints.
        space_chg = None

    return space_chg, lb, ub, constraints, prob_info


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
    The preprocessed `options`, user_options_fields and `method`.

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

    options = dict() if options is None else dict(options)
    fun_name = stack()[0][3]  # name of the current function

    if invoker not in invoker_list:
        raise SystemError('{}: {} serves only {}'.format(fun_name, fun_name, ', '.join(invoker_list)))

    # Which options' fields are provided?
    if options is not None:
        option_fields = list(options.keys())
    else:
        option_fields = []
    user_option_fields = list(option_fields)

    # Default values for each options.
    maxfev = 500 * lenx0
    rhobeg = 1  # The default rhobeg and rhoend will be revised for BOBYQA
    rhoend = 1e-6
    ftarget = -np.inf
    classical = False  # call the classical Powell code?
    eliminate_lin_eq = True
    scale = False  # scale the problem according to bounds?
    honour_x0 = False  # Respect the user-defined x0? Needed by BOBYQA
    quiet = True
    debugflag = False
    chkfunval = False

    # DO NOT REMOVE THE FOLLOWING!! Scale only if all variables are with finite lower and upper bounds.
    scale = scale and np.all(np.logical_not(np.isinf(np.r_[lb, ub])))

    # What is the solver? We need this information to decide which fields are 'known' (e.g., expected), and also to set
    # npt and rhobeg, rhoend.
    if invoker == 'pdfo':
        if method is not None and method.lower() not in invoker_list:
            warn_message = '{}: unknown solver specified; {} will select one automatically.'.format(invoker, invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
            method = None
        elif method is not None:
            method = method.lower()
    else:  # invoker is in {'uobyqa', ..., 'cobyla'}
        if method is not None and method.lower() != invoker:
            warn_message = '{}: a solver different from {} is specified; it is ignored.'.format(invoker, invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        method = invoker

    # Check whether the used provided any unknown option.
    known_field = ['maxfev', 'rhobeg', 'rhoend', 'ftarget', 'classical', 'eliminate_lin_eq', 'quiet', 'debug',
                   'chkfunval']
    if method is None or method.lower() in ['bobyqa', 'lincoa', 'newuoa']:
        known_field.append('npt')
    if method is None or method.lower() in ['bobyqa', 'cobyla', 'lincoa']:
        known_field.append('scale')
    if method is None or method.lower() == 'bobyqa':
        known_field.append('honour_x0')
    unknown_field = list(set(option_fields).difference(set(known_field)))

    # Remove the unknown fields. If we do not removed unknown fields, we may still complain later if an unknown field is
    # not properly set (e.g., options.npt is not a number) even though we have declared that this field will be ignored.
    for key in unknown_field:
        options.pop(key)
    if len(unknown_field) > 0:
        warn_message = '{}: unknown option(s): {}; they are ignored.'.format(invoker, ', '.join(unknown_field))
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)
    option_fields = list(options.keys())  # The fields may have been revised

    # Set default npt according to solver.
    # If method == '' or None, then invoker must be pdfo, and a solver will be selected later; when the solver is
    # chosen, a valid npt will be defined. So we do not need to consider the case with method == '' here. Note we have
    # to take maxfev into consideration when selecting the solver, because npt is at most maxfev-1!
    if method is None:
        # npt is set to NaN to be decided later, depending on the selected solver.
        npt = np.nan
    else:
        if method.lower() in ['bobyqa', 'lincoa', 'newuoa']:
            npt = 2 * lenx0 + 1
        elif method.lower() == 'uobyqa':
            npt = (lenx0 + 1) * (lenx0 + 2) / 2
        else:
            # The method is necessarily COBYLA.
            npt = lenx0 + 1

    # Validate options['scale'] at first. It needs to be done here since the trust-region radii revision required for
    # BOBYQA asks for the true scaling flag.
    validated = False
    if 'scale' in option_fields:
        if not isinstance(options['scale'], (bool, np.bool_)):
            warn_message = '{}: invalid scale flag; it should be True or False; it is set to {}.'.format(invoker, scale)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif options['scale'] and np.any(np.isinf(ub - lb)):
            warn_message = \
                '{}: problem cannot be scaled because not all variables have both lower and upper ' \
                'bounds.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
            options['scale'] = False  # it must be set to False not to be performed
            validated = True
        else:
            validated = True

    if not validated:  # options['scale'] has not got a valid value yet.
        options['scale'] = scale
    options['scale'] = bool(options['scale'])

    # Revise default rhobeg and rhoend if the solver is BOBYQA.
    if options['scale']:
        rhobeg = .5  # this value cannot be bigger than 1. Otherwise, BOBYQA will complain
        rhoend = 1e-6
    if method is not None and method.lower() == 'bobyqa' and not options['scale']:
        lb_mod, ub_mod = lb.copy(), ub.copy()
        lb_mod[np.logical_and(lb_mod > 0, np.isinf(lb_mod))] = -np.inf
        ub_mod[np.logical_and(ub_mod < 0, np.isinf(ub_mod))] = np.inf
        rhobeg_bobyqa = min(rhobeg, np.min(ub_mod - lb_mod) / 4)
        rhobeg = max(eps, rhobeg_bobyqa)
        rhoend = max(eps, min(.1 * rhobeg, rhoend))

    # Validate the user-specified options (except scale that has already been done); adopt the default values if needed.
    # Validate options['npt'].
    # There are the following possibilities.
    # 1. The user specifies options['npt']
    # 1.1. The solver is yet to decide (method=None): we keep options['npt'] if it is a positive integer; otherwise,
    #      raise a warning and set options['npt'] to NaN;
    # 1.2. The user has chosen a valid solver: we keep options['npt'] if it is compatible with the solver; otherwise,
    #      raise a warning and set options['npt'] to the default value according to the solver.
    # 2. The user does not specify options['npt']
    # 2.1. The solver is yet to decide (solver=None): we set options['npt'] to NaN;
    # 2.2. The user has chosen a valid solver: we set options.npt to the default value according to the solver. After
    #      this process, options['npt'] is either a positive integer (compatible with method if it is specified by the
    #      user) or NaN (only if the user does not specify a valid solver while options.npt is either unspecified or not
    #      a positive integer).
    validated = False
    if 'npt' in option_fields:
        integer_types = (int, np.int, np.int8, np.int16, np.int32, np.int64)
        if method is not None and \
                (not isinstance(options['npt'], integer_types) or options['npt'] < 1 or np.isnan(options['npt'])):
            warn_message = '{}: invalid npt. It should be a positive integer.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif method is not None and method.lower() in ['bobyqa', 'lincoa', 'newuoa'] and \
                (not isinstance(options['npt'], integer_types) or np.isnan(options['npt']) or
                 options['npt'] < lenx0 + 2 or options['npt'] > (lenx0 + 1) * (lenx0 + 2) / 2):
            warn_message = \
                '{}: invalid npt. {} requires it to be an integer and n+2 <= npt <= (n+1)*(n+2)/2; it is set to ' \
                '2n+1.'.format(invoker, method)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['npt'] has not got a valid value yet
        # For uobyqa and cobyla or empty solver, we adopt the 'default npt' defined above, although it will NOT be used
        # by the solver.
        options['npt'] = npt

    # If prepdfo is called by pdfo, and if pdfo has been called without precising the method, npt will be set to np.nan.
    # If we do not check whether this value is np.nan before casting it into np.int32, it will lead to an error.
    if not np.isnan(options['npt']):
        options['npt'] = np.int32(options['npt'])

    # Validate options['maxfev'].
    validated = False
    if 'maxfev' in option_fields:
        if not isinstance(options['maxfev'], scalar_types) or options['maxfev'] <= 0 or np.isnan(options['maxfev']) or \
                np.isinf(options['maxfev']):
            warn_message = \
                '{}: invalid maxfev; it should be a positive integer; it is set to {}.'.format(invoker, maxfev)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif method is None and options['maxfev'] <= lenx0 + 1:
            options['maxfev'] = lenx0 + 2  # The smallest possible value
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
                # The method is necessarily COBYLA.
                warn_message = '{}: invalid maxfev; {} requires maxfev > n+1; it is set to n+2.'.format(invoker, method)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['maxfev'] has not got a valid value yet.
        options['maxfev'] = max(maxfev, npt+1)
    if not np.isnan(options['maxfev']):
        options['maxfev'] = np.int32(options['maxfev'])

    # Validate options['rhobeg'].
    # NOTE: if the problem is to be scaled, then options['rhobeg'] and options['rhoend'] will be used as the intial and
    # final trust-region radii for the scaled problem.
    validated = False
    if 'rhobeg' in option_fields:
        if not isinstance(options['rhobeg'], scalar_types) or options['rhobeg'] <= 0 or np.isnan(options['rhobeg']) or \
                np.isinf(options['rhobeg']):
            warn_message = \
                '{}: invalid rhobeg; it should be a positive number; it is set to ' \
                'max({}, rhoend).'.format(invoker, rhobeg)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif method is not None and method.lower() == 'bobyqa':
            if options['scale'] and options['rhobeg'] > 1:
                warn_message = \
                    '{}: invalid rhobeg; {} requires rhobeg <= 1 when the problem is scaled; it is set to ' \
                    '0.5.'.format(invoker, method)
                warnings.warn(warn_message, Warning)
                list_warnings.append(warn_message)
                options['rhobeg'] = 0.5
            elif not options['scale'] and options['rhobeg'] > np.min(ub - lb) / 2:
                warn_message = \
                    '{}: invalid rhobeg; {} requires rhobeg <= min(ub-lb)/2; it is set to ' \
                    'min(ub-lb)/4.'.format(invoker, method)
                warnings.warn(warn_message, Warning)
                list_warnings.append(warn_message)
                options['rhobeg'] = np.min(ub - lb) / 4
            validated = True  # The value in options['rhobeg'] needs to be used.
        else:
            validated = True

    if not validated:  # options['rhobeg'] has not got a valid value yet.
        if 'rhoend' in option_fields and isinstance(options['rhoend'], scalar_types):
            options['rhobeg'] = max(rhobeg, 1e1 * options['rhoend'])
        else:
            options['rhobeg'] = rhobeg
    options['rhobeg'] = np.float64(max(options['rhobeg'], eps))

    # Validate options['rhoend'].
    validated = False
    if 'rhoend' in option_fields:
        if not isinstance(options['rhoend'], scalar_types) or options['rhoend'] > options['rhobeg'] or \
                0 >= options['rhoend'] or np.isnan(options['rhobeg']):
            warn_message = \
                '{}: invalid rhoend; we should have rhobeg >= rhoend > 0; it is set to ' \
                'min(0.1*rhobeg, {}).'.format(invoker, rhoend)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['rhoend'] has not got a valid value yet.
        options['rhoend'] = min(0.1 * options['rhobeg'], rhoend)
    options['rhoend'] = np.float64(np.nanmax((options['rhoend'], eps)))

    # Validate options['ftarget'].
    validated = False
    if 'ftarget' in option_fields:
        if not isinstance(options['ftarget'], scalar_types) or np.isnan(options['ftarget']):
            warn_message = '{}: invalid ftarget; it should be real number; it is set to {}.'.format(invoker, ftarget)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['ftarget'] has not got a valid value yet.
        options['ftarget'] = ftarget
    options['ftarget'] = np.float64(options['ftarget'])

    # Validate options['classical'].
    validated = False
    if 'classical' in option_fields:
        if not isinstance(options['classical'], (bool, np.bool_)):
            warn_message = \
                '{}: invalid scale flag; it should be True or False; it is set to {}.'.format(invoker, classical)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['classical'] has not got a valid value yet.
        options['classical'] = classical
    options['classical'] = bool(options['classical'])
    if options['classical']:
        warn_message = \
            "{}: in classical mode, which is recommended only for research purpose; set options['classical']=False " \
            "to disable classical mode.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    # Validate options['eliminate_lin_eq'].
    validated = False
    if 'eliminate_lin_eq' in option_fields:
        if not isinstance(options['eliminate_lin_eq'], (bool, np.bool_)):
            warn_message = \
                '{}: invalid eliminate_lin_eq flag; it should be True or False; it is set to ' \
                '{}.'.format(invoker, eliminate_lin_eq)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['eliminate_lin_eq'] has not got a valid value yet.
        options['eliminate_lin_eq'] = eliminate_lin_eq
    options['eliminate_lin_eq'] = bool(options['eliminate_lin_eq'])

    # Validate options['honour_x0'].
    validated = False
    if 'honour_x0' in option_fields:
        if not isinstance(options['honour_x0'], (bool, np.bool_)):
            warn_message = \
                '{}: invalid honour_x0 flag; it should be True or False; it is set to {}.'.format(invoker, honour_x0)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['honour_x0'] has not got a valid value yet.
        options['honour_x0'] = honour_x0
    options['honour_x0'] = bool(options['honour_x0'])

    # Validate options['quiet'].
    validated = False
    if 'quiet' in option_fields:
        if not isinstance(options['quiet'], (bool, np.bool_)):
            warn_message = '{}: invalid quiet flag; it should be True or False; it is set to {}.'.format(invoker, quiet)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['quiet'] has not got a valid value yet.
        options['quiet'] = quiet
    options['quiet'] = bool(options['quiet'])

    # Validate options['debug'].
    validated = False
    if 'debug' in option_fields:
        if not isinstance(options['debug'], (bool, np.bool_)):
            warn_message = \
                '{}: invalid debug flag; it should be True or False; it is set to {}.'.format(invoker, debugflag)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['debug'] has not got a valid value yet.
        options['debug'] = debugflag
    options['debug'] = bool(options['debug'])
    if options['debug']:
        warn_message = "{}: in debug mode; set options['debug']=False to disable debug.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)
        if options['quiet']:
            options['quiet'] = False
            warn_message = "options['quiet'] is set to False because options['debug'] = True.".format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    # Validate options['chkfunval'].
    validated = False
    if 'chkfunval' in option_fields:
        if not isinstance(options['chkfunval'], (bool, np.bool_)):
            warn_message = \
                '{}: invalid chkfunval flag; it should be True or False; it is set to {}.'.format(invoker, chkfunval)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        elif options['chkfunval'] and not options['debug']:
            warn_message = \
                '{}: chkfunval = True but debug = False; chkfunval is set to false; set both flags to true to check ' \
                'function values.'.format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)
        else:
            validated = True

    if not validated:  # options['chkfunval'] has not got a valid value yet.
        options['chkfunval'] = chkfunval
    options['chkfunval'] = bool(options['chkfunval'])
    if options['chkfunval']:
        warn_message = \
            "{}: checking whether fx = fun(x) and possibly conval = con(x) at exit, which will cost an extra " \
            "function/constraint evaluation; set options['chkfunval'] = False to disable the check.".format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    return options, user_option_fields, method


def _constr_violation(invoker, x, lb, ub, constraints, prob_info):
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
    prob_info: dict
        The same as in prepdfo.

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

    # Compute the constraint violation as the absolute distance to the bounds.
    constr_violation = np.max((0, np.nanmax(np.r_[lb - x, x - ub])))

    if constraints['linear'] is not None:
        a, b = _linear_constraints_constr(constraints['linear'])

        # Compute the constraint violation as the largest absolute distance to the bounds and the linear constraints.
        if b.size > 0:
            constr_violation = np.max((constr_violation, np.nanmax(np.dot(a, x) - b)))

    if constraints['nonlinear'] is not None:
        nonlinear = constraints['nonlinear']

        if not isinstance(nonlinear, dict) or not ({'type', 'fun'} <= set(nonlinear.keys())) or \
                nonlinear['type'] != 'ineq' or not callable(nonlinear['fun']) or \
                (python_version >= 3 and len(signature(nonlinear['fun']).parameters) == 0):
            raise ValueError('{}: UNEXPECTED ERROR: the nonlinear constraints are ill-defined.'.format(invoker))

        nlc = np.asarray(nonlinear['fun'](x), dtype=np.float64)
        if prob_info['infeasible']:
            # The variable nlc should be augmented to calculate the true constraint violation, with respect to the
            # bounds. In this case, the variable nlc_aug contains the evaluations of the constraint functions, in order.
            k_nonlinear = 0
            nlc_aug = np.asarray([], dtype=np.float64)
            try:
                nonlinear_indices = prob_info['constr_meta']['nonlinear_indices']
                for metadata in [prob_info['constr_meta']['data'][k] for k in nonlinear_indices]:
                    # When this section of the code is executed, the nonlinearly-constrained problem turned infeasible
                    # and the constraint functions have been evaluated. Therefore, all the metadata have been set,
                    # possibly with the value -inf or +inf with a lower or a upper was not defined. Thus, since we
                    # compute the constraint violation as the max of all constraint violations, we can augment the
                    # constraints on every variable.
                    nlc_current = nlc[k_nonlinear:k_nonlinear + metadata['len']]
                    nlc_aug = np.r_[nlc_aug, metadata['lb'] - nlc_current, nlc_current - metadata['ub']]
                    k_nonlinear += metadata['len']
            except (KeyError, IndexError):
                raise ValueError('{}: UNEXPECTED ERROR: the constraints metadata are ill-defined.'.format(invoker))
        else:
            nlc_aug = nlc

        # Compute the constraint violation as the largest relative distance to the bounds, the linear constraints and
        # the nonlinear constraints.
        constr_violation = max(constr_violation, np.max(nlc_aug))
    else:
        nlc = np.asarray([], dtype=np.float64)

    return np.float64(constr_violation), nlc


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

    # Define the type of the given problem.
    if constraints['nonlinear'] is not None:
        return 'nonlinearly-constrained'
    elif constraints['linear'] is not None and (constraints['linear'].lb.size > 0 or constraints['linear'].ub.size > 0):
        return 'linearly-constrained'
    elif (len(lb) > 0 and np.max(lb) > -np.inf) or (len(ub) > 0 and np.min(ub) < np.inf):
        return 'bound-constrained'
    else:
        return 'unconstrained'


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

    # Convert {lb_i <= A_i*x <= ub_i}_i into one constraint A*x <= b.
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
    fun_name = stack()[0][3]  # name of the current function

    if not (hasattr(freex_value, '__len__') and hasattr(fixedx_value, '__len__') and
            hasattr(freex, '__len__') and hasattr(fixedx, '__len__')):
        raise ValueError('{}: UNEXPECTED ERROR: the variable arrays have wrong types.'.format(fun_name))

    freex_value = np.asarray(freex_value, dtype=np.float64)
    fixedx_value = np.asarray(fixedx_value, dtype=np.float64)
    freex = np.asarray(freex, dtype=bool)
    fixedx = np.asarray(fixedx, dtype=bool)
    if freex.size != fixedx.size or freex_value.size + fixedx_value.size != freex.size:
        raise ValueError('{}: UNEXPECTED ERROR: the variable vector lengths are inconsistent'.format(fun_name))

    # Build the complete vector from the fixed and free values.
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
        # Essentially do nothing. DO NOT remove this case. Otherwise, the case would be included in 'else', which is not
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
            '{}: UNEXPECTED ERROR: the sizes of the initial guess and the bounds are inconsistent.'.format(invoker))

    if not isinstance(list_warnings, list):
        raise ValueError('{}: UNEXPECTED ERROR: the list of warnings is ill-defined.'.format(invoker))

    # x_before_scaling = scaling_factor*x_after_scaling + shift.

    # Question: What about scaling according to the magnitude of x0, lb, ub, x0-lb, ub-x0?
    # This can be useful if lb and ub reflect the nature of the problem well, and x0 is a reasonable approximation to
    # the optimal solution. Otherwise, it may be a bad idea.
    substantially_scaled_threshold = 2

    # We consider the problem substantially scaled_threshold if
    # max([1; scaling_factor])/min([1; scaling_factor]) > substantially_scaled_threshold

    # The strategy for scaling has changed with v1.0: do not scale the problem unless all variables have both lower and
    # upper bounds
    # lenx0 = x0.size
    # index_lub = np.logical_and(lb > -np.inf, ub < np.inf)
    # scaling_factor = np.ones(lenx0, dtype=np.float64)
    # shift = np.zeros(lenx0, dtype=np.float64)
    if np.any(np.isinf(ub - lb)):
        raise ValueError(
            '{}: UNEXPECTED ERROR: at least one of [-lb; ub] is infinity. Scaling should not be '
            'performed.'.format(invoker))

    # Define the scaling factor and the shift according to the bounds.
    scaling_factor = (ub - lb) / 2
    shift = (ub + lb) / 2

    # Build the scaled objective function.
    def fun_c(x):
        return np.float64(fun(scaling_factor * x + shift))

    # Scale the initial guess and the bounds.
    x0_c = (x0 - shift) / scaling_factor
    lb_c = (lb - shift) / scaling_factor
    ub_c = (ub - shift) / scaling_factor

    constraints_c = {'linear': None, 'nonlinear': None}

    # Scale the linear constraints.
    if constraints['linear'] is not None:
        a = constraints['linear'].A
        lb = constraints['linear'].lb
        ub = constraints['linear'].ub
        ashift = np.dot(a, shift)
        constraints_c['linear'] = LinearConstraint(np.dot(a, np.diag(scaling_factor)), lb=lb - ashift, ub=ub - ashift)

    # Scale the nonlinear constraints.
    if constraints['nonlinear'] is not None:
        constraints_c['nonlinear'] = \
            {'type': 'ineq', 'fun': lambda x: constraints['nonlinear']['fun'](scaling_factor * x + shift)}

    # From v1.0, we do not warn about scaling anymore. Scaling works well in several real problems.
    # if any(scaling_factor != 1):
    #     warn_message = \
    #         "{}: problem scaled according to bound constraints; do this only if the bounds reflect the scaling of " \
    #         "variables; if not, set options['scale']=False to disable scaling.".format(invoker)
    #     warnings.warn(warn_message, Warning)
    #     list_warnings.append(warn_message)

    substantially_scaled = False

    # Check whether the scaling is substantial. If this is true, the options rhobeg and rhoend may be updated for
    # BOBYQA.
    # if (max([scaling_factor; 1./scaling_factor]) > substantially_scaled_threshold)
    if np.max(np.r_[1, scaling_factor]) / np.min(np.r_[1, scaling_factor]) > substantially_scaled_threshold:
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

    # Validate invoker.
    if not isinstance(invoker, str):
        raise ValueError('unknown: UNEXPECTED ERROR: invoker should be a string.')

    # Validate method.
    if method is not None and not isinstance(method, str):
        raise ValueError('{}: UNEXPECTED ERROR: method should be a string.'.format(invoker))

    # Validate options.
    option_fields = {'maxfev', 'rhobeg', 'rhoend'}
    if options is None or not isinstance(options, dict) or not (option_fields <= set(options.keys())) or \
            not isinstance(options['maxfev'], scalar_types) or not isinstance(options['rhobeg'], scalar_types) or \
            not isinstance(options['rhoend'], scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: options should be a valid dictionary.'.format(invoker))

    # Validate prob_info.
    prob_info_fields = {'refined_type', 'refined_dim'}
    if prob_info is None or not isinstance(prob_info, dict) or not (prob_info_fields <= set(prob_info.keys())) or \
            not isinstance(prob_info['refined_type'], str) or not isinstance(prob_info['refined_dim'], scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: prob_info should be a valid dictionary.'.format(invoker))

    # Validate list_warnings.
    if not hasattr(list_warnings, '__len__'):
        raise ValueError('{}: UNEXPECTED ERROR: list_warnings should be a list.'.format(invoker))

    solver = method
    ptype = prob_info['refined_type']
    n = prob_info['refined_dim']

    if solver is None or not _prob_solv_match(ptype, solver):
        if solver is not None:
            # Do not complain if solver is None, i.e., if it has not been provided.
            warn_message = \
                '{}: {} cannot solve a {} problem; {} will select a solver ' \
                'automatically.'.format(invoker, solver, ptype.replace('-', ' '), invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

        # Define the solver depending on the problem characteristics.
        if ptype == 'unconstrained':
            if 2 <= n <= 8 and options['maxfev'] >= (n + 1) * (n + 2) / 2:
                solver = 'uobyqa'  # does not need options['npt']
            elif options['maxfev'] <= n + 2:  # options['maxfev'] == n + 2
                solver = 'cobyla'  # does not need options['npt']
            else:
                # Interestingly, we note in our tests that LINCOA may outperformed NEWUOA on unconstrained CUTEst
                # problems when the dimension was not large (i.e., <= 50) or the precision requirement was not high
                # (i.e., >= 1e-5). Therefore, it is worthwhile to try LINCOA when an unconstrained problem is given.
                # Nevertheless, for the moment, we set the default solver for unconstrained problems to be newuoa, since
                # the solver was intended to solve these problems.
                solver = 'newuoa'
        elif ptype == 'bound-constrained':
            if options['maxfev'] <= n + 2:
                solver = 'cobyla'  # Does not need options['npt']
            else:
                solver = 'bobyqa'
        elif ptype == 'linearly-constrained':
            if options['maxfev'] <= n + 2:
                solver = 'cobyla'  # Does not need options['npt']
            else:
                solver = 'lincoa'
        elif ptype == 'nonlinearly-constrained':
            solver = 'cobyla'  # does not need options['npt']
        else:
            raise SystemError("{}: UNEXPECTED ERROR: unknown problem type '{}' received.".format(fun_name, ptype))

    # Revise options['npt'] according to the selected solver.
    if solver.lower() in ['bobyqa', 'lincoa', 'newuoa'] and \
            (np.isnan(options['npt']) or options['npt'] < n + 2 or
             options['npt'] > min((n + 1)*(n + 2) // 2, options['maxfev'] - 1)):
        options['npt'] = min(2 * n + 1, options['maxfev'] - 1)
        if 'npt' in prob_info['user_options_fields']:
            warn_message = \
                '{}: npt is set to {} according to the selected solver {}, which requires n+2 <= npt <= ' \
                '(n+1)*(n+2)/2.'.format(invoker, options['npt'], solver)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    # Revise options['rhobeg'] and options['rhoend'] according to the selected solver.
    # For the moment, only BOBYQA needs such a revision.
    if solver.lower() == 'bobyqa' and \
            options['rhobeg'] > np.min(prob_info['refined_data']['ub'] - prob_info['refined_data']['lb']) / 2:
        lb_mod, ub_mod = prob_info['refined_data']['lb'].copy(), prob_info['refined_data']['ub'].copy()
        lb_mod[np.logical_and(lb_mod > 0, np.isinf(lb_mod))] = -np.inf
        ub_mod[np.logical_and(ub_mod < 0, np.isinf(ub_mod))] = np.inf
        options['rhobeg'] = max(eps, np.min(ub_mod - lb_mod) / 4)
        options['rhoend'] = max(eps, min(0.1 * options['rhobeg'], options['rhoend']))
        if 'rhobeg' in prob_info['user_options_fields'] or 'rhoend' in prob_info['user_options_fields']:
            warn_message = \
                '{}: rhobeg is set to {} and rhoend to {} acccording to the selected solver bobyqa, which requires ' \
                'rhoend <= rhobeg <= min(ub-lb)/2.'.format(invoker, options['rhobeg'], options['rhoend'])
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    if solver not in solver_list or not _prob_solv_match(ptype, solver):
        raise SystemError("{}: UNEXPECTED ERROR: invalid solver '{}' selected.".format(fun_name, solver))

    return solver


def _pre_rhobeg_x0(invoker, x0, lb, ub, user_options_fields, options, list_warnings):
    # possible solvers
    fun_name = stack()[0][3]  # name of the current function
    local_invoker_list = ['pdfo', 'bobyqa']

    if invoker not in local_invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(local_invoker_list)))
    invoker = stack()[1][3].lower()

    # Validate invoker.
    if not isinstance(invoker, str):
        raise ValueError('unknown: UNEXPECTED ERROR: invoker should be a string.')

    if not (hasattr(x0, '__len__') and hasattr(lb, '__len__') and hasattr(ub, '__len__')):
        raise TypeError('{}: UNEXPECTED ERROR: the initial guess and the bounds should be vectors.'.format(invoker))

    if not hasattr(user_options_fields, '__len__'):
        raise TypeError('{}: UNEXPECTED ERROR: the user defined option fields should be a list.'.format(invoker))

    # Validate options.
    option_fields = {'honour_x0', 'rhobeg', 'rhoend'}
    if options is None or not isinstance(options, dict) or not (option_fields <= set(options.keys())) or \
            not isinstance(options['honour_x0'], (bool, np.bool_)) or \
            not isinstance(options['rhobeg'], scalar_types) or not isinstance(options['rhoend'], scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: options should be a valid dictionary.'.format(invoker))

    # Validate list_warnings.
    if not hasattr(list_warnings, '__len__'):
        raise ValueError('{}: UNEXPECTED ERROR: list_warnings should be a list.'.format(invoker))

    if options['honour_x0']:
        # In this case, we respect the user-defined x0 and revise rhobeg.
        rhobeg_old = options['rhobeg']
        lbx = np.logical_and(np.logical_or(np.logical_not(np.isinf(lb)), lb > 0),
                             x0 - lb <= eps * np.max(np.r_[1, np.abs(lb)]))
        nlbx = np.logical_not(lbx)
        ubx = np.logical_and(np.logical_or(np.logical_not(np.isinf(ub)), ub < 0),
                             x0 - ub >= -eps * np.max(np.r_[1, np.abs(ub)]))
        nubx = np.logical_not(ubx)
        options['rhobeg'] = max(eps, np.min(np.r_[options['rhobeg'], x0[nlbx] - lb[nlbx], ub[nubx] - x0[nubx]]))
        options['rhoend'] = min(options['rhoend'], options['rhobeg'])
        x0[lbx] = lb[lbx]
        x0[ubx] = ub[ubx]
        if rhobeg_old - options['rhobeg'] > eps * max(1, rhobeg_old):
            # If the user does not specify rhobeg, no warning should be raised.
            options['rhoend'] = max(eps, min(options['rhoend'], .1 * options['rhobeg']))
            if 'rhobeg' in user_options_fields or 'rhoend' in user_options_fields:
                warn_message = \
                    '{}: rhobeg is revised so that the distance between x0 and the inactive bounds is at least ' \
                    'rhobeg.'.format(invoker)
                warnings.warn(warn_message, Warning)
                list_warnings.append(warn_message)
    else:
        x0_old = x0.copy()
        lbx = x0 <= lb + options['rhobeg'] / 2
        lbx_plus = np.logical_and(x0 > lb + options['rhobeg'] / 2, x0 < lb + options['rhobeg'])
        ubx = x0 >= ub - options['rhobeg'] / 2
        ubx_minus = np.logical_and(x0 < ub - options['rhobeg'] / 2, x0 > ub - options['rhobeg'])
        x0[lbx] = lb[lbx]
        x0[lbx_plus] = lb[lbx_plus] +options['rhobeg']
        x0[ubx_minus] = ub[ubx_minus] - options['rhobeg']
        x0[ubx] = ub[ubx]
        if np.linalg.norm(x0_old - x0) > eps * max(1, np.linalg.norm(x0_old)):
            warn_message = \
                "{}: x0 is revised so that the distance between x0 and the inactive bounds is at least rhobeg; set " \
                "options['honour_x0']=True if you prefer to keep x0.".format(invoker)
            warnings.warn(warn_message, Warning)
            list_warnings.append(warn_message)

    return x0, options


def _pre_npt_elimination(invoker, n, user_options_fields, options, list_warnings):
    npt_old = options['npt']
    options['npt'] = max(min(npt_old, (n + 1) * (n + 2) // 2), n + 2)

    if npt_old != options['npt'] and 'npt' in user_options_fields:
        warn_message = \
            '{}: the dimension of the searching space has been reduced, and npt must be revised to satisfy the ' \
            'requirements.'.format(invoker)
        warnings.warn(warn_message, Warning)
        list_warnings.append(warn_message)

    return options


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

    # Validate x0.
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

    # Validate lb.
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

    # Validate ub.
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

    # Validate constraints.
    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())) or \
            not (isinstance(constraints['linear'], LinearConstraint) or constraints['linear'] is None):
        # the nonlinear constraints will not be taken into account in this function and are, therefore, not validated
        raise ValueError('{}: UNEXPECTED ERROR: The constraints are ill-defined.'.format(invoker))

    # Validate options
    if options is not None and not isinstance(options, dict):
        raise ValueError('{}: UNEXPECTED ERROR: The options should be a dictionary.'.format(invoker))

    max_con = 1e20  # Decide whether an inequality constraint can be ignored

    # Projecte onto the feasible set.
    if constraints['linear'] is None:
        # Direct projection onto the bound constraints
        x_proj = np.nanmin((np.nanmax((x0_c, lb_c), axis=0), ub_c), axis=0)
        return OptimizeResult(x=x_proj)
    elif all(np.less_equal(np.abs(constraints['linear'].ub - constraints['linear'].lb), eps)) and \
            np.max(lb_c) <= -max_con and np.min(ub_c) >= max_con:
        # The linear constraints are all equality constraints. The projection can therefore be done by solving the
        # least-squares problem: min ||A*x - (b - A*x_0)||.
        a = constraints['linear'].A
        b = (constraints['linear'].lb + constraints['linear'].ub) / 2
        xi, _, _, _ = np.linalg.lstsq(a, b - np.dot(a, x0_c), rcond=None)

        # The problem is not bounded. However, if the least-square solver returned values bigger in absolute value
        # than max_con, they will be reduced to this bound.
        x_proj = np.nanmin((np.nanmax((x0_c + xi, lb_c), axis=0), ub_c), axis=0)

        return OptimizeResult(x=x_proj)

    if constraints['linear'] is not None:
        try:
            # Project the initial guess onto the linear constraints via SciPy.
            from scipy.optimize import minimize
            from scipy.optimize import Bounds as ScipyBounds
            from scipy.optimize import LinearConstraint as ScipyLinearConstraint

            linear = constraints['linear']

            # To be more efficient, SciPy asks to separate the equality and the inequality constraints into two
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

            if pc_args_ineq['A'].size > 0 and pc_args_ineq['lb'].size > 0 and pc_args_eq['lb'].size > 0:
                project_constraints = [ScipyLinearConstraint(**pc_args_ineq), ScipyLinearConstraint(**pc_args_eq)]
            elif pc_args_ineq['A'].size > 0 and pc_args_ineq['lb'].size > 0:
                project_constraints = ScipyLinearConstraint(**pc_args_ineq)
            elif pc_args_ineq['A'].size > 0:
                project_constraints = ScipyLinearConstraint(**pc_args_eq)
            else:
                project_constraints = ()

            # Perform the actual projection.
            ax_ineq = np.dot(pc_args_ineq['A'], x0_c)
            ax_eq = np.dot(pc_args_eq['A'], x0_c)
            if np.greater(ax_ineq, pc_args_ineq['ub']).any() or np.greater(pc_args_ineq['lb'], ax_ineq).any() or \
                    np.not_equal(ax_eq, pc_args_eq['lb']).any() or \
                    np.greater(x0_c, ub_c).any() or np.greater(lb_c, x0_c).any():
                return minimize(lambda x: np.dot(x - x0_c, x - x0_c) / 2, x0_c, jac=lambda x: (x - x0_c),
                                bounds=ScipyBounds(lb_c, ub_c), constraints=project_constraints)
            else:
                # Do not perform any projection if the initial guess is feasible.
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

    # Construct the linear constraints that refers to the bounds.
    idmatrix = np.eye(n)
    lb, ub = bounds['lb'], bounds['ub']
    lb_kept_indices = np.logical_not(np.logical_and(np.isinf(lb), lb < 0))
    ub_kept_indices = np.logical_not(np.logical_and(np.isinf(ub), ub > 0))
    alb = idmatrix[lb_kept_indices, :]
    aub = idmatrix[ub_kept_indices, :]

    # Reshape the empty matrices to avoid concatenate exception.
    if aub.size == 0:
        aub = aub.reshape(0, n)
    if alb.size == 0:
        alb = alb.reshape(0, n)

    # Remove infinite bounds.
    lb, ub = lb[lb_kept_indices], ub[ub_kept_indices]

    # Construct of the actual augmented matrices.
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

    fun_name = stack()[0][3]  # name of the current function

    if len(stack()) < 3 or stack()[1][3].lower() not in invoker_list:
        raise SystemError('`{}` should only be called by {}'.format(fun_name, ', '.join(invoker_list)))
    invoker = stack()[1][3].lower()

    # Validate x.
    if not hasattr(x, '__len__') and \
            not isinstance(x, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: x should be a scalar or a vector.'.format(invoker))
    try:
        x_c = np.asarray(x, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: x should contain only scalars.'.format(invoker))
    if len(x_c.shape) > 1:
        raise ValueError('{}: UNEXPECTED ERROR: x should be a vector.'.format(invoker))

    # Validate fx.
    if hasattr(fx, '__len__') and len(fx) == 1:
        fx_c = np.float64(fx[0])
    else:
        fx_c = np.float64(fx)
    if not isinstance(fx_c, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: fx should be a scalar.'.format(invoker))

    # Validate exitflag.
    if not isinstance(exitflag, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: exitflag should be a scalar.'.format(invoker))
    exitflag_c = np.int32(exitflag)
    if exitflag_c != exitflag:
        raise ValueError('{}: UNEXPECTED ERROR: exitflag should not be a floating number.'.format(invoker))

    # Validate output.
    if output is None or not isinstance(output, dict):
        raise ValueError('{}: UNEXPECTED ERROR: output should be a valid dictionary.'.format(invoker))

    # Validate method.
    if method is None or not isinstance(method, str):
        raise ValueError('{}: UNEXPECTED ERROR: method should be a string.'.format(invoker))

    # Validate nf.
    if not isinstance(nf, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: nf should be a scalar.'.format(invoker))
    nf_c = np.int32(nf)
    if nf_c != nf:
        raise ValueError('{}: UNEXPECTED ERROR: nf should not be a floating number.'.format(invoker))

    # Validate fhist.
    if not hasattr(fhist, '__len__') and not isinstance(fhist, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: fhist should be a scalar of a vector.'.format(invoker))
    try:
        fhist_c = np.asarray(fhist[:nf], dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: fhist should contain nf scalars.'.format(invoker))
    if len(fhist_c.shape) != 1:
        raise ValueError('{}: UNEXPECTED ERROR: fhist should be a vector.'.format(invoker))

    # Validate constrviolation.
    if not np.isnan(constrviolation) and not isinstance(constrviolation, scalar_types):
        raise ValueError('{}: UNEXPECTED ERROR: constrviolation should be a scalar.'.format(invoker))
    if np.isnan(constrviolation):
        constrviolation_c = constrviolation
    else:
        constrviolation_c = np.float64(constrviolation)

    # Validate chist.
    if not (chist is None and method in ['pdfo', 'newuoa', 'uobyqa']) and \
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

    # If the invoker is a solver called by pdfo, then let pdfo do the post-processing.
    output['x'] = x_c
    output['fun'] = fx_c
    output['status'] = exitflag_c
    output['success'] = (exitflag_c in [0, 1, 14]) or (exitflag_c == 13 and abs(constrviolation_c) <= 1e-15)
    if len(stack()) >= 4 and stack()[2][3].lower() == 'pdfo':
        output['nfev'] = nf_c
        output['constrviolation'] = constrviolation_c
        output['fhist'] = fhist_c
        output['chist'] = chist_c

        return OptimizeResult(**output)

    # If the solver is not called by pdfo (can be pdfo directly), perform the post-processing.
    option_fields = {'quiet', 'debug', 'classical', 'chkfunval'}
    if options is None or not isinstance(options, dict) or not (option_fields <= set(options.keys())) or \
            not isinstance(options['quiet'], (bool, np.bool_)) or not isinstance(options['debug'], (bool, np.bool_)) or \
            not isinstance(options['classical'], (bool, np.bool_)) or \
            not isinstance(options['chkfunval'], (bool, np.bool_)):
        raise ValueError('{}: UNEXPECTED ERROR: options should be a valid dictionary.'.format(invoker))

    # Validate prob_info.
    prob_info_fields = \
        {'infeasible', 'nofreex', 'warnings', 'scaled', 'reduced', 'space_chg', 'fixedx', 'fixedx_value',
         'refined_type', 'raw_type', 'infeasible_linear', 'infeasible_bound', 'feasibility_problem'}
    if prob_info is None or not isinstance(prob_info, dict) or not (prob_info_fields <= set(prob_info.keys())) or \
            not isinstance(prob_info['infeasible'], (bool, np.bool_)) or \
            not isinstance(prob_info['nofreex'], (bool, np.bool_)) or \
            not hasattr(prob_info['warnings'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, str), prob_info['warnings'])) or \
            not isinstance(prob_info['scaled'], (bool, np.bool_)) or \
            not isinstance(prob_info['reduced'], (bool, np.bool_)) or \
            not (prob_info['space_chg'] is None or callable(prob_info['space_chg'])) or \
            not hasattr(prob_info['fixedx'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, (bool, np.bool_)), prob_info['fixedx'])) or \
            not hasattr(prob_info['fixedx_value'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, scalar_types), prob_info['fixedx_value'])) or \
            not hasattr(prob_info['infeasible_linear'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, (bool, np.bool_)), prob_info['infeasible_linear'])) or \
            not hasattr(prob_info['infeasible_bound'], '__len__') or \
            not all(map(lambda pi: isinstance(pi, (bool, np.bool_)), prob_info['infeasible_bound'])) or \
            not isinstance(prob_info['refined_type'], str) or not isinstance(prob_info['raw_type'], str) or \
            not isinstance(prob_info['feasibility_problem'], (bool, np.bool_)):
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
    prob_info_c = dict(prob_info)

    # Manage the extreme barriers.
    if not options['classical']:
        if ((fhist_c > hugefun).any() or np.isnan(fhist_c).any()) and not prob_info_c['infeasible'] and \
                not prob_info_c['nofreex']:
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
                ((chist_c > hugecon).any() or np.isnan(chist_c).any()) and not prob_info_c['infeasible'] and \
                not prob_info_c['nofreex']:
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} returns an chist with NaN or values larger than hugecon={}; this is '
                'impossible with extreme barrier.'.format(invoker, method, hugecon))
        elif chist_c is not None and chist_c.size > 0 and np.max(chist_c) == hugecon:
            warn_message = '{}: extreme barrier is invoked; function values that are NaN or larger than hugecon={} ' \
                           'are replaced by hugecon.'.format(invoker, hugecon)
            warnings.warn(warn_message, Warning)
            output['warnings'].append(warn_message)

    # Validate the value of the inputs.
    if nf_c <= 0:
        raise ValueError(
            '{}: UNEXPECTED ERROR: {} returns a nf <= 0 unexpectedly with exitflag '
            '{}'.format(invoker, method, exitflag_c))

    # The scaling affects constrviolation when there are bound constraint. Hence constrviolation has to be recalculated
    # so that it equals the constraint violation of the returned x with respect to the original problem.  Ideally, chist
    # should also be recalculated. However, it is impossible because we do not save the history of x. Therefore, when
    # prob_info['scaled'] == True, chist is not the history of constraint violation of the original problem but the
    # scaled one. It it not consistent with constrviolation. Without saving of history of x, we cannot do better.
    # Before recalculating constrviolation, save the one returned by the solver, because it will be used in debug mode
    # when checking whether fx is consistent with fhist and chist. See the definition of fhist for details.
    constrv_returned = constrviolation_c
    if prob_info_c['scaled']:
        # First calculate the residuals of the linear constraints. This must be calculated before x is scaled back.
        # Otherwise, we would have to scale also the linear constraints back to get the correct residuals. When
        # descaling the bound information contained in the linear constraints, we should not apply any shift, since this
        # information is already contained in the evaluation.
        linear = prob_info_c['refined_data']['constraints']['linear']
        if linear is not None:
            ax = np.dot(linear.A, x_c)
            r1 = linear.lb - ax
            r2 = ax - linear.ub
            if prob_info_c['space_chg'] is not None:
                to_scale = prob_info_c['bounds_in_lin_eq']
                if to_scale.size > 0:
                    r1[-to_scale.size:] *= prob_info_c['scaling_factor'][to_scale]
                    r2[-to_scale.size:] *= prob_info_c['scaling_factor'][to_scale]
            r = np.r_[r1, r2]
        else:
            r = np.asarray([np.nan])

        if prob_info_c['space_chg'] is not None:
            # The space was (possibly) changed, change it back.
            x_c = prob_info_c['space_chg'](x_c)
            lb = prob_info_c['space_chg'](prob_info_c['refined_data']['lb'])
            ub = prob_info_c['space_chg'](prob_info_c['refined_data']['ub'])
        else:
            lb = prob_info_c['refined_data']['lb']
            ub = prob_info_c['refined_data']['ub']

        # Scale x back.
        x_c = prob_info_c['scaling_factor'] * x_c + prob_info_c['shift']

        # Scale the bounds back.
        lb *= prob_info_c['scaling_factor']
        ub *= prob_info_c['scaling_factor']
        lb += prob_info_c['shift']
        ub += prob_info_c['shift']

        # We only need to calculate constrviolation for lincoa and cobyla, because uobyqa and newuoa do not handle
        # constrained problems, while bobyqa is a feasible method and should return constrviolation = 0 regardless of
        # the scaling unless something goes wrong.
        if prob_info_c['space_chg'] is not None:
            not_to_scale = prob_info_c['bounds_in_lin_eq']
            to_scale = np.array([i for i in np.arange(x_c.size) if i not in not_to_scale])
            if to_scale.size > 0:
                conv_bounds = np.r_[lb[to_scale] - x_c[to_scale], x_c[to_scale] - ub[to_scale]]
            else:
                conv_bounds = [-np.inf]
        elif not prob_info_c['nofreex']:
            conv_bounds = np.r_[lb - x_c, x_c - ub]
        else:
            conv_bounds = [-np.inf]
        if method == 'lincoa':
            conv_n = np.concatenate((r, conv_bounds))
            conv_n = np.nanmax((np.zeros_like(conv_n), conv_n), axis=0)
            constrviolation_c = np.max(conv_n)
        else:
            # Compute the constraint violation as usual.
            nlc = np.asarray([-np.inf], dtype=np.float64)
            if 'constr_value' in output.keys():
                nlc = np.asarray(output['constr_value'], dtype=np.float64)
            conv = np.concatenate((r, conv_bounds, nlc))
            if np.isnan(conv).all():
                constrviolation_c = np.nan
            else:
                constrviolation_c = np.nanmax(np.append(conv, 0))
    elif prob_info_c['space_chg'] is not None:
        # The space was (possibly) changed, change it back.
        x_c = prob_info_c['space_chg'](x_c)

    # The problem was (possibly) reduced, get the full x.
    if prob_info_c['reduced'] and not prob_info['nofreex']:
        x_c = _fullx(x_c, prob_info_c['fixedx_value'], np.logical_not(prob_info_c['fixedx']), prob_info_c['fixedx'])
    output['x'] = x_c

    # Set output.{nf, constrviolation, fhist, chist, method}.
    output['nfev'] = nf_c
    output['constrviolation'] = constrviolation_c
    output['fhist'] = fhist_c
    output['chist'] = chist_c
    output['method'] = method

    # If the problem is a feasibility problem, set fx to an empty array and remove fhist from the output
    if prob_info['feasibility_problem']:
        output['fun'] = None
        del output['fhist']

        if prob_info['refined_type'] != 'nonlinearly-constrained':
            # No function evaluation involved when solving a linear feasibility problem. By "function evaluation", we
            # mean the evaluation of the objective function and nonlinear constraint functions, which do not exist in
            # this case. For nonlinear feasibility problems, funcCount is positive.
            output['nfev'] = 0

    # Revise constrviolation and chist according to problem type.
    # max_c = 0 if chist_c is None or chist_c.size == 0 else np.nanmax(chist_c)
    # if prob_info_c['refined_type'] == 'unconstrained' and (constrviolation_c > 0 or max_c > 0):
    #     raise ValueError(
    #         '{}: UNEXPECTED ERROR: {} returns positive constrviolations for an unconstrained '
    #         'problem.'.format(invoker, method))

    if prob_info_c['raw_type'] == 'unconstrained':
        if 'constrviolation' in output.keys():
            del output['constrviolation']
        if 'chist' in output.keys():
            del output['chist']
    elif prob_info_c['refined_type'] == 'unconstrained' and prob_info_c['raw_type'] != 'unconstrained':
        output['constrviolation'] = np.float64(0)
        output['chist'] = np.zeros(nf_c)

    # Revise output['constr_value'] according to problem type.
    if prob_info_c['refined_type'] != 'nonlinearly-constrained' and 'constr_value' in output.keys() and \
            output['constr_value'].size > 0:
        raise ValueError(
            '{}: UNEXPECTED ERROR: {} returns values of nonlinear constraints for a problem that does not admit '
            'such constraints.'.format(invoker, method))

    if prob_info_c['raw_type'] != 'nonlinearly-constrained' and 'constr_value' in output.keys():
        del output['constr_value']

    # Record the returned message.
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
        output['message'] = 'Return from {} because all the variables are fixed by the constraints.'.format(method)
    elif exitflag_c == 14:
        output['message'] = '{} receives a linear feasibility problem and finds a feasible point.'.format(method)
    elif exitflag_c == 15:
        output['message'] = \
            '{} receives a linear feasibility problem but does not find a feasible point.'.format(method)
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
        if np.any(prob_info['infeasible_nonlinear']):
            output['InfeasibleNonlinear'] = np.where(prob_info['infeasible_nonlinear'])[0]

        if np.any(prob_info['infeasible_linear']):
            output['InfeasibleLinear'] = np.where(prob_info['infeasible_linear'])[0]

        if np.any(prob_info['infeasible_bound']):
            output['InfeasibleBound'] = np.where(prob_info['infeasible_bound'])[0]

        output['message'] = 'Return from {} because the constraints are infeasible.'.format(method)
    else:
        raise ValueError('{}: UNEXPECTED ERROR: {} returns an invalid exitflag {}.'.format(invoker, method, exitflag_c))

    # Get the warnings memorized in output.
    if 'warnings' in output.keys():
        warning_list_output = output['warnings']
        del output['warnings']

        if not hasattr(warning_list_output, '__len__'):
            raise SystemError('{}: UNEXPECTED ERROR: the warnings should be defined as a list.'.format(invoker))
    else:
        warning_list_output = []

    # Get the warnings memorized in prob_info.
    if 'warnings' in prob_info_c.keys():
        warning_list = list(prob_info_c['warnings'])
        warning_list.extend(list(warning_list_output))
    else:
        warning_list = []

    # More careful checks about fx, constrviolation, fhist and chist.
    # We do this only if the coe is in debug mode but not in classical mode. The classical mode cannot pass these
    # checks.
    if options['debug'] and not options['classical']:
        if 'raw_data' not in prob_info_keys:
            raise ValueError("{}: UNEXPECTED ERROR: 'raw_data' should be a field of prob_info".format(invoker))

        # Check whether fx is 'optimal'.
        fhistf = fhist_c
        if method in ['bobyqa', 'lincoa', 'cobyla']:
            fhistf = fhistf[chist_c <= max(constrv_returned, 0)]

        if np.isnan(fhistf).all():
            min_f = np.nan
        else:
            min_f = np.nanmin((fx_c, np.nanmin(fhistf)))

        # Tom 2021-05-26: The following test is disabled for lincoa for the moment.
        # if fx != min_f and not (np.isnan(fx) and np.isnan(min_f)) and method != 'lincoa' and \
        #         'constr_modified' in output.keys() and output['constr_modified']:
        if fx != min_f and not (np.isnan(fx) and np.isnan(min_f)) and method != 'lincoa':
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} returns an fhist that does not match nf or fx'.format(invoker, method))

        # Check whether constrviolation is correct.
        cobyla_prec = np.float64(1e-10)
        lincoa_prec = np.float64(1e-12)
        gen_prec = np.float64(1e-12)

        # COBYLA cannot ensure fx=fun(x) or conval=con(x) due to rounding errors. Instead of checking the equality, we
        # check whether the relative error is within cobyla_prec. There can also be a difference between constrviolation
        # and conv, especially if the problem is scaled.
        constrviolation = np.float64(0)
        if 'constrviolation' in output.keys():
            constrviolation = output['constrviolation']
        if method == 'bobyqa' and np.nanmax((constrviolation, np.nanmax(chist_c))) > 0 and \
                not prob_info_c['infeasible'] and not prob_info_c['fixedx']:
            raise ValueError(
                '{}: UNEXPECTED ERROR: {} is a feasible solver yet it returns positive '
                'constrviolations.'.format(invoker, method))

        if (method == 'lincoa' and not output['constr_modified']) or method == 'cobyla':
            linear = prob_info_c['raw_data']['constraints']['linear']
            lb, ub = prob_info_c['raw_data']['bounds']

            # Compute the linear constraint value.
            if linear is not None:
                try:
                    ax = np.dot(linear.A, x_c)
                    r = np.r_[linear.lb - ax, ax - linear.ub]
                except ValueError:
                    raise ValueError(
                        '{}: UNEXPECTED ERROR: the linear constraints are no more consistent'.format(invoker))
            else:
                r = np.asarray([np.nan])

            if method == 'lincoa':
                conv_n = np.concatenate((r, lb - x_c, x_c - ub))
                conv_n = np.nanmax((np.zeros_like(conv_n), conv_n), axis=0)
                conv = np.max(conv_n)
            else:
                nlc = np.asarray([-np.inf], dtype=np.float64)
                if 'constr_value' in output.keys():
                    nlc = np.asarray(output['constr_value'], dtype=np.float64)
                conv = np.concatenate(([0], r, lb - x_c, x_c - ub, nlc))
                conv = np.nanmax(conv)

            if not prob_info_c['infeasible'] and not (np.isnan(conv) and np.isnan(constrviolation_c)) and \
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
            # Check whether fx = fun(x).
            if prob_info_c['raw_data']['objective'] is not None:
                fun_x = prob_info_c['raw_data']['objective'](x_c, *prob_info_c['raw_data']['args'])
            else:
                fun_x = np.float64(0)
            if np.isnan(fun_x) or (fun_x > hugefun):
                fun_x = hugefun
                # Due to extreme barrier (implemented when options['classical']=False), all the function values that are
                # NaN or larger than hugefun are replaced by hugefun.

            # It seems that COBYLA can return fx~=fun(x) due to rounding errors. Moreover, in the general case, the
            # objective function may be noisy, in which case the exact comparison makes no sense.
            if not (np.isnan(fx_c) and np.isnan(fun_x)) and \
                    not (abs(fun_x - fx_c) <= gen_prec * max(1, abs(fx_c)) or
                         (abs(fun_x - fx_c) <= cobyla_prec * max(1, abs(fx_c))) and method == 'cobyla'):
                raise ValueError(
                    '{}: UNEXPECTED ERROR: {} returns an fx that does not match x.'.format(invoker, method))

            # Check whether nlc = nonlinear(x) (true equality).
            nonlinear = prob_info_c['raw_data']['constraints']['nonlinear']
            if nonlinear is not None:
                if 'constr_value' in output.keys():
                    nlc = output['constr_value']
                else:
                    nlc = np.array([], dtype=np.float64)

                if not hasattr(nlc, '__len__'):
                    raise ValueError('{}: UNEXPECTED ERROR: nlc should be recorded as a ndarray.'.format(invoker))

                nlcx = nonlinear['fun'](x_c)
                nlcx[np.logical_or(np.isnan(nlcx), nlcx > hugecon)] = hugecon

                # This part is NOT extreme barrier. We replace extremely negative values of cineq (which leads to no
                # constraint violation) by -hugecon. Otherwise, NaN or Inf may occur in the interpolation models.
                nlcx[nlcx < -hugecon] = -hugecon

                max_x = 0 if nlcx.size == 0 else np.nanmax(nlcx)

                if nlc.size != nlcx.size or not np.array_equal(np.isnan(nlc), np.isnan(nlcx)) or \
                        (np.isnan(nlc).any() and (np.abs(nlc - nlcx) > cobyla_prec * max(1, max_x)).any()):
                    raise ValueError(
                        '{}: UNEXPECTED ERROR: {} returns a con(x) that does not match x.'.format(invoker, method))

    if 'constr_modified' in output.keys():
        del output['constr_modified']

    # Create the 'constr_value' field in the output, with respect to the structure of the constraints in input.
    if 'constr_meta' not in prob_info.keys():
        raise ValueError('{}: UNEXPECTED ERROR: the constraints metadata are not defined.'.format(invoker))

    # If any constraint was provided by the user, the structure constr_value should be added in the output
    if len(prob_info['constr_meta']['linear_indices']) + len(prob_info['constr_meta']['nonlinear_indices']) > 0:
        constr_value = []  # reconstructed list
        k_nonlinear = 0  # index of the current nonlinear constraint in the output array
        try:
            for i_meta, metadata in enumerate(prob_info['constr_meta']['data']):
                if prob_info['infeasible']:
                    # If the problem turned infeasible, the raw constraint values have been recorded, they just need to
                    # be read in order.
                    if i_meta in prob_info['constr_meta']['linear_indices']:
                        constr_value.append(np.dot(metadata['A'], output['x']))
                    else:
                        constr_value.append(output['constr_value'][k_nonlinear:k_nonlinear + metadata['len']])
                        k_nonlinear += metadata['len']
                elif not metadata['trivial']:
                    # The constraint is a non-trivial constraint, which has therefore some components that have been
                    # evaluated.
                    if i_meta in prob_info['constr_meta']['linear_indices']:
                        # The constraint is a linear constraint: we just need to compute the product Ax. Since the
                        # computation of the value of the linear constraint is considered low, we do not built the
                        # global evaluation by using the computation already done by the Fortran code.
                        constr_value.append(np.dot(metadata['A'], output['x']))
                    else:
                        # The current constraint is a nonlinear constraint: since we should absolutely not re-evaluated
                        # the nonlinear constraint function, we decode the values contain in the constraint array,
                        # forwarded by the Fortran code. Note that some evaluations may have been dropped.
                        n_lb = sum(np.logical_not(metadata['dropped_indices_lb']))
                        n_ub = sum(np.logical_not(metadata['dropped_indices_ub']))
                        missing_lb, missing_ub = 0, 0  # the number of missing values so far
                        values = np.full(metadata['len'], np.nan, dtype=np.float64)
                        for i in range(metadata['len']):
                            if not metadata['dropped_indices_lb'][i] and not metadata['dropped_indices_ub'][i]:
                                # If both upper and lower bound are significative, the considered constraint value is
                                # the average of both computed values, to increase the precision.
                                vlb = metadata['lb'][i] - output['constr_value'][k_nonlinear + i - missing_lb]
                                vub = output['constr_value'][k_nonlinear + i + n_lb - missing_ub] + metadata['ub'][i]
                                values[i] = (vlb + vub) / 2
                            elif not metadata['dropped_indices_lb'][i]:
                                # The upper bound of the current component was set to +inf but the lower bound was a
                                # true scalar: we can build the constraint value from it.
                                vlb = metadata['lb'][i] - output['constr_value'][k_nonlinear + i - missing_lb]
                                values[i] = vlb
                                missing_ub += 1  # the upper bound was not significative
                            elif not metadata['dropped_indices_ub'][i]:
                                # The lower bound of the current component was set to -inf but the upper bound was a
                                # true scalar: we can build the constraint value from it.
                                vub = output['constr_value'][k_nonlinear + i + n_lb - missing_ub] + metadata['ub'][i]
                                values[i] = vub
                                missing_lb += 1  # the lower bound was not significative
                            else:
                                # Both lower and upper bound are missing, the value should be set to NaN.
                                missing_lb += 1
                                missing_ub += 1
                        constr_value.append(values)
                        k_nonlinear += n_lb + n_ub
                else:
                    # The constraint has been considered trivial, which led to no evaluation of it. Moreover, it is
                    # necessarily defined as a NonlinearConstraint structure, which provides its length.
                    constr_value.append(np.full(metadata['len'], np.nan, dtype=np.float64))
        except (KeyError, IndexError):
            raise ValueError('{}: UNEXPECTED ERROR: the constraints metadata are ill-defined.'.format(invoker))

        if any(map(lambda a: np.any(np.isnan(a)), constr_value)):
            # The list of constraint values contains some NaN values because some constraints were not considered by the
            # code: the user should be informed.
            w_message = \
                '{}: some constraints components are trivial. They are not evaluated during the computation, and ' \
                'their values are represented by NaN in constr_value.'.format(invoker)
            warnings.warn(w_message, Warning)
            warning_list.append(w_message)
        if any(map(lambda a: a.size == 0, constr_value)):
            w_message = '{}: some constraints are trivial. They are not evaluated during the computation, and they ' \
                        'are represented by empty arrays in constr_value.'.format(invoker)
            warnings.warn(w_message, Warning)
            warning_list.append(w_message)

        if prob_info['constr_meta']['is_list']:
            output['constr_value'] = constr_value
        else:
            output['constr_value'] = constr_value[0]

    # Give back all the warning messages to the user.
    if len(warning_list) > 0:
        output['warnings'] = warning_list

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

    if missing_file is None:
        missing_file = invoker

    raise ImportError('{} is missing. Please reinstall {} and ensure the '
                      'fulfilment of the requirements '
                      '(numpy>=1.20.0).'.format(missing_file, __package__))
