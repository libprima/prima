import ctypes as ct
import ctypes.util
from dataclasses import dataclass
import numpy as np
import os
import sys


_lib_name = "primac"
try:
    # import from current dir
    import numpy as np

    _path = os.path.abspath(os.path.dirname(__file__))
    _handle = np.ctypeslib.load_library(_lib_name, _path)
except Exception:
    # else rely on loader env variables
    _lib_loc = ctypes.util.find_library(_lib_name)
    if _lib_loc is None:
        raise ValueError(f"Cannot find {_lib_name} library")
    _handle = ct.cdll.LoadLibrary(_lib_loc)

# primac callback form
_prima_function_callback = ct.CFUNCTYPE(
    None, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_void_p
)
_prima_function_con_callback = ct.CFUNCTYPE(
    None, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double), ct.c_void_p
)

# prima_get_rc_string returns char*
_handle.prima_get_rc_string.restype = ct.c_char_p


def _prima_pure_callback(f_obj, n, args):
    """ Convert a Python objective function into a C function."""
    def _inner(c_x, c_f, c_data):
        x = [c_x[i] for i in range(n)]
        try:
            f = f_obj(x, *args)
        except Exception:
            # prima treats nan as error indicator
            f = float("nan")
        c_f[0] = f

    return _prima_function_callback(_inner)


def _prima_pure_con_callback(f_obj, n, f_con, m, args):
    """ Convert Python objective/constraints functions into a C function."""
    def _inner(c_x, c_f, c_con, c_data):
        x = [c_x[i] for i in range(n)]
        try:
            f = f_obj(x, *args)
            con = f_con(x, *args)
        except Exception:
            # prima treats nan as error indicator
            f = float("nan")
            con = [float("nan")] * m
        c_f[0] = f
        for i in range(m):
            c_con[i] = con[i]

    return _prima_function_con_callback(_inner)


@dataclass
class Result:
    """Optimization result."""
    x: list
    success: bool
    status: int
    message: str
    fun: float
    n_fev: int


def _convert_bounds(bounds, n):
    """
    xl, xu : sequence or float
    n : int
    """
    if bounds is None:
        xl = [-sys.float_info.max] * n
        xu = [sys.float_info.max] * n
    else:
        xl = bounds.lb
        xu = bounds.ub
    assert len(xl) == n
    assert len(xu) == n
    c_xl = (ct.c_double * n)(*xl)
    c_xu = (ct.c_double * n)(*xu)
    return c_xl, c_xu


def bobyqa(fun, x0, args=(), bounds=None,
           rhobeg=1.0, rhoend=1e-3,
           ftarget=float('-inf'), maxfun=1000, iprint=0):
    """
    BOBYQA solver.

    Parameters
    ----------
    fun : callable
        objective function, returns a float
    x0 : sequence of float
        starting point
    args : tuple, optional
        Extra arguments passed to the objective function
    bounds, scipy.optimize.Bounds
        Bounds on variables
    rhobeg, rhoend : float, optional
        trust region radius
    ftarget : float, optional
        function value target
    maxfun : int, optional
        max function evaluation
    iprint : int, optional
        verbosity level

    Returns
    -------
    res : Result
        Result class
    """

    # convert input arguments into C types
    n = len(x0)
    cn = ct.c_int(n)
    c_x = (ct.c_double * n)(*x0)
    c_f = ct.pointer(ct.c_double())
    c_xl, c_xu = _convert_bounds(bounds, n)
    c_nf = ct.pointer(ct.c_int())
    c_rhobeg = ct.c_double(rhobeg)
    c_rhoend = ct.c_double(rhoend)
    c_ftarget = ct.c_double(ftarget)
    c_maxfun = ct.c_int(maxfun)
    c_npt = ct.c_int(2*n+1)
    c_iprint = ct.c_int(iprint)
    c_fun = _prima_pure_callback(fun, n, args)
    c_data = ct.c_void_p(0)

    # call the PRIMA library
    status = _handle.prima_bobyqa(c_fun, c_data, cn, c_x, c_f, c_xl, c_xu,
                                  c_nf, c_rhobeg, c_rhoend,
                                  c_ftarget, c_maxfun, c_npt, c_iprint)

    # convert back output arguments
    x = [c_x[i] for i in range(n)]
    success = (status == 0)
    c_status = ct.c_int(status)
    message = _handle.prima_get_rc_string(c_status).decode()
    n_fev = c_nf[0]

    return Result(x, success, status, message, c_f[0], n_fev)


def cobyla(fun, x0, f_con, m_nlcon, args=(),
           Aineq=[], bineq=[], Aeq=[], beq=[], bounds=None,
           rhobeg=1.0, rhoend=1e-3,
           ftarget=float('-inf'), maxfun=1000, iprint=0):
    """
    COBYLA solver.

    Parameters
    ----------
    fun : callable
        objective function, returns a float
    x0 : sequence of float
        starting point
    args : tuple, optional
        Extra arguments passed to the objective function
    f_con : callable
         return the constraint values as a sequence of float of size m_nlcon
    m_nlcon : int
        non-linear constraint dimension
    Aineq, bineq : 2-d/1-d arrays of float
        Inequality constraints Aineq*x<=bineq
    Aeq, beq : 2-d/1-d arrays of float
        Equality constraints Aeq*x=beq
    bounds, scipy.optimize.Bounds
        Bounds on variables
    rhobeg, rhoend : float, optional
        trust region radius
    ftarget : float, optional
        function value target
    maxfun : int, optional
        max function evaluation
    iprint : int, optional
        verbosity level

    Returns
    -------
    res : Result
        Result class
    """

    # convert input arguments into C types
    n = len(x0)
    c_m_nlcon = ct.c_int(m_nlcon)
    c_fun_con = _prima_pure_con_callback(fun, n, f_con, m_nlcon, args)
    c_data = ct.c_void_p(0)
    cn = ct.c_int(n)
    c_x = (ct.c_double * n)(*x0)
    c_f = ct.pointer(ct.c_double())
    c_cstrv = ct.pointer(ct.c_double())
    nl = [0.0] * m_nlcon
    c_nlconstr = (ct.c_double * m_nlcon)(*nl)
    bineq = np.atleast_1d(bineq).flatten()
    m_ineq = len(bineq)
    cm_ineq = ct.c_int(m_ineq)
    Aineq = np.atleast_2d(Aineq).flatten()
    assert len(Aineq) == n * m_ineq
    c_Aineq = (ct.c_double * (m_ineq * n))(*Aineq)
    c_bineq = (ct.c_double * n)(*bineq)
    beq = np.atleast_1d(beq).flatten()
    m_eq = len(beq)
    cm_eq = ct.c_int(m_eq)
    Aeq = np.atleast_2d(Aeq).flatten()
    assert len(Aeq) == n * m_eq
    c_Aeq = (ct.c_double * (m_eq * n))(*Aeq)
    c_beq = (ct.c_double * n)(*beq)
    c_xl, c_xu = _convert_bounds(bounds, n)
    c_nf = ct.pointer(ct.c_int())
    c_rhobeg = ct.c_double(rhobeg)
    c_rhoend = ct.c_double(rhoend)
    c_ftarget = ct.c_double(ftarget)
    c_maxfun = ct.c_int(maxfun)
    c_iprint = ct.c_int(iprint)

    # call the PRIMA library
    status = _handle.prima_cobyla(c_m_nlcon, c_fun_con, c_data, cn, c_x, c_f,
                                  c_cstrv, c_nlconstr,
                                  cm_ineq, c_Aineq, c_bineq,
                                  cm_eq, c_Aeq, c_beq,
                                  c_xl, c_xu,
                                  c_nf, c_rhobeg, c_rhoend,
                                  c_ftarget, c_maxfun, c_iprint)

    # convert back output arguments
    x = [c_x[i] for i in range(n)]
    n_fev = c_nf[0]
    success = (status == 0)
    c_status = ct.c_int(status)
    message = _handle.prima_get_rc_string(c_status).decode()
    return Result(x, success, status, message, c_f[0], n_fev)


def newuoa(fun, x0, args=(), rhobeg=1.0, rhoend=1e-3,
           ftarget=float('-inf'), maxfun=1000, iprint=0):
    """
    NEWUOA solver.

    Parameters
    ----------
    fun : callable
        objective function, returns a float
    x0 : sequence of float
        starting point
    args : tuple, optional
        Extra arguments passed to the objective function
    rhobeg, rhoend : float, optional
        trust region radius
    ftarget : float, optional
        function value target
    maxfun : int, optional
        max function evaluation
    iprint : int, optional
        verbosity level

    Returns
    -------
    res : Result
        Result class
    """

    # convert input arguments into C types
    n = len(x0)
    cn = ct.c_int(n)
    c_x = (ct.c_double * n)(*x0)
    c_f = ct.pointer(ct.c_double())
    c_nf = ct.pointer(ct.c_int())
    c_rhobeg = ct.c_double(rhobeg)
    c_rhoend = ct.c_double(rhoend)
    c_ftarget = ct.c_double(ftarget)
    c_maxfun = ct.c_int(maxfun)
    c_npt = ct.c_int(2*n+1)
    c_iprint = ct.c_int(iprint)
    c_fun = _prima_pure_callback(fun, n, args)
    c_data = ct.c_void_p(0)

    # call the PRIMA library
    status = _handle.prima_newuoa(c_fun, c_data, cn, c_x, c_f,
                                  c_nf, c_rhobeg, c_rhoend,
                                  c_ftarget, c_maxfun, c_npt, c_iprint)

    # convert back output arguments
    x = [c_x[i] for i in range(n)]
    n_fev = c_nf[0]
    success = (status == 0)
    c_status = ct.c_int(status)
    message = _handle.prima_get_rc_string(c_status).decode()
    return Result(x, success, status, message, c_f[0], n_fev)


def uobyqa(fun, x0, args=(), rhobeg=1.0, rhoend=1e-3,
           ftarget=float('-inf'), maxfun=1000, iprint=0):
    """
    UOBYQA solver.

    Parameters
    ----------
    fun : callable
        objective function, returns a float
    x0 : sequence of float
        starting point
    args : tuple, optional
        Extra arguments passed to the objective function
    rhobeg, rhoend : float, optional
        trust region radius
    ftarget : float, optional
        function value target
    maxfun : int, optional
        max function evaluation
    iprint : int, optional
        verbosity level

    Returns
    -------
    res : Result
        Result class
    """

    # convert input arguments into C types
    n = len(x0)
    cn = ct.c_int(n)
    c_x = (ct.c_double * n)(*x0)
    c_f = ct.pointer(ct.c_double())
    c_nf = ct.pointer(ct.c_int())
    c_rhobeg = ct.c_double(rhobeg)
    c_rhoend = ct.c_double(rhoend)
    c_ftarget = ct.c_double(ftarget)
    c_maxfun = ct.c_int(maxfun)
    c_iprint = ct.c_int(iprint)
    c_fun = _prima_pure_callback(fun, n, args)
    c_data = ct.c_void_p(0)

    # call the PRIMA library
    status = _handle.prima_uobyqa(c_fun, c_data, cn, c_x, c_f,
                                  c_nf, c_rhobeg, c_rhoend,
                                  c_ftarget, c_maxfun, c_iprint)

    # convert back output arguments
    x = [c_x[i] for i in range(n)]
    n_fev = c_nf[0]
    success = (status == 0)
    c_status = ct.c_int(status)
    message = _handle.prima_get_rc_string(c_status).decode()
    return Result(x, success, status, message, c_f[0], n_fev)


def lincoa(fun, x0, args=(), Aineq=[], bineq=[], Aeq=[], beq=[],
           bounds=None, rhobeg=1.0, rhoend=1e-3,
           ftarget=float('-inf'), maxfun=1000, iprint=0):
    """
    LINCOA solver.

    Parameters
    ----------
    fun : callable
        objective function, returns a float
    x0 : sequence of float
        starting point
    args : tuple, optional
        Extra arguments passed to the objective function
    Aineq, bineq : 2-d/1-d arrays of float
        Inequality constraints Aineq*x<=bineq
    Aeq, beq : 2-d/1-d arrays of float
        Equality constraints Aeq*x=beq
    bounds, scipy.optimize.Bounds
        Bounds on variables
    rhobeg, rhoend : float, optional
        trust region radius
    ftarget : float, optional
        function value target
    maxfun : int, optional
        max function evaluation
    iprint : int, optional
        verbosity level

    Returns
    -------
    res : Result
        Result class
    """

    # convert input arguments into C types
    n = len(x0)
    cn = ct.c_int(n)
    c_x = (ct.c_double * n)(*x0)
    c_f = ct.pointer(ct.c_double())
    c_cstrv = ct.pointer(ct.c_double())
    bineq = np.atleast_1d(bineq).flatten()
    m_ineq = len(bineq)
    cm_ineq = ct.c_int(m_ineq)
    Aineq = np.atleast_2d(Aineq).flatten()
    assert len(Aineq) == n * m_ineq
    c_Aineq = (ct.c_double * (m_ineq * n))(*Aineq)
    c_bineq = (ct.c_double * m_ineq)(*bineq)
    beq = np.atleast_1d(beq).flatten()
    m_eq = len(beq)
    cm_eq = ct.c_int(m_eq)
    Aeq = np.atleast_2d(Aeq).flatten()
    assert len(Aeq) == n * m_eq
    c_Aeq = (ct.c_double * (m_eq * n))(*Aeq)
    c_beq = (ct.c_double * n)(*beq)
    c_xl, c_xu = _convert_bounds(bounds, n)
    c_nf = ct.pointer(ct.c_int())
    c_rhobeg = ct.c_double(rhobeg)
    c_rhoend = ct.c_double(rhoend)
    c_ftarget = ct.c_double(ftarget)
    c_maxfun = ct.c_int(maxfun)
    c_npt = ct.c_int(2*n+1)
    c_iprint = ct.c_int(iprint)
    c_fun = _prima_pure_callback(fun, n, args)
    c_data = ct.c_void_p(0)

    # call the PRIMA library
    status = _handle.prima_lincoa(c_fun, c_data, cn, c_x, c_f,
                                  c_cstrv,
                                  cm_ineq, c_Aineq, c_bineq,
                                  cm_eq, c_Aeq, c_beq,
                                  c_xl, c_xu,
                                  c_nf, c_rhobeg, c_rhoend,
                                  c_ftarget, c_maxfun, c_npt, c_iprint)

    # convert back output arguments
    x = [c_x[i] for i in range(n)]
    n_fev = c_nf[0]
    success = (status == 0)
    c_status = ct.c_int(status)
    message = _handle.prima_get_rc_string(c_status).decode()
    return Result(x, success, status, message, c_f[0], n_fev)
