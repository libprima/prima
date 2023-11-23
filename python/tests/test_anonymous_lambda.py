# test with anonymous lambda as: objective function, callback function, nonlinear constraint function.
from prima import minimize, NonlinearConstraint as NLC
import numpy as np
import pytest

def test_anonymous_lambda():
    '''
    In previous iterations of these bindings, memory was not handled correctly in C++ and using anonymous lambdas
    would cause issues including, but not limited to, segfaults and infinite hangs. This test is to ensure that
    anonymous lambdas can be used without issue.
    '''
    myNLC = NLC(lambda x: x[0]**2 - 9, [-np.inf], [0])
    res = minimize(lambda x: (x[0] - 5)**2 + (x[1] - 4)**2, [0.0] * 2, method='COBYLA', constraints=myNLC, callback=lambda x, *args: print(x))
    assert abs(res.x[0] - 3) < 1e-2 and abs(res.x[1] - 4) < 1e-2 and abs(res.fun - 4) < 1e-2


def test_anonymous_lambda_unclean_exit():
    '''
    Another potential issue with memory management is when the minimize function throws an exception instead of terminating
    normally. This can be triggered by providing an invalid method. Another option might be to raise an exception within any
    of the anonymous lambdas.
    '''
    myNLC = NLC(lambda x: x[0]**2 - 9, [-np.inf], [0])
    with pytest.raises(ValueError):
        minimize(lambda x: (x[0] - 5)**2 + (x[1] - 4)**2, [0.0] * 2, method='GARBAGE', constraints=myNLC, callback=lambda x, *args: print(x))
    
