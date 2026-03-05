import prima as p
import numpy as np
import os
import pytest

@pytest.mark.skipif(os.getenv("SKBUILD_WHEEL_CMAKE") != "0",
  reason="We only run this when we skip building bindings")
def test_without_pybindings():
  '''
  This test makes sure we can still run the prima package even if the bindings are
  missing. This is important for the integration of the pure python backend with scipy.
  As an example, this test will catch errors like trying to import the bindings at the
  top level of the package, which would break the pure python backend.

  '''
  result = p.minimize(lambda x: x[0]**2, np.array([5]), method='cobyla', options={'backend': "Python"})
  assert result.x[0] == 0
