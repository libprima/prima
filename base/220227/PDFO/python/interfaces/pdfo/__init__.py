# -*- coding: utf-8 -*-
"""Management of the importable functions of pdfo.

Authors
-------
Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
Department of Applied Mathematics,
The Hong Kong Polytechnic University.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
"""
from __future__ import division, print_function, absolute_import

from datetime import datetime

from ._dependencies import OptimizeResult, Bounds, LinearConstraint, NonlinearConstraint

from ._bobyqa import bobyqa
from ._cobyla import cobyla
from ._lincoa import lincoa
from ._newuoa import newuoa
from ._uobyqa import uobyqa
from ._pdfo import pdfo
from . import tests
from .tests import test_pdfo as testpdfo

# Definition of the metadata of PDFO for Python. It is accessible via:
# >>> import pdfo
# >>> print(pdfo.__author__)
# >>> ...
__all__ = ['OptimizeResult', 'Bounds', 'LinearConstraint', 'NonlinearConstraint', 'bobyqa', 'cobyla', 'lincoa',
           'newuoa', 'uobyqa', 'pdfo', 'tests', 'testpdfo']
__author__ = 'Tom M. Ragonneau and Zaikun Zhang'
__copyright__ = 'Copyright 2020--{}, Tom M. Ragonneau and Zaikun Zhang'.format(datetime.now().year)
__credits__ = ['Tom M. Ragonneau', 'Zaikun Zhang', 'Antoine Dechaume']
__license__ = 'LGPLv3+'
__version__ = '1.1'
__date__ = 'August, 2021'
__maintainer__ = 'Tom M. Ragonneau and Zaikun Zhang'
__email__ = 'tom.ragonneau@connect.polyu.hk and zaikun.zhang@polyu.edu.hk'
__status__ = 'Production'
