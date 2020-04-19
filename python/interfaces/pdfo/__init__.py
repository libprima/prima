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

from numpy import set_printoptions

from ._dependencies import OptimizeResult, Bounds, LinearConstraint, NonlinearConstraint

from ._bobyqa import bobyqa
from ._cobyla import cobyla
from ._lincoa import lincoa
from ._newuoa import newuoa
from ._uobyqa import uobyqa
from ._pdfo import pdfo

# Automatic truncation of the arrays after more than 10 iterations in the printing statements (for fhist and chist)
set_printoptions(threshold=10)

__author__ = 'Tom M. Ragonneau and Zaikun Zhang'
if datetime.now().year == 2020:
    __copyright__ = 'Copyright {}, Tom M. Ragonneau and Zaikun Zhang'.format(datetime.now().year)
else:
    __copyright__ = 'Copyright 2020--{}, Tom M. Ragonneau and Zaikun Zhang'.format(datetime.now().year)
__credits__ = ['Tom M. Ragonneau', 'Zaikun Zhang']
__license__ = 'GPLv3'
__version__ = '0.9'
__date__ = 'March, 2020'
__maintainer__ = 'Tom M. Ragonneau and Zaikun Zhang'
__email__ = 'tom.ragonneau@connect.polyu.hk and zaikun.zhang@polyu.edu.hk'
__status__ = 'Prototype'  # Production
