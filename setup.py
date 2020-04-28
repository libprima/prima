#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PDFO: Powell's Derivative-Free Optimization solvers.

Authors
-------
Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
Department of Applied Mathematics,
The Hong Kong Polytechnic University.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
"""
import setuptools

import re
from os import listdir
from os.path import dirname, abspath, join, relpath

try:
    from numpy.distutils.core import setup, Extension
except ModuleNotFoundError:
    print('NumPy package is missing.')
    print('Please install a version greater than or equal to 1.10.0.\n')
    raise

# Set the paths to all the folders that will be used to build PDFO
CURRENT_WD = dirname(abspath(__file__))
PDFO_WD = join(CURRENT_WD, 'python')
FSRC_WD = join(dirname(PDFO_WD), 'fsrc')
FSRC_CLASSICAL_WD = join(FSRC_WD, 'classical')
GATEWAYS_WD = join(PDFO_WD, 'py_gateways')
GATEWAYS_CLASSICAL_WD = join(GATEWAYS_WD, 'classical')
INTERFACES_WD = join(PDFO_WD, 'interfaces')

# Set the options that will be given to F2PY to build PDFO
OPTIONS = ['--quiet']


def build_solver(solver):
    """Generates the extensions that should be compiled, i.e. Powell's code."""
    # Default version: this extension corresponds to the Fortran code modified by Tom M. RAGONNEAU and Zaikun ZHANG.
    sources = [
        relpath(join(GATEWAYS_WD, '{}-interface.pyf'.format(solver))),
        relpath(join(GATEWAYS_WD, '{}.f90'.format(solver))),
    ]
    solver_src = relpath(join(FSRC_WD, solver))
    f77_extension = re.compile('.*\.f$')
    sources.extend([join(solver_src, fortran) for fortran in filter(f77_extension.match, listdir(solver_src))])
    ext_solver = Extension(name='pdfo.f{}'.format(solver), sources=sources, f2py_options=OPTIONS)

    # Classical version: this extension corresponds to the original Fortran code written by late Prof. Powell. They are
    # called when the user provides the option classical=True to any of the Python interface.
    sources = [
        relpath(join(GATEWAYS_CLASSICAL_WD, '{}-interface.pyf'.format(solver))),
        relpath(join(GATEWAYS_CLASSICAL_WD, '{}.f90'.format(solver))),
    ]
    solver_src = relpath(join(FSRC_CLASSICAL_WD, solver))
    f77_extension = re.compile('.*\.f$')
    sources.extend([join(solver_src, fortran) for fortran in filter(f77_extension.match, listdir(solver_src))])
    ext_solver_classical = Extension(name='pdfo.f{}_classical'.format(solver), sources=sources, f2py_options=OPTIONS)

    return ext_solver, ext_solver_classical


# Define the PDFO's extensions
#    * PDFOCONST holds all the general constants.
#    * GETHUGE gives the maximal values that the Fortran code can handle for each type.
EXT_MODULES = [
    Extension(
        name='pdfo.pdfoconst',
        sources=[
            relpath(join(GATEWAYS_WD, 'pdfoconst-interface.pyf')),
            relpath(join(FSRC_WD, 'pdfoconst.F'))
        ],
        f2py_options=OPTIONS),
    Extension(
        name='pdfo.gethuge',
        sources=[
            relpath(join(GATEWAYS_WD, 'gethuge-interface.pyf')),
            relpath(join(GATEWAYS_WD, 'gethuge.f90'))
        ],
        f2py_options=OPTIONS),
]
for algorithm in ['uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla']:
    # Generate the default and classical extensions to the Powell's Fortran code
    EXT_MODULES.extend(build_solver(algorithm))


if __name__ == '__main__':
    setup(
        name='pdfo',
        version=open('VERSION.txt').read().rstrip(),
        description="Powell's Derivative-Free Optimization solvers",
        long_description=open('README.txt').read(),
        long_description_content_type='text/plain',
        author='Tom M. Ragonneau and Zaikun Zhang',
        author_email='pdfocode@gmail.com',
        url='https://www.pdfo.net',
        packages=['pdfo'],
        package_dir={'': INTERFACES_WD},
        include_package_data=True,
        ext_modules=EXT_MODULES,
        license='GNU Lesser General Public License v3 or later (LGPLv3+)',
        keywords='Fortran Powell DFO optimization',
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Framework :: IDLE',
            'Framework :: IPython',
            'Framework :: Jupyter',
            'Framework :: Pytest',
            'Framework :: tox',
            'Intended Audience :: Developers',
            'Intended Audience :: Education',
            'Intended Audience :: End Users/Desktop',
            'Intended Audience :: Information Technology',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
            'Natural Language :: English',
            'Operating System :: MacOS',
            'Operating System :: MacOS :: MacOS 9',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: Microsoft :: Windows :: Windows 10',
            'Operating System :: Microsoft :: Windows :: Windows 7',
            'Operating System :: Microsoft :: Windows :: Windows 8',
            'Operating System :: Microsoft :: Windows :: Windows 8.1',
            'Operating System :: OS Independent',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Unix',
            'Programming Language :: C',
            'Programming Language :: Fortran',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.0',
            'Programming Language :: Python :: 3.1',
            'Programming Language :: Python :: 3.2',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: Implementation',
            'Topic :: Education',
            'Topic :: Scientific/Engineering',
            'Topic :: Software Development :: Libraries',
            'Topic :: Utilities',
        ],
        install_requires=['numpy'],
        setup_requires=['pytest-runner'],
        tests_require=['pytest'],
        test_suite='pdfo.tests',
        python_requires='>=2.7',
        zip_safe=True,
    )
