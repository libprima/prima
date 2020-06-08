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

from fnmatch import fnmatch

import glob
import platform
import shutil

import re
from os import listdir, remove, walk
from os.path import dirname, abspath, join, relpath

from numpy.distutils.core import setup, Extension

# Set the paths to all the folders that will be used to build PDFO
CURRENT_WD = dirname(abspath(__file__))
PDFO_WD = join(CURRENT_WD, 'python')
FSRC_WD = join(dirname(PDFO_WD), 'fsrc')
FSRC_CLASSICAL_WD = join(FSRC_WD, 'classical')
GATEWAYS_WD = join(PDFO_WD, 'py_gateways')
GATEWAYS_CLASSICAL_WD = join(GATEWAYS_WD, 'classical')
INTERFACES_WD = join(PDFO_WD, 'interfaces', 'pdfo')

# Set the options that will be given to F2PY to build PDFO
OPTIONS = ['--quiet']


def clean():
    """Clean up the arborescence of the package."""
    folders_to_delete = ['build', 'develop-eggs', 'dist', 'eggs', '.eggs', 'sdist', 'wheels']
    if platform.system() == 'windows':
        files_to_delete = glob.glob('**/*.pyd')  # Windows-based binary files
    else:
        files_to_delete = glob.glob('**/*.so')  # Unix-based binary files
    files_to_delete.extend(glob.glob('*.mod') + glob.glob('**/*.mod'))  # Fortran module files
    files_to_delete.extend(glob.glob('**/*.o'))  # C and Fortran object files
    for root, dirs, files in walk(CURRENT_WD):
        for folder in dirs:
            if fnmatch(folder, '*.egg-info'):
                folders_to_delete.append(join(root, folder))
        for file in files:
            if fnmatch(file, '*.egg'):
                files_to_delete.append(join(root, file))

    for folder in folders_to_delete:
        shutil.rmtree(folder, ignore_errors=True)
    for tmp_file in files_to_delete:
        try:
            remove(tmp_file)
        except FileNotFoundError:
            # This exception should never occur, except if a user delete a file manually after running this script and
            # before the execution of the remove command. Since it is almost impossible, the exception is just caught.
            pass


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
    clean()

    setup(
        name='pdfo',
        version=open('VERSION.txt').read().rstrip(),
        description="Powell's Derivative-Free Optimization solvers",
        long_description=open('README.txt').read(),
        long_description_content_type='text/plain',
        author='Tom M. Ragonneau and Zaikun Zhang',
        author_email='pdfocode@gmail.com',
        url='https://www.pdfo.net',
        packages=['pdfo', 'pdfo.tests'],
        package_dir={
            'pdfo': relpath(INTERFACES_WD),
            'pdfo.tests': relpath(join(INTERFACES_WD, 'tests')),
        },
        include_package_data=True,
        ext_modules=EXT_MODULES,
        license='GNU Lesser General Public License v3 or later (LGPLv3+)',
        keywords='Powell Derivative-Free Optimization Fortran UOBYQA NEWUOA BOBYQA LINCOA COBYLA',
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Framework :: IDLE',
            'Framework :: IPython',
            'Framework :: Jupyter',
            'Framework :: Pytest',
            'Intended Audience :: Developers',
            'Intended Audience :: Education',
            'Intended Audience :: End Users/Desktop',
            'Intended Audience :: Information Technology',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS',
            'Operating System :: MacOS :: MacOS 9',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: Microsoft :: Windows :: Windows 10',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Unix',
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
