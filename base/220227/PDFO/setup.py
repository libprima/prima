#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PDFO - Powell's Derivative-Free Optimization solvers

PDFO (Powell's Derivative-Free Optimization solvers) is a cross-platform package providing interfaces for using late
Professor M. J. D. Powell's derivative-free optimization solvers, including UOBYQA, NEWUOA, BOBYQA, LINCOA, and COBYLA.

See https://www.pdfo.net for more information.
"""
import setuptools  # noqa

from fnmatch import fnmatch

import glob
import platform
import shutil

import re
import sys
from os import listdir, remove, walk
from os.path import dirname, abspath, join, relpath

if sys.version_info < (3, 7):
    raise RuntimeError('Python version >= 3.7 required.')

try:
    from numpy.distutils.core import setup, Extension
except:
    raise Exception('\nPlease install NumPy before installing PDFO.\n')

if sys.platform == "win32":
    # Fix build with gcc under windows.
    # See https://github.com/jameskermode/f90wrap/issues/96
    from numpy.f2py.cfuncs import includes0
    includes0["setjmp.h"] = '#include <setjmpex.h>'

# Set the paths to all the folders that will be used to build PDFO.
CURRENT_WD = dirname(abspath(__file__))
PDFO_WD = join(CURRENT_WD, 'python')
FSRC_WD = join(dirname(PDFO_WD), 'fsrc')
FSRC_CLASSICAL_WD = join(FSRC_WD, 'classical')
GATEWAYS_WD = join(PDFO_WD, 'py_gateways')
GATEWAYS_CLASSICAL_WD = join(GATEWAYS_WD, 'classical')
INTERFACES_WD = join(PDFO_WD, 'interfaces', 'pdfo')

# Set the options that will be given to F2PY to build PDFO.
OPTIONS = ['--quiet']

# Set the descriptions of the package.
DOCLINES = (__doc__ or '').split('\n')


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
        except:
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
        description=DOCLINES[0],
        long_description='\n'.join(DOCLINES[2:]),
        long_description_content_type='text/plain',
        author='Tom M. Ragonneau and Zaikun Zhang',
        author_email='pdfocode@gmail.com',
        maintainer='Tom M. Ragonneau and Zaikun Zhang',
        maintainer_email='pdfocode@gmail.com',
        url='https://www.pdfo.net',
        download_url='https://www.pdfo.net/docs.html#releases',
        project_urls={
            'Bug Tracker': 'https://github.com/pdfo/pdfo/issues',
            'Documentation': 'https://www.pdfo.net',
            'Source Code': 'https://github.com/pdfo/pdfo',
        },
        packages=['pdfo', 'pdfo.tests'],
        package_dir={
            'pdfo': relpath(INTERFACES_WD),
            'pdfo.tests': relpath(join(INTERFACES_WD, 'tests')),
        },
        include_package_data=True,
        ext_modules=EXT_MODULES,
        license='GNU Lesser General Public License v3 or later (LGPLv3+)',
        keywords='Powell Derivative-Free Optimization UOBYQA NEWUOA BOBYQA LINCOA COBYLA',
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Developers',
            'Intended Audience :: Education',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
            'Operating System :: Unix',
            'Operating System :: MacOS',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Microsoft :: Windows',
            'Programming Language :: Fortran',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Software Development :: Libraries',
            'Topic :: Software Development :: Libraries :: Python Modules',
        ],
        install_requires=['numpy>=1.20.0'],
        setup_requires=['pytest-runner'],
        tests_require=['pytest'],
        test_suite='pdfo.tests',
        python_requires='>=3.7',
        zip_safe=True,
    )
