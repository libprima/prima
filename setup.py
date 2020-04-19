#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compilation of the pdfo package.

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
import os
import platform
import sys
from os.path import abspath, dirname, expanduser, isfile, join, relpath
from os import listdir, remove, rename
import re
from subprocess import Popen, PIPE

python_cmd = sys.executable  # current Python interpreter
c_pwd = dirname(abspath(__file__))
pdfo_pwd = join(c_pwd, 'python')
fsrc = join(dirname(pdfo_pwd), 'fsrc')
fsrc_classical = join(fsrc, 'classical')
gateways = join(pdfo_pwd, 'py_gateways')
gateways_classical = join(gateways, 'classical')
interfaces_path = join(pdfo_pwd, 'interfaces')
interfaces = join(interfaces_path, 'pdfo')

# define the current os
system_os = platform.system().lower()
system_known = True
if system_os in ['darwin', 'linux']:
    system_os = 'unix'
elif system_os == 'windows':
    system_os = 'win'
else:
    system_known = False

# define the current version of Python
python_version = sys.version_info[0]
python_version_minor = sys.version_info[1]


def remove_temporary_compilation(folder, include_bin=False):
    """Removes the compilation temporary files.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    assert isinstance(folder, str)

    if include_bin:
        removed_extensions = ['mod', 'o', 'pyd', 'so']
    else:
        removed_extensions = ['mod', 'o']

    for extension in removed_extensions:
        filter_extension = re.compile('.*\.{}$'.format(extension))
        for del_files in filter(filter_extension.match, listdir(folder)):
            remove(join(folder, del_files))


def system_stdout(commands):
    """Executes the compilation command in argument.

    Authors
    -------
    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    Department of Applied Mathematics,
    The Hong Kong Polytechnic University.

    Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    """
    assert isinstance(commands, list)

    p = Popen(commands, stdout=PIPE, stderr=PIPE)
    _, error = p.communicate()

    if p.returncode != 0:
        raise OSError(error.decode('utf-8'))


if __name__ == '__main__':
    # define which solvers should be compiled
    all_solvers = ['uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla']
    if len(sys.argv) > 1:
        solvers = sys.argv[1:]
        unknown = []
        for solver in solvers:
            if solver not in all_solvers:
                unknown.append(solver)

        if len(unknown) == 1:
            raise NameError('The following solver is unknown and cannot be compiled: {}.'.format(unknown[0]))
        elif len(unknown) > 1:
            raise NameError('The following solvers are unknown and cannot be compiled: {}.'.format(', '.join(unknown)))

        solvers = list(dict.fromkeys(solvers))
    else:
        solvers = all_solvers

    # check the >=2.7 version requirement
    success = True
    python_success = True
    if python_version < 2 or (python_version == 2 and python_version_minor < 7):
        print('\nThe installed version of Python is not compatible with PDFO.')
        print('Please install a version greater than or equal to 2.7.', end='\n\n')
        success = False
        python_success = False

    # check the 32-bit requirement on Windows
    if success:
        if system_os == 'win' and platform.architecture()[0].lower() == '64bit':
            print('\nThe current version of PDFO does not support 64-bit Python on Windows.')
            success = False
            python_success = False

    # check the NumPy requirements
    if success:
        print("\nVerifying the set-up of NumPy and F2PY for '{}' ...".format(python_cmd))
        try:
            import numpy as np
        except:
            print('\nNumPy package is missing.')
            print('Please install a version greater than or equal to 1.10.0.')
            success = False

    if success:
        numpy_version = np.__version__.split('.')
        if len(numpy_version) < 2 or int(numpy_version[0]) < 1 or int(numpy_version[1]) < 10:
            print('\nThe installed version of NumPy is not compatible with PDFO.')
            print('Please install a version greater than or equal to 1.10.0.')
            success = False

    if success:
        try:
            import numpy.f2py
        except:
            print('\nThe NumPy module F2PY cannot be found.')
            print('Please check the installation of your NumPy.')
            success = False
    
    if success:
        # compiler command
        compiler = '{} -m numpy.f2py -c --quiet -I{}'.format(python_cmd, relpath(fsrc)).split(' ')

        # if the compilation fails on Windows, the compiler will try to precise the expected compiler, i.e., MinGW
        compiler_mingw32 = '{} -m numpy.f2py -c --compiler=mingw32 --quiet -I{}'.format(python_cmd, relpath(fsrc)).split(' ')
    else:
        compiler = ''
        compiler_mingw32 = ''

    if success:
        # check that a fortran compiler has been installed and configured
        remove_temporary_compilation(c_pwd)
        remove_temporary_compilation(pdfo_pwd)
        remove_temporary_compilation(fsrc)
        remove_temporary_compilation(fsrc_classical)
        remove_temporary_compilation(gateways)
        remove_temporary_compilation(gateways_classical)
        signature = relpath(join(gateways, 'pdfoconst-interface.pyf'))
        fortran = relpath(join(fsrc, 'pdfoconst.F'))
        try:
            system_stdout(compiler + [signature, fortran])
        except OSError:
            remove_temporary_compilation(c_pwd, include_bin=True)
            if system_os == 'win':
                try:
                    compiler = compiler_mingw32
                    system_stdout(compiler + [signature, fortran])
                except OSError:
                    success = False
            else:
                success = False

        if success:
            try:
                from pdfoconst import use_pdfoconst

                # check (very roughly) whether use_pdfoconst calculates 2*1 correctly
                if not (1.9 < use_pdfoconst(1) < 2.1):
                    success = False
            except:
                success = False

        if success:
            print('NumPy and F2PY are correctly set up.', end='\n\n')
        else:
            remove_temporary_compilation(c_pwd, include_bin=True)
            print('\nNumPy and F2PY are installed but F2PY is not ready to compile Fortran.')
            if system_os == 'win':
                print('It may be because F2PY cannot find gcc-fortran.')
            else:
                print('It may be because F2PY cannot find gfortran.')
            if system_known:
                print('Please see `README_py_{}.txt`.'.format(system_os))
            else:
                print('Please see the README file corresponding to your system.')

    if python_success and not success:
        print('\nVerification failed.', end='\n\n')
        print('Your Python is not properly configured for compiling Fortran using F2PY.')
        if system_os == 'win':
            print('Please check the installation and configuration of Python, NumPy, F2PY,\nand gcc-fortran before using this package.', end='\n\n')
        else:
            print('Please check the installation and configuration of Python, NumPy, F2PY,\nand gfortran before using this package.', end='\n\n')
    elif python_success:
        # start the compilation
        print('Compilation starts. It may take some time ...')

        try:
            # compile gethuge
            signature = relpath(join(gateways, 'gethuge-interface.pyf'))
            fortran = relpath(join(gateways, 'gethuge.f90'))
            system_stdout(compiler + [signature, fortran])

            # compile the fortran source files
            for algorithm in solvers:
                print('Compiling {} ...'.format(algorithm))

                # compilation of the default version
                signature = relpath(join(gateways, '{}-interface.pyf'.format(algorithm)))
                module = relpath(join(gateways, '{}.f90'.format(algorithm)))
                command = compiler + [signature, module]

                sub_fsrc = relpath(join(fsrc, algorithm))
                f77_extension = re.compile('.*\.f$')
                f_sources = [join(sub_fsrc, fortran) for fortran in filter(f77_extension.match, listdir(sub_fsrc))]
                command.extend(f_sources)
                system_stdout(command)

                # compilation of the classical version
                signature = relpath(join(gateways_classical, '{}-interface.pyf'.format(algorithm)))
                module = relpath(join(gateways_classical, '{}.f90'.format(algorithm)))
                command = compiler + [signature, module]

                sub_fsrc = relpath(join(fsrc_classical, algorithm))
                f_sources = [join(sub_fsrc, fortran) for fortran in filter(f77_extension.match, listdir(sub_fsrc))]
                command.extend(f_sources)
                system_stdout(command)
                print('Done.')

            # clean up the working directory
            remove_temporary_compilation(c_pwd)

            # move the binary-file to the 'interfaces' folder
            so_extension = re.compile('.*\.(so|pyd)$')
            for move_file in filter(so_extension.match, listdir(c_pwd)):
                pdfo_move_file = join(interfaces, move_file)
                if isfile(pdfo_move_file):
                    remove(pdfo_move_file)
                rename(join(c_pwd, move_file), join(interfaces, move_file))

            print('\nPackage compiled successfully!', end='\n\n')

            # add pdfo to the PYTHONPATH environment variable in the shell initialization scripts
            pdfo_to_path = False
            rc_files = []
            rc_files_failed = []
            print('Adding PDFO to PYTHONPATH will enable you to use the package in other shell sessions.')
            msg = 'Add PDFO to the PYTHONPATH in the shell initialization script ([Y]/n)? '
            nb_test = 0
            correct_input = False
            pdfo_to_path_choice = False
            while nb_test < 3 and not correct_input:
                nb_test += 1
                if python_version >= 3:
                    get_input = input(msg)
                else:
                    get_input = raw_input(msg)

                if get_input.lower() in ['', 'y', 'yes']:
                    correct_input = True
                    pdfo_to_path_choice = True
                elif get_input.lower() in ['n', 'no']:
                    correct_input = True
                elif nb_test < 3:
                    print('Sorry, we did not understand your input. Try again.')
            print()
            if nb_test == 3 and not correct_input:
                print('Sorry, we did not understand your inputs.', end='\n\n')
            elif system_os == 'unix' and pdfo_to_path_choice:
                # the user chose to add PDFO to the PYTHONPATH in his shell initialization script
                rc_list = ['~/.bashrc', '~/.zshrc', '~/.fishrc', '~/.kshrc', '~/.tcshrc']
                for rc in rc_list:
                    if isfile(expanduser(rc)):
                        rc_files.append(expanduser(rc))

                if len(rc_files) == 0:
                    # no common shell initialization script has been found
                    msg = 'Path to your shell initialization script? '
                    nb_test = 0
                    correct_input = False
                    while nb_test < 3 and not correct_input:
                        nb_test += 1
                        if python_version >= 3:
                            get_input = input(msg)
                        else:
                            get_input = raw_input(msg)

                        if isfile(expanduser(get_input)):
                            correct_input = True
                            rc_files.append(expanduser(get_input))
                        elif nb_test < 3:
                            print('Sorry, try again.')
                    print()

                    if nb_test == 3 and not correct_input:
                        print('3 incorrect answer attempts.', end='\n\n')

                # add PDFO to PYTHONPATH in every shell initialization script found
                for expand_rc in rc_files:
                    try:
                        export_stmt = 'export PYTHONPATH={}:$PYTHONPATH\n'.format(interfaces_path)
                        already_exported = False
                        fd = open(expand_rc, 'r')
                        for line in fd:
                            already_exported = already_exported or (export_stmt in line)
                        fd.close()
                        if not already_exported:
                            fd = open(expand_rc, 'a')
                            fd.write('\n# Added by PDFO ({})\n'.format(datetime.now().strftime('%d/%m/%Y %H:%M')))
                            fd.write(export_stmt)
                            fd.close()
                        pdfo_to_path = True
                    except:
                        rc_files_failed.append(expand_rc)
                        print('Failed to add PDFO to the PYTHONPATH in {}.'.format(expand_rc), end='\n\n')
            elif system_os == 'win' and pdfo_to_path_choice:
                # the user chose to add PDFO to the PYTHONPATH
                if 'PYTHONPATH' in dict(os.environ).keys():
                    python_path = os.environ['PYTHONPATH']
                else:
                    python_path = ''

                if interfaces_path not in python_path:
                    if 'PYTHONPATH' in dict(os.environ).keys():
                        cmd = ['setx', 'PYTHONPATH', '{};{}'.format(interfaces_path, python_path)]
                    else:
                        cmd = ['setx', 'PYTHONPATH', interfaces_path]
                    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
                    _, error = p.communicate()

                    if p.returncode == 0:
                        pdfo_to_path = True
                    else:
                        print('Failed to add PDFO to the PYTHONPATH.')
                else:
                    pdfo_to_path = True

            # print success message
            try:
                columns_term, _ = os.get_terminal_size(0)
            except:
                columns_term = 72
            star_width_env = min(columns_term, max(49, len(interfaces_path)))
            star_width_fail = min(columns_term, max(37 + len(', '.join(rc_files_failed)), 23 + len(interfaces_path)))
            star_width_source = min(columns_term, max(66, 7 + len(' '.join(rc_files))))
            star_width_set = min(columns_term, max(66, 28 + len(interfaces)))

            if not pdfo_to_path:
                print('-' * star_width_env)
                print('Please add the following path to your PYTHONPATH:')
                print(interfaces_path)
                print('-' * star_width_env, end='\n\n')
            elif len(rc_files_failed) != 0:
                print('-' * star_width_fail)
                print('Please add the following command to {}:'.format(', '.join(rc_files_failed)))
                print('PYTHONPATH={}:$PYTHONPATH'.format(interfaces_path))
                print('-' * star_width_fail, end='\n\n')
            elif system_os != 'win' or pdfo_to_path:
                print('PDFO added successfully to PYTHONPATH.', end='\n\n')
            if pdfo_to_path and system_os == 'unix':
                print('-' * star_width_source)
                print('To use PDFO in the current session, execute the following command:')
                print('source {}'.format(' '.join(rc_files)))
                print('-' * star_width_source, end='\n\n')
            elif pdfo_to_path and system_os == 'win':
                print('-' * star_width_set)
                print('To use PDFO in the current session, execute the following command:')
                print('set PYTHONPATH={};%PYTHONPATH%'.format(interfaces_path))
                print('-' * star_width_set, end='\n\n')
            print('You may afterwards try the following command in Python for the usage of the package:')
            print('>>> import pdfo')
            print('>>> help(pdfo.pdfo)', end='\n\n')
            print("You may also test the package by executing the script 'testpdfo.py' via")
            print('{} {}'.format(python_cmd, join(pdfo_pwd, 'examples', 'testpdfo.py')), end='\n\n')
        except:
            remove_temporary_compilation(c_pwd, include_bin=True)
            print('\nNumPy and F2PY are installed but the compilation failed.')
            if system_os == 'win':
                print('It is possibly because gcc-fortran is not properly installed.')
            else:
                print('It is possibly because gfortran is not properly installed.')
            if system_known:
                print('Please see `README_py_{}.txt`.'.format(system_os))
            else:
                print('Please see the README file corresponding to your system.')
