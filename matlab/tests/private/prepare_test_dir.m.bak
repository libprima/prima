function test_dir = prepare_test_dir(solver, test_name, test_options)
%PREPARE_TEST_DIR prepares a copy of the package in `test_dir` for the test named `test_name`.

if verLessThan('matlab', '9.10') && ispc
    callstack = dbstack;
    funname = callstack(1).name; % Name of the current function
    warning('%s may FAIL because ''copyfile'' of MATLAB R2020b or earlier may raise an error when handling unix symbolic links under Windows.', funname);
end

mfilepath = fileparts(mfilename('fullpath'));  % Directory where this .m file resides.

% In the sequel, 's_' indicates source, 'd_' indicates destination.

% `root_dir`: root directory of the project
s_root_dir = fileparts(fileparts(fileparts(mfilepath)));

% Define the name for the test directory. We do not add a rand stamp to `test_dir_name`, so that we
% can reuse it easily. It is needed if we want to repeat a test without recompiling the solvers.
root_dir_name = strrep(s_root_dir, [fileparts(s_root_dir), filesep],'');
test_dir_name = [solver, '_', test_name, '_', root_dir_name];
test_dir = fullfile(tempdir, test_dir_name);  % Full path to the test directory

% Check that `test_dir` is not a substring of `mfilepath`, so that `mfilepath` is not subdirectory of
% `test_dir`. Also check that `s_root_dir` is not a subdirectory of `test_dir`. Without doing this
% we may mistakenly remove important data.
assert(~contains(mfilepath, test_dir));
assert(~contains(s_root_dir, test_dir));

d_root_dir = test_dir;
% If test_options.compile = true, then remove `d_root_dir` to make sure that we are testing the latest
% code in `s_root_dir`
if test_options.compile
    if exist(d_root_dir, 'dir')
        rmdir(d_root_dir, 's');
    end
end
if ~exist(d_root_dir, 'dir')
    mkdir(d_root_dir);
end

% `matlab_dir`: the `matlab` directory
s_matlab_dir = fullfile(s_root_dir, 'matlab');
d_matlab_dir = fullfile(d_root_dir, 'matlab');
if ~exist(d_matlab_dir, 'dir')
    mkdir(d_matlab_dir);
end

% The following lines get the path for SOLVER_archiva
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% `dev_arch` a subdirectory of fullfile(s_root_dir, '.development', 'archiva'). It contains the "archiva"
% version of solvers used as a benchmark for the development of the current version of the solvers.
archiva_dir_name = 'dev_arch';
dev_arch_dir = fullfile(s_root_dir, '.development', 'archiva', archiva_dir_name);
if isunix && ~ismac
    [~, dev_arch_dir] = system(['realpath ', dev_arch_dir]);
end
fprintf('\nThe archiva version: %s\n', dev_arch_dir);
% Define `archiva_dir` as the "norma" directory under the `dev_arch` directory. Indeed, the solvers in
% fullfile(s_root_dir, '.development', 'archiva', archiva_dir_name) and
% fullfile(s_root_dir, '.development', 'archiva', archiva_dir_name, 'norma')
% are equivalent. We use the latter because the name of the solver there is SOLVER_norma, which is
% convenient for the test.
s_archiva_dir = fullfile(s_root_dir, '.development', 'archiva', archiva_dir_name, 'norma');
d_archiva_dir = fullfile(d_root_dir, '.development', 'archiva', archiva_dir_name, 'norma');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%!!!----------------------------------------------------------------------------------------!!!%%%
% Remove the existing compiled MEX files in the source directories. This is IMPORTANT!
% Without, MATLAB will crash due to a bug of MATLAB under Linux. This was observed on
% 2022-02-16 and took two days to fix. See https://github.com/zequipe/test_matlab/blob/master/crash/.
% It must be done BEFORE copying the files.
setup_tools = fullfile(s_root_dir, 'matlab', 'setup_tools');
addpath(setup_tools);  % We use `clean_mex` from `setup_tools` to clean up the compiled MEX files.
s_mexdir = fullfile(s_matlab_dir, 'interfaces', 'private');
s_norma_mexdir = fullfile(s_root_dir, '.development', 'norma', 'matlab', 'interfaces', 'private');
s_archiva_mexdir = fullfile(s_archiva_dir, 'matlab', 'interfaces', 'private');
mexdir_list = {s_mexdir, s_norma_mexdir, s_archiva_mexdir};
cellfun(@clean_mex, mexdir_list);
rmpath(setup_tools);  % Remove `setup_tools` from path since it has finishes its job.
%%%!!!----------------------------------------------------------------------------------------!!!%%%


% Remove the interform Fortran source code under `s_root_dir`. This is needed if we would like to
% repeat a debugging test without recompiling the solvers; without doing this, the interform Fortran
% code in the destination directory would be overwritten, which makes debugging impossible.
s_fortran_interform = fullfile(s_root_dir, 'fortran', '.interform');
if exist(s_fortran_interform, 'dir')
    rmdir(s_fortran_interform, 's');
end
s_mex_interform = fullfile(s_root_dir, 'matlab', 'mex_gateways', '.interform');
if exist(s_mex_interform, 'dir')
    rmdir(s_mex_interform, 's');
end

% `root_sub_list: directories/files to be copied under `root_dir`
root_sub_list = {'fortran', '.development/norma', 'setup.m'};
for il = 1 : length(root_sub_list)
    copyfile(fullfile(s_root_dir, root_sub_list{il}), fullfile(d_root_dir, root_sub_list{il}));
end
%!------------------------------------------------------------------------------------------------!%
% Remove the .development/norma/matlab/tests directories under d_root_dir. This is IMPORTANT! Without
% this, the test will mistakenly call scripts from there, because norma/matlab/tests may be added to
% the path by the corresponding setup.m. This problem occurred on 2022-02-18 and took two days to debug.
d_norma_matlab_test = fullfile(d_root_dir, '.development', 'norma', 'matlab', 'tests');
if exist(d_norma_matlab_test, 'dir')
    rmdir(d_norma_matlab_test, 's');
end
mkdir(d_norma_matlab_test);  % Necessary, because `setup.m` may try adding this directory to path
%!------------------------------------------------------------------------------------------------!%

% `matlab_sub_list`: directories/files to be copied under `matlab_dir`
matlab_sub_list = {'setup_tools', 'mex_gateways', 'interfaces'};
for il = 1 : length(matlab_sub_list)
    copyfile(fullfile(s_matlab_dir, matlab_sub_list{il}), fullfile(d_matlab_dir, matlab_sub_list{il}));
end
%!------------------------------------------------------------------------------------------------!%
% Remove the tests directories under d_matlab_dir.
prima_matlab_test = fullfile(d_matlab_dir, 'tests');
if exist(prima_matlab_test, 'dir')
    rmdir(prima_matlab_test, 's');
end
mkdir(prima_matlab_test);  % Necessary, because `setup.m` may try adding this directory to path
%!------------------------------------------------------------------------------------------------!%
% `testprima` is needed by `pdv`.
copyfile(fullfile(s_matlab_dir, 'tests', 'testprima.m'), fullfile(prima_matlab_test, 'testprima.m'));

% `archiva_sub_list`: directories/files to be copied under `archiva_dir`
archiva_sub_list = {'matlab', 'fortran', 'setup.m'};
for il = 1 : length(archiva_sub_list)
    copyfile(fullfile(s_archiva_dir, archiva_sub_list{il}), fullfile(d_archiva_dir, archiva_sub_list{il}));
end
%!------------------------------------------------------------------------------------------------!%
% Remove the matlab/tests directories under d_archiva_dir.
d_archiva_matlab_test = fullfile(d_archiva_dir, 'matlab', 'tests');
if exist(d_archiva_matlab_test, 'dir')
    rmdir(d_archiva_matlab_test, 's');
end
mkdir(d_archiva_matlab_test);  % Necessary, because `setup.m` may try adding this directory to path
%!------------------------------------------------------------------------------------------------!%
