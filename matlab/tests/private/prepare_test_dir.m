function test_dir = prepare_test_dir(solver, test_name, test_options)
%PREPARE_TEST_DIR prepares a copy of the package in `test_dir` for the test named `test_name`.

if verLessThan('matlab', '9.10') && ispc
    callstack = dbstack;
    funname = callstack(1).name; % Name of the current function
    warning('%s may FAIL because ''copyfile'' of MATLAB R2020b or earlier may raise an error when handling unix symbolic links under Windows.', funname);
end

fprintf('\nPreparing the testing directory for test ''%s'' ...\n', test_name);
fprintf('\nThe solver: %s\n', solver);

competitor = '';
if isfield(test_options, 'competitor')
    competitor = test_options.competitor;
    fprintf('\nThe competitor: %s\n', test_options.competitor);
end

compile = ~isfield(test_options, 'compile') || test_options.compile;
fprintf('\nCompile: %d\n', compile);

% mfilepath: directory where this .m file resides.
mfilepath = fileparts(mfilename('fullpath'));

% root_dir: root directory of the project
root_dir = fileparts(fileparts(fileparts(mfilepath)));
fprintf('\nThe root directory: %s\n', root_dir);

% test_dir: the test directory. We do not add a randomized stamp to `test_dir_name`, so that we
% can reuse it easily. It is needed if we want to repeat a test without recompiling the solvers.
[~, root_dir_name] = fileparts(root_dir);
test_dir_name = [solver, '_', test_name, '_', root_dir_name];
test_dir = fullfile(tempdir, test_dir_name);  % Full path to the test directory
fprintf('\nThe test directory: %s\n', test_dir);

% Check that `test_dir` is not a substring of `mfilepath`, so that `mfilepath` is not subdirectory of
% `test_dir`. Also check that `root_dir` is not a subdirectory of `test_dir`. Without doing this
% we may mistakenly remove important data.
assert(~contains(mfilepath, test_dir));
assert(~contains(root_dir, test_dir));

% If compile = true, then remove `test_dir` to make sure that it will contain the latest code in `root_dir`
if compile
    if exist(test_dir, 'dir')
        rmdir(test_dir, 's');
    end
end
if ~exist(test_dir, 'dir')
    mkdir(test_dir);
end

% norma_dir: subdirectory of that contains the "norma" version of solvers used as the reference
% implementation for the development of the current version of the solvers.
norma_dir = fullfile(root_dir, '.development', 'norma');

% dev_arch: a subdirectory of fullfile(root_dir, '.development', 'archiva'). It contains the "archiva"
% version of solvers used as a benchmark for the development of the current version of the solvers.
archiva_dir = fullfile(root_dir, '.development', 'archiva', 'dev_arch', 'norma');

% s_dir_list: a list of directories to be copied to `test_dir`, "s" for "source".
s_dir_list = {root_dir};
% d_dir_list: a list of directories to receive the directories copied from `s_dir_list`, "d" for
% "destination".
d_dir_list = {fullfile(test_dir, root_dir_name)};
switch competitor
case 'norma'
    s_dir_list = [s_dir_list, {norma_dir}];
    fprintf('\nThe norma directory: %s\n', norma_dir);
    d_dir_list = [d_dir_list, {fullfile(test_dir, 'norma')}];
case 'archiva'
    s_dir_list = [s_dir_list, {archiva_dir}];
    if isunix && ~ismac
        [~, archiva_dir] = system(['realpath ', archiva_dir]);  % archiva_dir may be a symlink.
    end
    fprintf('\nThe archiva directory: %s\n', archiva_dir);
    d_dir_list = [d_dir_list, {fullfile(test_dir, 'archiva')}];
otherwise
    % Do nothing
end

% fdlist: a list of files/directories to be copied to `test_dir` from each directory in `s_dir_list`.
fdlist = {'matlab', 'setup.m'};
if compile
    fdlist = [fdlist, {'fortran'}];
end

% We use `clean_mex` from `setup_tools` to clean up the compiled MEX files.
setup_tools = fullfile(root_dir, 'matlab', 'setup_tools');
addpath(setup_tools);

for idir = 1 : length(s_dir_list)
    s_dir = s_dir_list{idir};
    d_dir = d_dir_list{idir};
    if ~exist(d_dir, 'dir')
        mkdir(d_dir);
    end

    %%%!!!------------------------------------------------------------------------------------!!!%%%
    % Remove the existing compiled MEX files in the source directories. This is IMPORTANT!
    % Without doing so, MATLAB will crash due to a bug of MATLAB under Linux. It must be done
    % before copying the files. This was observed on 2022-02-16 and took two days to fix.
    % See https://github.com/zequipe/test_matlab/blob/master/crash/.
    s_mexdir = fullfile(s_dir,  'matlab', 'interfaces', 'private');
    clean_mex(s_mexdir);
    %%%!!!------------------------------------------------------------------------------------!!!%%%

    % Copy the files in `s_dir` to `d_dir`
    for ifd = 1 : length(fdlist)
        fd = fdlist{ifd};
        copyfile(fullfile(s_dir, fd), fullfile(d_dir, fd));
    end

    % Empty fullfile(d_dir, 'matlab', 'tests') to make sure that we will call the testing scripts in
    % `root_dir`, not those in `test_dir`.
    d_matlab_tests = fullfile(d_dir, 'matlab', 'tests');
    if exist(d_matlab_tests, 'dir')
        rmdir(d_matlab_tests, 's');
    end

    % It is necessary to create an empty "tests" directory, as `setup.m` may try adding it to path.
    mkdir(d_matlab_tests);

    % `testprima` is needed by `pdv`.
    if strcmp(s_dir, root_dir) && strcmp(test_name, 'pdv')
        copyfile(fullfile(s_dir, 'matlab', 'tests', 'testprima.m'), fullfile(d_matlab_tests, 'testprima.m'));
    end
end

% Remove `setup_tools` from path since it has finishes its job.
rmpath(setup_tools);
