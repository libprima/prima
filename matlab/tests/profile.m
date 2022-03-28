function output = profile(varargin)
%This function profiles the modernized version of Powell's solver against Powell's version.
%
% Usage:
%
%   profile(solver, dimrange, nocompile_flag, sequential_flag, reverse_flag, problem_type, options)
%   profile(solver, dimrange, reload_flag, reverse_flag, problem_type, options)
%
% where
% - `solver` is the name of the solver to test
% - `dimrange` is the vector [mindim, maxdim], or "small", or "big", or "large"
% - `nocompile_flag` is either 'nocompile' or 'ncp', indicating not to compile the solves
% - `reload_flag` is either 'reload' or 'load', indicating to load the data directly from the .mat
%   file corresponding to `solver` and `dimrange`
% - `sequential_flag` (optional) is either 'sequential' or 'seq', which means to test the problems sequentially
% - `reverse_flag` (optional) is either 'reverse' or 'rev', which means to test the solvers in the reverse order
% - `problem_type` can be any of {'u', 'b', 'l', 'n', 'ub', 'ubl', 'ubln', 'bl', 'bln', 'ln'},
%   indicating the problem type to test
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: October 2021
%
% Last Modified: Monday, February 12, 2022 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldpath = path();  % Record the current path.
restoredefaultpath;  % Restore the "right out of the box" path of MATLAB

olddir = pwd();  % Record the current directory.

% Parse the inputs.
[solver, options] = parse_input(varargin);


% Set up the directory to save the testing data, i.e., `data_dir`.
if isfield(options, 'data_dir')
    data_dir = options.data_dir;
else
    mfiledir = fileparts(mfilename('fullpath')); % Directory where this .m file resides.
    data_dir = fullfile(mfiledir, 'testdata');
end
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

% Prepare the test directory, i.e., `test_dir`.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
test_dir = prepare_test_dir(solver, funname, options);

exception = [];

try

    % Specify where to store the test data.
    options.data_dir = data_dir;

    % Test feature and time
    test_feature = '';
    if isfield(options, 'compiler_options') && (isa(options.compiler_options, 'char') || ...
          isa(options.compiler_options, 'string')) && ~isempty(options.compiler_options)
        test_feature = [test_feature, '.', regexprep(options.compiler_options, '\s*','_')];
    end
    if isfield(options, 'perm') && islogical(options.perm) && isscalar(options.perm) && options.perm
        test_feature = [test_feature, '.', 'perm'];
    end
    if isfield(options, 'randomizex0') && isnumeric(options.randomizex0) && isscalar(options.randomizex0) && abs(options.randomizex0) > 0
        test_feature = [test_feature, '.', 'randomizex0', sprintf('%g', options.randomizex0)];
    end
    if isfield(options, 'eval_options') && isstruct(options.eval_options) && ~isempty(fieldnames(options.eval_options))
        test_feature = [test_feature, '.', strjoin(fieldnames(options.eval_options), '_')];
        if isfield(options.eval_options, 'dnoise')
            if isnumeric(options.eval_options.dnoise) && isscalar(options.eval_options)
                dnoise_level = abs(options.eval_options.dnoise);
            elseif isstruct(options.eval_options.dnoise) && isfield(options.eval_options.dnoise, 'level')
                dnoise_level = abs(options.eval_options.dnoise.level);
            else
                dnoise_level = 0;
            end
            if dnoise_level > 0
                test_feature = regexprep(test_feature, 'dnoise', ['dnoise', sprintf('%g', dnoise_level)]);
            end
        end
        if isfield(options.eval_options, 'noise')
            if isnumeric(options.eval_options.noise) && isscalar(options.eval_options)
                noise_level = abs(options.eval_options.noise);
            elseif isstruct(options.eval_options.noise) && isfield(options.eval_options.noise, 'level')
                noise_level = abs(options.eval_options.noise.level);
            else
                noise_level = 0;
            end
            if noise_level > 0
                test_feature = regexprep(test_feature, 'noise', ['noise', sprintf('%g', noise_level)]);
            end
        end
        if isfield(options.eval_options, 'signif')
            test_feature = regexprep(test_feature, 'signif', ['signif', sprintf('%d', options.eval_options.signif)]);
        end
    end
    test_feature = regexprep(test_feature, '^\.', '');
    if isempty(test_feature)
        test_feature = 'plain';
    end
    options.test_feature = test_feature;

    if ~isfield(options, 'time')
        options.time = datestr(datetime(), 'yymmdd_HHMM');
    end

    % Make the solvers available. Note that the solvers are under `test_dir`.
    get_solvers(solver, test_dir, options);

    % Tell MATLAB where to find CUTEST.
    locate_cutest();

    % Go to the test directory. This is not really necessary. It will not affect the test, but any
    % output (e.g., NEWUOA_output.txt, fort.6) will be dumped to `test_dir`.
    cd(test_dir);

    % Define the solvers to test.
    if isfield(options, 'reverse') && options.reverse
        solvers = {[solver, 'n_classical'], [solver, 'n']};  % Reverse order: first run 'SOLVERn_classical', and then run 'SOLVERn'
    else
        solvers = {[solver, 'n'], [solver, 'n_classical']};  % % Default order: first run 'SOLVERn', and then run 'SOLVERn_classical'.
    end

    % The following line can be used for testing the single-precision version. If such a test is
    % intended, remember to set mexopt.single = true in `get_solvers`.
    %solvers = {[solver, 'n'], [solver, 'n_single']};

    % Show current path information.
    showpath(solvers);

    % Profile the solvers.
    tic;
    output = struct(options.type, perfdata(solvers, options));
    problem_type=options.type;
    if length(problem_type) > 1
        options.reload = true;
        if strcmpi(solver, 'cobyla') && contains(problem_type, 'n')
            options.type = 'n';
            output.n = perfdata(solvers, options);
        end
        if ismember(solver, {'cobyla', 'lincoa'}) && contains(problem_type, 'l')
            options.type = 'l';
            output.l = perfdata(solvers, options);
        end
        if ismember(solver, {'cobyla', 'lincoa', 'bobyqa'}) && contains(problem_type, 'b')
            options.type = 'b';
            output.b = perfdata(solvers, options);
        end
        if ismember(solver, {'cobyla', 'lincoa', 'bobyqa'}) && contains(problem_type, 'u')
            options.type = 'u';
            output.u = perfdata(solvers, options);
        end
    end
    toc;

    % Show current path information again at the end of test.
    showpath(solvers);

catch exception

    % Do nothing for the moment.

end

setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.
fprintf('\nCurrently in %s\n\n', pwd());

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end
