function [solver, options] = parse_input(argin)
%This function parses the input to a testing function, returning the name of the solver to test and
% the testing options. The testing function can have the signature
%
%   test(solver, dimrange, nocompile_flag)
%   test(solver, dimrange, reload_flag)
%
% where
% - `solver` is the name of solver to test
% - `dimrange` (optional) is the vector [mindim, maxdim], or "small", or "big", or "large", or "all"
% - `nocompile` (optional) is either 'nocompile' or 'ncp', which means not to compile the solvers
% - `reload_flag` is either 'reload' or 'load', indicating not load the data directly from the .mat
% file corresponding to `solver` and `dimrange`
%
% If the testing function is `verify`, then the following signatures are also supported:
%
%   verify(solver, problem, nocompile_flag)
%   verify(solver, {problem, ir}, nocompile_flag)
%
% where
% - `problem` is a problem to test
% - `ir` is the index of the random run in `verify`
%
% If the testing function is `profile`, then the following signature is also supported
%
%   profile(solver, dimrange, reload_flag)
%
% where
% - `reload_flag` is either 'reload' or 'load', indicating not load the data directly from the .mat
% file corresponding to `solver` and `dimrange`
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: July 2020
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

callstack = dbstack;
invoker = callstack(2).name;  % The function that calls this function.
known_solvers = {'cobyla','uobyqa','newuoa','bobyqa','lincoa'};
known_sizes = {'small', 'big', 'large', 'all'};
nocomplie_flags = {'nocompile', 'ncp'};
reload_flags = {'reload', 'load'};

% Default values.
solver = '';
prob = '';
ir = NaN;
mindim = 1;
maxdim = 50;
compile = true;
reload = false;

wrong_input = (length(argin) < 1 || length(argin) > 3);


if length(argin) == 3
    if ischstr(argin{1}) && ismember(argin{1}, nocomplie_flags)
        compile = false;
        argin = argin([2, 3]);
    elseif ischstr(argin{2}) && ismember(argin{2}, nocomplie_flags)
        compile = false;
        argin = argin([1, 3]);
    elseif ischstr(argin{3}) && ismember(argin{3}, nocomplie_flags)
        compile = false;
        argin = argin([1, 2]);
    else
        wrong_input = true;
    end
end

if length(argin) == 3 && strcmp(invoker, 'profile')
    if ischstr(argin{1}) && ismember(argin{1}, reload_flags)
        reload = true;
        compile = false;
        argin = argin([2, 3]);
    elseif ischstr(argin{2}) && ismember(argin{2}, reload_flags)
        reload = true;
        compile = false;
        argin = argin([1, 3]);
    elseif ischstr(argin{3}) && ismember(argin{3}, reload_flags)
        reload = true;
        compile = false;
        argin = argin([1, 2]);
    else
        wrong_input = true;
    end
end

if length(argin) == 2
    if ischstr(argin{1}) && isnumvec(argin{2}) && length(argin{2}) == 2
        solver = argin{1};
        mindim = min(argin{2});
        maxdim = max(argin{2});
    elseif ischstr(argin{2}) && isnumvec(argin{1}) && length(argin{1}) == 2
        solver = argin{2};
        mindim = min(argin{1});
        maxdim = max(argin{1});
    elseif ischstr(argin{1}) && ischstr(argin{2})
        argin = lower(argin);
        if length(intersect(argin, known_solvers)) == 1
            solver = intersect(argin, known_solvers);
            solver = solver{1};
            sdiff = setdiff(argin, known_solvers);
            if length(intersect(sdiff, known_sizes)) == 1
                [mindim, maxdim, wrong_input] = parse_dim(sdiff{1});
            elseif strcmp(invoker, 'verify')
                prob = lower(sdiff{1});
            else
                wrong_input = true;
            end
        else
            wrong_input = true;
        end
    elseif ischstr(argin{1}) && iscell(argin{2}) && strcmp(invoker, 'verify')
        solver = argin{1};
        [prob, ir, wrong_input] = parse_prob_ir(argin{2});
    elseif ischstr(argin{2}) && iscell(argin{1}) && strcmp(invoker, 'verify')
        solver = argin{2};
        [prob, ir, wrong_input] = parse_prob_ir(argin{1});
    else
        wrong_input = true;
    end
end

if length(argin) == 1
    if ischstr(argin{1})
        solver = argin{1};
    else
        wrong_input = true;
    end
end

solver = lower(solver);

wrong_input = wrong_input || ~(ismember(solver, known_solvers) && (mindim <= maxdim) && mindim >= 1 && isint(mindim) && isint(maxdim));

if wrong_input
    if (strcmp(invoker, 'verify'))
        errmsg = sprintf('\nUsage:\n\n\t%s(solver, dimrange, nocomplie_flag), or %s(solver, {problem, ir}, nocomplie_flag).\n', invoker, invoker);
    elseif (strcmp(invoker, 'profile'))
        errmsg = sprintf('\nUsage:\n\n\t%s(solver, dimrange, nocompile_flag), or %s(solver, dimrange, reload_flag).\n', invoker, invoker);
    else
        errmsg = sprintf('\nUsage:\n\n\t%s(solver, dimrange, nocompile_flag).\n', invoker);
    end
    error(errmsg)
end

% Define the testing options.
options = struct();
options.compile = compile;
options.reload = reload;
if isempty(prob)
    options.mindim = mindim;
    options.maxdim = maxdim;
    options.nr = 20;
    % Define the problem type(s) to test.
    switch solver
    case {'uobyqa', 'newuoa'}
        options.type = 'u';
    case 'bobyqa'
        options.type = 'bu';
    case 'lincoa'
        options.type = 'lbu';
    otherwise
        options.type = 'nlbu';
    end
else
    options.list = {prob};
    if (~isnan(ir))
        options.ir = ir;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mindim, maxdim, wrong_input] = parse_dim(dimrange)

wrong_input = false;
switch lower(dimrange)
case 'small'
    mindim = 1;
    maxdim = 50;
case 'big'
    mindim = 51;
    maxdim = 100;
case 'large'
    mindim = 101;
    maxdim = 200;
case 'all'
    mindim = 1;
    maxdim = 200;
otherwise
    wrong_input = true;
end


function [prob, ir, wrong_input] = parse_prob_ir(prob_ir)

wrong_input = false;
if ischstr(prob_ir{1}) && isnumvec(prob_ir{2})
    prob = prob_ir{1};
    ir = prob_ir{2};
elseif ischstr(prob_ir{2}) && isnumvec(prob_ir{1})
    prob = prob_ir{2};
    ir = prob_ir{1};
else
    wrong_input = true;
end

wrong_input = wrong_input || ~(isint(ir) && ir >= 0);
