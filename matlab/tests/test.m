function test(varargin)
%This function tests a modernized version of Powell's solver against Powell's version, verifying
% whether they produce the same results on CUTEst problems.
%
% Usage: test(solver, dimrange) , where `solver` is the name of the solver to test, while `dimrange`
% is the vector [mindim, maxdim], or 'small', or 'big'"
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: July 2020
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wrong_input = false;
if nargin == 1
    solver = varargin{1};
    mindim = 1;
    maxdim = 100;
elseif nargin == 2
    if ischstr(varargin{1}) && isnumvec(varargin{2})
        solver = varargin{1};
        mindim = min(varargin{2});
        maxdim = max(varargin{2});
    elseif ischstr(varargin{2}) && isnumvec(varargin{1})
        solver = varargin{2};
        mindim = min(varargin{1});
        maxdim = max(varargin{1});
    elseif ischstr(varargin{1}) && ischstr(varargin{2})
        if strcmpi(varargin{1}, 'small')
            solver = varargin{2};
            mindim = 1;
            maxdim = 50;
        elseif strcmpi(varargin{1}, 'big')
            solver = varargin{2};
            mindim = 51;
            maxdim = 100;
        elseif strcmpi(varargin{2}, 'small')
            solver = varargin{1};
            mindim = 1;
            maxdim = 50;
        elseif strcmpi(varargin{2}, 'big')
            solver = varargin{1};
            mindim = 51;
            maxdim = 100;
        else
            wrong_input = true;
        end
    else
        wrong_input = true;
    end
else
    wrong_input = true;
end

if wrong_input
    error(sprintf("\nUsage:\n\n\ttest(solver, dimrange) ,\n\nwhere `solver` is the name of the solver to test, while `dimrange` is the vector [mindim, maxdim], or 'small', or 'big'.\n"));
end


% Make CUTEst available. The following three lines should be configured to fit the installation of
% CUTEst on the current machine.
setenv('CUTEST', '~/local/cutesif/cutest');
setenv('MASTSIF', '~/local/cutesif/sif');
addpath('~/local/cutesif/mtools/msrc');

current_dir = cd();
neupdfo_dir = fileparts(fileparts(current_dir));
opdfo_dir = fullfile(neupdfo_dir, 'OPDFO');

mexopt = struct();
mexopt.debug = true;

solvern = [solver, 'n'];
cd(neupdfo_dir);
setup(solvern, mexopt);
cd(opdfo_dir);
setup(solver, mexopt)
cd(current_dir);

solvers = {[solver, 'n'], solver};
options = struct();
options.mindim = mindim;
options.maxdim = maxdim;
options.nr = 20;
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

assert(isequiv(solvers, options));

fprintf('\n\nThe test on %s is successful!\n\n', solver);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%
function is_numvec = isnumvec(x)
is_numvec = isnumeric(x) && isvector(x);

function is_chstr = ischstr(x)
is_chstr = isa(x, 'char') || isa(x, 'string');
