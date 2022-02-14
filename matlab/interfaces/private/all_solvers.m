function solver_list = all_solvers(solver_type)
%ALL_SOLVERS returns a cell array containing the names of all the solvers of `solver_type` contained
% in this package.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

if nargin > 1 || (nargin == 1 && ~ischarstr(solver_type))
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid input received', funname);
end

all_solvers_available = {'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan'};
solvers_without_constraints = {'uobyqan', 'newuoan'};
solvers_with_constraints = setdiff(all_solvers_available, solvers_without_constraints);
unconstrained_solvers = solvers_without_constraints;
bound_constrained_solvers = {'bobyqan', 'lincoan', 'cobylan'};
linearly_constrained_solvers = {'lincoan', 'cobylan'};
nonlinearly_constrained_solvers = {'cobylan'};
internal_solvers = {'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan'};

if nargin < 1
    solver_list = all_solvers_available;
else
    switch lower(solver_type)
    case 'all'
        solver_list = all_solvers_available;
    case 'without_constraints'
        solver_list = solvers_without_constraints;
    case 'with_constraints'
        solver_list = solvers_with_constraints;
    case 'unconstrained_solvers'
        solver_list = unconstrained_solvers;
    case 'bound_constrained_solvers'
        solver_list = bound_constrained_solvers;
    case 'linearly_constrained_solvers'
        solver_list = linearly_constrained_solvers;
    case 'nonlinearly_constrained_solvers'
        solver_list = nonlinearly_constrained_solvers;
    case 'internal'
        solver_list = internal_solvers;
    otherwise
         % Private/unexpected error
         error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: unknown solver type ''%s'' received', funname, solver_type);
    end
end

return
