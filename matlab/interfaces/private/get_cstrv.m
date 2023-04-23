function cstrv = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq)
%GET_CSTRV calculates the (absolute) constraint violation of x. Note that `nlcineq` and `nlceq` are
% values of the nonlinear inequality and equality constraints, the inequality constraint being
% nlcineq <= 0. `nlcineq` and `nlceq` are optional, while all the other inputs are obligatory.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if nargin < 7
    % Private/unexpected error
    error(sprintf('%s:InvalidNargin', funname), ...
    '%s: UNEXPECTED ERROR: at least 8 inputs expected.', funname);
end
if nargin < 8
    nlcineq = [];
end
if nargin < 9
    nlceq = [];
end

rlb = [];
if ~isempty(lb)
    rlb = lb - x;
    rlb(lb <= x) = 0;  % Prevent NaN in case lb = x = +/-Inf; OR: rlb = rlb(~(lb <= x))
end

rub = [];
if ~isempty(ub)
    rub = x - ub;
    rub(x <= ub) = 0;  % Prevent NaN in case ub = x = +/-Inf; OR: rub = rub(~(x <= ub))
end

rineq = [];
if ~isempty(Aineq)
    rineq = Aineq*x - bineq;
    % Do not write `rineq(Aineq*x <= bineq) = 0`, as we want to get NaN if Aineq*x = bineq = +/-Inf.
    % What if Aineq*x = bineq = +Inf? Shouldn't such a constraint has a zero violation? Assume that
    % Aineq*x = +Inf implies either Aineq or x contains infinite values, which is true
    % mathematically but false numerically. Then we take the view that it is better to return NaN
    % than 0, because either the constraint or x is problematic.
    % This function is called only by preprima and postprima. With `rineq(Aineq*x <= bineq) = 0`,
    % the same constraint may not return the same constraint violation for preprima and postprima
    % due to the reduction of linear constraints by pre_lcon if some variables are fixed by bound
    % constraints, particularly when both sides of a linear constraint contain infinite values. A
    % well-defined problem should not have such constraints anyway. We consider them for robustness.
end

req = [];
if ~isempty(Aeq)
    req = Aeq*x - beq;
end

% max(X, [], 'includenan') returns NaN if X contains NaN, and maximum of X otherwise
cstrv = max([0; rineq; abs(req); rlb; rub; nlcineq; abs(nlceq)], [], 'includenan');
return
