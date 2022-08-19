function cstrv = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq)
%GET_CSTRV calculates the (absolute) constraint violation of x. Note that `nlcineq` and `nlceq` are
% values of the nonlinear inequality and equality constraints, the inequality constraint being
% nlcineq <= 0. `nlcineq` and `nlceq` are optional, while all the other inputs are obligatory.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if nargin < 7
    % Private/unexpcted error
    error(sprintf('%s:InvalidNargin', funname), ...
    '%s: UNEXPECTED ERROR: at least 8 inputs expected.', funname);
end
if nargin < 8
    nlcineq = 0;
end
if nargin < 9
    nlceq = 0;
end
lb(isnan(lb)) = -inf; % Replace the NaN in lb with -inf
ub(isnan(ub)) = inf; % Replace the NaN in ub with inf
bineq(isnan(bineq)) = inf; % Replace the NaN in bineq with inf
if ~isempty(Aeq)
    nan_eq = isnan(sum(abs(Aeq), 2)) & isnan(beq); % NaN equality constraints
    Aeq = Aeq(~nan_eq, :); % Remove NaN equality constraints
    beq = beq(~nan_eq);
end
if isempty(lb)
    lb = -inf(size(x));
end
if isempty(ub)
    ub = inf(size(x));
end
rineq = [];
req = [];
if ~isempty(Aineq)
    rineq = Aineq*x-bineq;
end
if ~isempty(Aeq)
    req = Aeq*x-beq;
end
% max(X, [], 'includenan') returns NaN if X contains NaN, and maximum of X otherwise
cstrv = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
return
