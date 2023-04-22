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
    nlcineq = 0;
end
if nargin < 9
    nlceq = 0;
end

if isempty(lb)
    lb = -inf(size(x));
end
if isempty(ub)
    ub = inf(size(x));
end

rlb = lb - x;
rlb(lb <= x) = 0;
rub = x - ub;
rub(x <= ub) = 0;

rineq = [];
req = [];
if ~isempty(Aineq)
    rineq = Aineq*x-bineq;
    %rineq(Aineq*x <= bineq) = 0;
    %rineq(bineq >= inf & ~isnan(Aineq*x )) = 0;
end
if ~isempty(Aeq)
    req = Aeq*x-beq;
    %req (Aeq*x == beq) = 0;
end
% max(X, [], 'includenan') returns NaN if X contains NaN, and maximum of X otherwise
cstrv = max([0; rineq; abs(req); rlb; rub; nlcineq; abs(nlceq)], [], 'includenan');
return
