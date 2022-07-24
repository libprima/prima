%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cstrv = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq)
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
if nargin < 8
    nlcineq = 0;
    nlceq = 0;
end
cstrv = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
% max(X, [], 'includenan') returns NaN if X contains NaN, and maximum of X otherwise
return
