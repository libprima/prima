function prob = permprob(prob, permutation)
%PERMPROB returns the permuted version of `prob`  with the variables permuted by `permutation`, i.e.,
%prob(x(permutation)).
inv_perm(permutation) = (1 : length(permutation));
prob.objective =  @(x) prob.objective(x(permutation));
prob.x0 = prob.x0(inv_perm);
if ~isempty(prob.lb)
    prob.lb = prob.lb(inv_perm);
end
if ~isempty(prob.ub)
    prob.ub = prob.ub(inv_perm);
end
if ~isempty(prob.Aineq)
    prob.Aineq = prob.Aineq(:, permutation);
end
if ~isempty(prob.Aeq)
    prob.Aeq = prob.Aeq(:, permutation);
end

if ~isempty(prob.nonlcon)
    prob.nonlcon = @(x) prob.nonlcon(x(permutation));
end
return
