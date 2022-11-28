function testperm(mindim, maxdim, type)

requirement.mindim = mindim;
requirement.maxdim = maxdim;
requirement.type = type;

switch type
case 'u'
    solvers = {'uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla'};
case 'b'
    solvers = {'bobyqa', 'lincoa', 'cobyla'};
case 'l'
    solvers = {'lincoa', 'cobyla'};
case 'n'
    solvers = {'cobyla'};
otherwise
end

locate_cutest();
plist = secup(requirement);

rng(42);
for ip = 1 : length(plist)
    pname = plist{ip}
    p = macup(pname);
    permutation = get_perms(1, length(p.x0));
    q = permprob(p, permutation);
    for is = 1 : length(solvers)
        solver = str2func([solvers{is}, 'n'])
        [x, fx] = solver(p);
        [y, fy] = solver(q);
        if (norm(x-y(permutation)) > 1e-3*max(1, norm(x)))
            solver
            pname
            x
            y(permutation)
            norm(x-y(permutation)) / max(1, norm(x))
            fx
            fy
            keyboard
        end
    end
end
