function [tlist, plist] = timing(solver, mindim, maxdim)

locate_cutest();

req.mindim = mindim;
req.maxdim = maxdim;
req.maxcon = 500*maxdim;
switch solver
    case {'uobyqa', 'newuoa'}
        req.type = 'u';
    case 'bobyqa'
        req.type = 'ub';
    case 'lincoa'
        req.type = 'ubl';
    case 'cobyla'
        req.type = 'ubln';
end

solver_fun = str2func([solver, 'n']);
plist = secup(req);
tlist= NaN(length(plist), 1);

for ip = 1 : length(plist)
    plist{ip}
    tic;
    prob = macup(plist{ip});
    prob.options.debug = true;
    solver_fun(prob);
    tlist(ip) = toc;
    tlist(ip)
end

save([solver, '_', int2str(mindim), '_', int2str(maxdim), '_', 'timing.mat'], 'plist', 'tlist');
