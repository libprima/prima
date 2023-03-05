function ndftpath = showpath(solvers)

fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
fprintf('\nCurrently in %s\n', pwd());

fprintf('\nNon-default paths:\n');
oldpath = path();
restoredefaultpath;
dftpath = path();
path(oldpath);
ndftpath = setdiff(split(oldpath, ':'), split(dftpath, ':'));
for ic = 1 : length(ndftpath)
    fprintf('\n%s', ndftpath{ic});
end

fprintf('\n\nSolver paths:\n');
for isol = 1 : length(solvers)
    solver = solvers{isol};
    % `regexprep` removes '_classical' in case 'solver' ends with it.
    solver = regexprep(solver, '_classical$', '');
    solver = regexprep(solver, '_single$', '');
    solver = regexprep(solver, '_quadruple$', '');
    fprintf('\n%s: %s', solvers{isol}, which(solver));
end
fprintf('\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n');
