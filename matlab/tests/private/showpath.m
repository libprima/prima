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
    fprintf('\n%s: %s', solver, which(regexprep(solver, '_classical$', '')));
    % `regexprep` removes '_classical' in case 'solver' ends with it.
end
fprintf('\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n');
