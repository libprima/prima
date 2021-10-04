function test(solver)

% Make CUTEst available. The following three lines should be configured to fit the installation of
% CUTEst on the current machine.
setenv('CUTEST', '~/local/cutesif/cutest');
setenv('MASTSIF', '~/local/cutesif/sif');
addpath('~/local/cutesif/mtools/msrc');

current_dir = cd();
neupdfo_dir = fileparts(fileparts(current_dir));
opdfo_dir = fullfile(neupdfo_dir, 'OPDFO');

mexopt = struct();
mexopt.debug = true;

solvern = [solver, 'n'];
cd(neupdfo_dir);
setup(solvern, mexopt);
cd(opdfo_dir);
setup(solver, mexopt)
cd(current_dir);

solvers = {[solver, 'n'], solver};
options = struct();
options.mindim = 1;
options.maxdim = 100;
options.nr = 20;
switch solver
    case {'uobyqa', 'newuoa'}
        options.type = 'u';
    case 'bobyqa'
        options.type = 'bu';
    case 'lincoa'
        options.type = 'lbu';
    otherwise
        options.type = 'nlbu';
end

assert(isequiv(solvers, options));

fprintf('\n\nThe test on %s is successful!\n\n', solver);
