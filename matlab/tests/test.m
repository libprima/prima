function test(solver)

solvern = [solver, 'n'];

current_dir = cd();
neupdfo_dir = fileparts(fileparts(current_dir));
opdfo_dir = fullfile(neupdfo_dir, 'OPDFO');

mexopt = struct();
mexopt.debug = true;

cd(neupdfo_dir);
setup(solvern, mexopt);
cd(opdfo_dir);
setup(solver, mexopt)
cd(current_dir);

solvers = {[solver, 'n'], solver};
options = struct();
options.mindim = 1;
options.maxdim = 5;
options.nr = 1;
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
