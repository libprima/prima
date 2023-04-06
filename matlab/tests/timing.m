function [tlist, plist] = timing(solver, mindim, maxdim, blacklist, minip)

% Find where MatCUTEst is installed
locate_matcutest();

% Mute the warnings
orig_warning_state = warnoff({solver});

% Prepare the solver
cpwd = pwd();  % The current directory
mfiledir = fileparts(mfilename('fullpath'));  % The directory containing this m file
prima_dir = fileparts(fileparts(mfiledir));  % The directory containing the setup script
cd(prima_dir);
clear('setup');  % Without this, the next line may not call the latest version of `setup`
setup(solver);
cd(cpwd);

% In case blacklist and minip are not present ...
if nargin <= 3
    blacklist = {};
end
if nargin <= 4
    minip = 1;
end

% Define the problem list
req.mindim = mindim;
req.maxdim = maxdim;
req.maxcon = 100*maxdim;
switch solver
    case {'uobyqa', 'newuoa'}
        req.type = 'u';
    case 'bobyqa'
        req.type = 'b';
    case 'lincoa'
        req.type = 'l';
    case 'cobyla'
        req.type = 'ln';
end

solver_fun = str2func(solver);
plist = setdiff(secup(req), blacklist);
tlist= NaN(length(plist), 1);

for ip = minip : length(plist)
    fprintf("%4d. %s ", ip, plist{ip});
    tic;
    prob = macup(plist{ip});
    prob.options.debug = true;
    solver_fun(prob);
    decup(prob);
    tlist(ip) = toc;
    fprintf("\t\t%g\n", tlist(ip));
end

% Restore the behavior of displaying warnings
warning(orig_warning_state);

% Save the data
data_dir = fullfile(mfiledir, 'testdata');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
matfile = fullfile(data_dir, [solver, '_', int2str(mindim), '_', int2str(maxdim), '_', 'timing.mat']);
save(matfile, 'plist', 'tlist', '-v7.3');
