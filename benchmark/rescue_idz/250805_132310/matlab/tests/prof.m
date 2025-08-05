function outputfiles = prof(varargin)

time = datestr(datetime(), 'yymmdd_HHMM');

% Set up the directory to save the testing data, i.e., `data_dir`.
mfiledir = fileparts(mfilename('fullpath')); % Directory where this .m file resides.
data_dir = fullfile(mfiledir, 'testdata');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

options = struct();
for iin = 1 : nargin
    if isstruct(varargin{iin})
        options = varargin{iin};
        varargin = {varargin{1:iin-1}, varargin{iin+1:end}};
        break;
    end
end

options.time = time;
options.data_dir = data_dir;

output = cell(1, 4);


timerVal = tic;

argin = [varargin, {options}];
output{1} = profile(argin{:});  % The plain profiling

if any(cellfun(@(x)strcmpi(x, 'load'), varargin))
    toc;
    return;
end

if ~any(cellfun(@(x)strcmpi(x, 'ncp'), varargin))
    varargin = [varargin, {'ncp'}];
end

options.randomizex0 = eps;
argin = [varargin, {options}];
output{2} = profile(argin{:});
options = rmfield(options, {'randomizex0'});

options.perm = true;
argin = [varargin, {options}];
output{3} = profile(argin{:});
options = rmfield(options, {'perm'});

% Precision of function evaluation ~ 1.0e-7 (eps('single') = 1.1921e-07)
options.eval_options = struct('single', true);
argin = [varargin, {options}];
output{4} = profile(argin{:});
options = rmfield(options, {'eval_options'});

% Precision of function evaluation ~ 1.0e-6
options.eval_options = struct('dnoise', 1e-6);
argin = [varargin, {options}];
output{5} = profile(argin{:});

% Precision of function evaluation ~ 1.0e-5
options.eval_options = struct('signif', 5);
argin = [varargin, {options}];
output{6} = profile(argin{:});
options = rmfield(options, {'eval_options'});

% Precision of function evaluation ~ 1.0e-4
options.eval_options = struct('noise', 1e-4);
argin = [varargin, {options}];
output{7} = profile(argin{:});
options = rmfield(options, {'eval_options'});

% Precision of function evaluation ~ 1.0e-3
options.eval_options = struct('signif', 3);
argin = [varargin, {options}];
output{8} = profile(argin{:});
options = rmfield(options, {'eval_options'});


outputfiles = struct();
prob_types = fieldnames(output{1});
for ipt = 1 : length(prob_types)
    ptype = prob_types{ipt};
    [~, outputfile] = fileparts(output{1}.(ptype));
    outputfile = strrep(outputfile, '.plain.', '.');
    outputfile = fullfile(data_dir, [outputfile, '.pdf']);
    delete(regexprep(outputfile, 'summary.*$','*.pdf'));
    inputfiles = '';
    for iout = 1 : length(output)
        inputfiles = [inputfiles, ' ', output{iout}.(ptype)];
    end
    system(['bash ', fullfile(mfiledir, 'private', 'compdf'), ' ', inputfiles, ' -o ', outputfile]);
    outputfiles.(ptype) = outputfile;
    fprintf('\nSummary for problem type %s:\n\n%s\n\n', ptype, outputfile);
end

toc(timerVal);  % Use `timerVal` to specify which call of `tic` we are comparing with.
