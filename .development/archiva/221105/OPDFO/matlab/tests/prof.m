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

options.eval_options = struct('single', true);
argin = [varargin, {options}];
output{4} = profile(argin{:});
options = rmfield(options, {'eval_options'});

%options.eval_options = struct('signif', 5);
%argin = [varargin, {options}];
%output{5} = profile(argin{:});
%options = rmfield(options, {'eval_options'});

%options.eval_options = struct('dnoise', 1e-5);
%argin = [varargin, {options}];
%output{6} = profile(argin{:});

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
    system(['bash ', fullfile(mfiledir, 'compdf'), ' ', inputfiles, ' -o ', outputfile]);
    outputfiles.(ptype) = outputfile;
    fprintf('\nSummary for problem type %s:\n\n%s\n\n', ptype, outputfile);
end

toc(timerVal);  % Use `timerVal` to specify which call of `tic` we are comparing with.
