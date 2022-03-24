function perfdata(solvers, options)

stamp = strcat(strjoin(solvers, '_'), '.', int2str(options.mindim), '_', int2str(options.maxdim), '.', options.type);
time = options.time;
if isfield(options, 'stamp') && (isa(options.stamp, 'char') || isa(options.stamp, 'string')) && ~isempty(options.stamp)
    time = strcat(time, '.', options.stamp);
end
data_dir = options.data_dir;
matfile = fullfile(data_dir, strcat(stamp, '.perfdata.mat'));

loaded = false;
% Check whether to reload data or calculate everything from scratch.
if (isfield(options, 'reload') && options.reload == true)
    % Use directly the data of the last computation for the same solvers and dimensions.
    stamp_with_wildcard = strcat(strjoin(solvers, '_'), '.', int2str(options.mindim), '_', int2str(options.maxdim), '.*', options.type, '*');
    matfile_with_wildcard = fullfile(data_dir, strcat(stamp_with_wildcard, '.perfdata.mat'));
    mlist = dir(matfile_with_wildcard);
    if ~isempty(mlist)
        matfile = fullfile(mlist(1).folder, mlist(1).name);
    end
    load(matfile, 'frec', 'fmin', 'pdim', 'plist');
    loaded = true;
    qlist = secup(options);
    fun = @(p) ismember(p, qlist);
    pind = find(cellfun(fun, plist));
    frec = frec(pind, :, :, :);
    fmin = fmin(pind);
    %pdim = pdim(pind);
    plist = plist(pind);
end
if ~loaded
    [frec, output] = testcu(solvers, options);
    fmin = min(min(min(frec, [], 4), [], 3), [], 2);
    pdim = output.pdim;
    plist = output.plist;
    save(matfile, 'frec', 'fmin', 'pdim', 'plist');
end

% `outdir` is the directory to contain the figures (.eps, .pdf) and problem list (problem.txt).
outdir = fullfile(data_dir, strcat(stamp, '.', time));
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Record the tested problems in `problems.txt`.
fprob = fullfile(outdir, strcat(stamp, '.', time, '.', 'problems.txt'));
fid = fopen(fprob, 'w');
if fid >= 3
    for k = 1 : length(plist)
        fprintf(fid, "%s\n", plist{k});
    end
    fclose(fid);
else
    error('\nFail to open file %s.\n', fprob);
end

% Plot the profiles.
prof_options = struct();
prof_options.solvers = solvers;
prof_options.outdir = outdir;
prof_options.stamp = stamp;
prof_options.time = time;
for tau = 10.^(-1:-1:-10)
    prof_options.tau = tau;
    perfprof(frec, fmin, prof_options);
    %dataprof(frec, fmin, pdim, prof_options);
end

% For convenience, save a copy of `problems.txt` and the figures in data_dir. They will be
% replaced in next test with the same `solvers` and `dimrange`.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the old files.
delete(fullfile(data_dir, strcat(stamp, '.*.problems.txt')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.pdf')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.eps')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
copyfile(fprob, data_dir);
epsfiles = dir(fullfile(outdir, '*.eps'));
for k = 1 : length(epsfiles)
    epsname = epsfiles(k).name;
    pdfname = strrep(epsname, '.eps', '.pdf');
    if exist(fullfile(outdir, pdfname), 'file')
        copyfile(fullfile(outdir, pdfname), data_dir);
    else
        copyfile(fullfile(outdir, epsname), data_dir);
    end
end
