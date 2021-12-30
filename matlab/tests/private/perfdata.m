function perfdata(solvers, options)

time = datestr(datetime(), 'yymmdd_HHMM');
stamp = strcat(strjoin(solvers, '_'), '.', int2str(options.mindim), '_', int2str(options.maxdim), '.', options.type);
datadir = options.datadir;
matfile = fullfile(datadir, strcat(stamp, '.perfdate', '.mat'));

% Check whether to reload data or calculate everything from scratch.
if (isfield(options, 'reload') && options.reload == true)
    % Use directly the data of the last computation for the same solvers and dimensions.
    load(matfile, 'frec', 'fmin', 'pdim', 'plist');
else
    [frec, output] = testcu(solvers, options);
    fmin = min(min(min(frec, [], 4), [], 3), [], 2);
    pdim = output.pdim;
    plist = output.plist;
    save(matfile, 'frec', 'fmin', 'pdim', 'plist');
end

% `outdir` is the directory to contain the figures (.eps, .pdf) and problem list (problem.txt).
outdir = fullfile(datadir, strcat(stamp, '.', time));
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Record the tested problems in `problems.txt`.
fprob = fullfile(outdir, strcat(stamp, '.', 'problems.txt'));
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
for tau = 10.^(-1:-1:-10)
    prof_options.tau = tau;
    perfprof(frec, fmin, prof_options);
    %dataprof(frec, fmin, pdim, prof_options);
end

% For convenience, save a copy of `problems.txt` and the figures in datadir. They will be
% overwritten in next test with the same `solvers` and `dimrange`.
copyfile(fprob, datadir);
epsfiles = dir(fullfile(outdir, '*.eps'));
for k = 1 : length(epsfiles)
    epsname = epsfiles(k).name;
    pdfname = strrep(epsname, '.eps', '.pdf');
    if exist(fullfile(outdir, pdfname), 'file')
        copyfile(fullfile(outdir, pdfname), datadir);
    else
        copyfile(fullfile(outdir, epsname), datadir);
    end
end
