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
    for ip = 1 : length(plist)
        fprintf(fid, "%s\n", plist{ip});
    end
    fclose(fid);
else
    error('\nFail to open file %s.\n', fprob);
end

% Plot the profiles.
ns = length(solvers);
for is = 1:ns
    solvers{is} = regexprep(solvers{is}, '_4test', '');
    solvers{is} = regexprep(solvers{is}, '_classical$', ' (classical)');
    solvers{is} = regexprep(solvers{is}, '_single$', ' (single)');
    solvers{is} = regexprep(solvers{is}, '_quadruple$', ' (quadruple)');
    %solvers{is} = regexprep(solvers{is}, 'newuoa', 'NEWUOA');
end

prof_options = struct();
prof_options.solvers = solvers;
prof_options.outdir = outdir;
prof_options.stamp = stamp;
prof_options.time = time;

prec = (1:10)
nprec = length(prec)
tau = 10.^(-prec);
prof_output = repmat({struct()}, 1, nprec)
size(prof_output)
for iprec = 1 : nprec
    iprec
    prof_options.tau = tau(iprec);
    prof_output{iprec} = perfprof(frec, fmin, prof_options);
end

hfig = figure("visible", false, 'DefaultAxesPosition', [0, 0, 1, 1]);
for iprec = 1 : nprec
    %dataprof(frec, fmin, pdim, prof_options);
    subplot(1, nprec, iprec);
    for is = 1 : ns
       plot(prof_output{iprec}.profile{is}(1,:), prof_output{iprec}.profile{is}(2,:));
       hold on;
    end
    xlabel(sprintf('%d', iprec));
    axis([0 prof_output{iprec}.cut_ratio 0 1]);
    grid on;
end
for iprec = 1 : nprec
    ha=get(gcf,'children');
    set(ha(iprec),'position', [0.01+(nprec-iprec)/nprec, 0.1, 0.9/nprec, 0.9]);
end

% The following appears only in the last subplot, but it is sufficient for our use.
if isfield(options, 'stamp')
    ylabel(options.stamp);
end
legend(solvers,'Location', 'southeast','Orientation','vertical');

% Save the figure as eps.
figname = strcat(stamp, '.', 'summary', '.', time);
epsname = fullfile(outdir, strcat(figname,'.eps'));
set(gcf,'position',[0, 0, 4600,460]);
saveas(hfig, epsname, 'epsc2');

% Try converting the eps to pdf.
try
    system(['epstopdf ',epsname]);
catch
    % Do nothing in case of failure.
end


% For convenience, save a copy of `problems.txt` and the figures in data_dir. They will be
% replaced in next test with the same `solvers` and `dimrange`.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the old files.
delete(fullfile(data_dir, strcat(stamp, '.*.problems.txt')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.pdf')));
delete(fullfile(data_dir, strcat(stamp, '.perf_*.eps')));
delete(fullfile(data_dir, strcat(stamp, '.summary.*.eps')));
delete(fullfile(data_dir, strcat(stamp, '.summary.*.pdf')));
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
