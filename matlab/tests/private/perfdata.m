function perfdata(solvers, options)

time = datestr(datetime(), 'yymmddHHMM');
stamp = strcat(strjoin(solvers, '_'), '.', int2str(options.mindim), '_', int2str(options.maxdim), '.', time);
matname = strcat(stamp, '.perfdate', '.mat');

if (isfield(options, 'load') && options.load == true)
    % Use directly the data of the last computation for the same n and feature.
    load(strcat('perfdata.', feature, '.', int2str(n), '.mat'), 'frec', 'fmin', 'n', 'feature');
else
    frec = testcu(solvers, options);
    [np, ~, ~, ~] = size(frec);
    fmin = NaN(np,1);
    for ip = 1:np
        fmin(ip) = min(min(min(frec(ip, :, :, :))));
    end
    save(strcat('perfdata.', feature, '.', int2str(n), '.mat'), 'frec', 'fmin', 'n', 'feature');
end

for tau = 10.^(-1:-1:-10)
    perfprof(frec, fmin, tau, solvers);
    %dataprof(frec, fmin, tau, n, feature);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test_options = struct();
%if (nargin == 2 && isfield(options, 'reproduce')) % Set options.reproduce = true for reproducing results of the last experiment.
%    test_options.reproduce = options.reproduce;
%else
%    test_options.reproduce = false;
%end

%if (nargin == 1)
%    options = [];
%end

%if (isfield(options, 'load') && options.load == true)
%    % Use directly the data of the last computation for the same n and feature.
%    load(strcat('perfdata.', feature, '.', int2str(n), '.mat'), 'frec', 'fmin', 'n', 'feature');
%else
%    switch feature
%    case {'Lq', 'Lh', 'L1'}
%        if (isfield(options, 'regul_lambda'))
%            test_options.regul.lambda = regul_lambda;
%        else
%            test_options.regul.lambda = 10;
%        end
%        if (strcmp(feature, 'Lq'))
%            test_options.regul.p = 0.25;
%        elseif (strcmp(feature, 'Lh'))
%            test_options.regul.p = 0.5;
%        elseif (strcmp(feature, 'L1'))
%            test_options.regul.p = 1;
%        end
%    case 'noisy'
%        if (isfield(options, 'noise_type'))
%            test_options.noise.type = options.noise_type;
%        else
%            test_options.noise.type = 'relative';
%        end
%        if (isfield(options, 'noise_level'))
%            test_options.noise.level = options.noise_level;
%        else
%            test_options.noise.level = 1e-3;
%        end
%        if (isfield(options, 'randrun'))
%            test_options.randrun = options.randrun;
%        end
%    case {'single', 's'}
%        test_options.prec = 'single';
%    case {'signif1', 'signif2', 'signif3', 'signif4', 'signif5', 'signif6', 'signif7', 'signif8'}
%        test_options.signif = str2num(feature(length(feature)));
%    case 'randomizex0'
%        if (isfield(options, 'randomizex0'))
%            test_options.randomizex0 = options.randomizex0;
%        else
%            test_options.randomizex0 = 1;
%        end
%        if (isfield(options, 'randrun'))
%            test_options.randrun = options.randrun;
%        end
%    otherwise
%        test_options = [];
%    end

%    test_options.debug=true;
%    test_options.chkfunval=true;

%    frec = testcu(solvers, 'ALL', test_options);

%    switch feature
%       case {'noisy', 'single', 's', 'signif1', 'signif2', 'signif3', 'signif4', 'signif5', 'signif6', 'signif7', 'signif8', 'randomizex0'}
%           display('=======================================================');
%           display('Auxiliary Computation:');
%           frec_aux = testcu(solvers, 'ALL', []);
%           [np, ~, ~, ~] = size(frec_aux);
%           fmin = NaN(np,1);
%           for ip = 1:np
%               fmin(ip) = min(min(min(min(frec(ip, :, :, :)))), min(min(min(frec_aux(ip, :, :, :)))));
%           end
%       otherwise
%           [np, ~, ~, ~] = size(frec);
%           fmin = NaN(np,1);
%           for ip = 1:np
%               fmin(ip) = min(min(min(frec(ip, :, :, :))));
%           end
%    end

%   % save(strcat('perfdata.', feature, '.', int2str(n), '.mat'), 'frec', 'fmin', 'n', 'feature');
%end

%for tau = 10.^(-1:-1:-10)
%    perfprof(frec, fmin, tau, solvers);
%    %dataprof(frec, fmin, tau, n, feature);
%end

%return;
