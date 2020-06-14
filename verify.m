function success = verify(solvers, options)

success = false;

if nargin < 1
    fprintf('\nSolvers must be specified.\n');
    return
end

if length(solvers) ~= 2 
    fprintf('\nThere should be two solvers.\n')
    return
end

if nargin == 1
    options = struct();
end

if isfield(options, 'prec')
    prec = options.prec;
else
    prec = 0;
end

test_options = struct();
test_options.debug = true;
test_options.chkfunval = true;

requirements = struct();
if (isfield(options, 'mindim'))
    requirements.mindim = options.mindim;
else
    requirements.mindim = 1;
end
if (isfield(options, 'maxdim'))
    requirements.maxdim = options.maxdim;
else
    requirements.maxdim = 20;
end
if (isfield(options, 'mincon'))
    requirements.mincon = options.mincon;
else
    requirements.mincon = 0;
end
if (isfield(options, 'maxcon'))
    requirements.maxcon = options.maxcon;
else
    requirements.maxcon = 100;
end
if (isfield(options, 'type'))
    requirements.type = options.type;
else
    requirements.type = 'u';
%    requirements.type = 'ubln';
end

% Supress the following warning
orig_warning_state = warning;
cellfun(@(solver) warning('off', [solver, ':Debug']), solvers);
cellfun(@(solver) warning('off', [solver, ':ChkFunval']), solvers);
cellfun(@(solver) warning('off', [solver, ':Classical']), solvers);
cellfun(@(solver) warning('off', [solver, ':ReviseX0']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownProbField']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownOption']), solvers);

plist = secup(requirements);

fprintf('\n')
for ip = 1 : length(plist)
    pname = plist{ip};
    prob = macup(pname);
    fprintf('%3d. \t%16s:\t', ip, pname);
    prob.options = test_options;
    [x1, fx1, exitflag1, output1] = feval(solvers{1}, prob);
    [x2, fx2, exitflag2, output2] = feval(solvers{1}, prob);
    decup(pname);
    if ~iseq(x1, fx1, exitflag1, output1, x2, fx2, exitflag2, output2, prec)
        fprintf('The solvers produce different results on %s.\n', pname);
        warning(orig_warning_state); % Restore the behavior of displaying warnings
        return
    else
        fprintf('Success\n');
    end
end

success = true;

return


function eq = iseq(x, f, exitflag, output, xx, ff, ee, oo, prec) 
eq = 1;

if ~isfield(output,'constrviolation')
    output.constrviolation = 0;
end
if ~isfield(oo,'constrviolation')
    oo.constrviolation = 0;
end

if (norm(xx-x)/(1+norm(x)) > prec || abs(ff-f)/(1+abs(f)) > prec || abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation)) > prec)
    eq = 0;
end
if (prec == 0 && (exitflag ~= ee|| oo.funcCount ~= output.funcCount))
    eq = 0;
end

%diff = max([abs(ff-f)/(1+abs(f)), norm(xx-x)/(1+norm(x)), abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation))]) 

return
