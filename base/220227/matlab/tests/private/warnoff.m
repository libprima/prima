function orig_warning_state = warnoff(solvers)
orig_warning_state = warning;

for isol = 1 : length(solvers)
    if endsWith(solvers{isol}, '_classical')
        solvers{isol} = regexprep(solvers{isol}, '_classical', '');
        warning('off', [solvers{isol}, ':Classical']);
    end
end

solvers = [solvers, 'pdfo', 'pdfon'];

cellfun(@(solver) warning('off', [solver, ':Debug']), solvers);
cellfun(@(solver) warning('off', [solver, ':ChkFunval']), solvers);
cellfun(@(solver) warning('off', [solver, ':ReviseX0']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownProbField']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownOption']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidMaxfun']), solvers);
cellfun(@(solver) warning('off', [solver, ':ExtremeBarrier']), solvers);
cellfun(@(solver) warning('off', [solver, ':IprintContradictFortran']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidMaxhist']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidNpt']), solvers);
cellfun(@(solver) warning('off', [solver, ':ObjectiveFailure']), solvers);
cellfun(@(solver) warning('off', [solver, ':ConstraintFailure']), solvers);

return
