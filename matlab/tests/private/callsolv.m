function [x, f, exitflag, output] = callsolv(solv, problem, thorough_test)

if (nargin <= 1)
    fprintf('solv and problem must be specified.\n');
    return;
end

[x, f, exitflag, output] = feval(solv, problem);

if (nargin == 3 && thorough_test == 1)
    type = problem.type;
    objective = problem.objective;
    x0 = problem.x0;
    A = problem.Aineq;
    b = problem.bineq;
    Aeq = problem.Aeq;
    beq = problem.beq;
    lb = problem.lb;
    ub = problem.ub;
    options = problem.options;
%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (strcmpi(solv, 'oldcobyla'))
        nonlcon = problem.nlc;
        prob.x0 = x0;
        prob.objective = objective;
        prob.Aineq = [];
        prob.bineq = [];
        prob.Aeq = [];
        prob.beq = [];
        prob.lb = [];
        prob.ub = [];
        prob.nonlcon = problem.nonlcon;
        prob.options = problem.options;
        problem.nonlcon = problem.nlc;
        problem = rmfield(problem, 'nlc');
%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        nonlcon = problem.nonlcon;
        prob = problem;
    end
    n = length(x0);
    m = 10*n;
    
    if (strcmpi(solv, 'uobyqa') || strcmpi(solv, 'newuoa') || strcmpi(solv, 'bobyqa') || strcmpi(solv, 'lincoa') || strcmpi(solv, 'cobyla') || strcmpi(solv, 'pdfo'))
        % Calls with options
        prob.solver = solv;
        [xx, ff, ee, oo] = pdfo(prob);
        oo
        if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        options.solver = solv;
        [xx, ff, ee, oo] = pdfo(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        oo
        if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        % Calls with empty options/without options
        prob.options = []; [x1, f1, exitflag1, output1] = feval(solv, prob);
        if (~iseq(x, f, exitflag, output, x1, f1, exitflag1, output1, 1e-5)) 
            where; keyboard;;
        end
        prob = rmfield(prob, 'options'); [xx, ff, ee, oo] = feval(solv, prob);
        oo
        if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        prob.options = []; [xx, ff, ee, oo] = pdfo(prob);
        oo
        if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        prob = rmfield(prob, 'options'); [xx, ff, ee, oo] = pdfo(prob);
        oo
        if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
    end
    
    % Solver-dependent calls 
    switch solv
    case {'uobyqa', 'newuoa'}
        switch type
        case 'u'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        otherwise
            fprintf('%s does not solve a problem of type %s.\n', solv, type);
        end
    case 'bobyqa'
        switch type
        case 'u'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set bounds to [] or inf
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, -inf(n,1), [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, -inf(n,1), inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, -inf(n,1), []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, -inf(n,1), inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'b'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, lb, ub, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0, lb, ub, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, lb, ub);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        otherwise
            fprintf('%s does not solve a problem of type %s.\n', solv, type);
        end
    case 'lincoa'
        switch type
        case 'u'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set A, b, Aeq, beq, lb, ub to [] (or inf for lb/ub)
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], [], [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], -inf(n,1), [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], [], inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], -inf(n,1), inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set A, Aeq to zero 
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), [], [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), -inf(n,1), [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), [], inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), -inf(n,1), inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set b to inf 
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), [], [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), -inf(n,1), [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), [], inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), -inf(n,1), inf(n,1), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], -inf(n,1), []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], [], inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], -inf(n,1), inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), -inf(n,1), []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), [], inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), -inf(n,1), inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), -inf(n,1), []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), [], inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), -inf(n,1), inf(n,1));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'b'
            % Calls with options
            % Set A, b, Aeq, beq to []
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], lb, ub, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set A, Aeq to zero 
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), lb, ub, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set b to inf
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), lb, ub, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            % Set A, b, Aeq, beq to []
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], lb, ub);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], lb, ub, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set A, Aeq to zero 
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), lb, ub);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, zeros(m,n), zeros(m,1), zeros(m,n), zeros(m,1), lb, ub, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set b to inf
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), lb, ub);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set b to inf
            [xx, ff, ee, oo] = feval(solv, objective, x0, rand(m,n), inf(m,1), zeros(m,n), zeros(m,1), lb, ub, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'l'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        otherwise
            fprintf('%s does not solve a problem of type %s.\n', solv, type);
        end
    case 'cobyla'
        problem.solver = 'cobyla';
        [xx, ff, ee, oo] = pdfo(problem);
        if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        problem.options = []; [xx, ff, ee, oo] = pdfo(problem); 
        if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        problem = rmfield(problem, 'options'); [xx, ff, ee, oo] = pdfo(problem);
        if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
            where; keyboard;;
        end
        switch type
        case 'u'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set nonlcon to [] 
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set nonlcon to "trivial constraints"
            [xx, ff, ee, oo] = feval(solv, objective, x0, @tricon, options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options 
            [xx, ff, ee, oo] = feval(solv, objective, x0, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set nonlcon to [] 
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Set nonlcon to "trivial constraints"
            [xx, ff, ee, oo] = feval(solv, objective, x0, @tricon, []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, @tricon);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'b'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, @(x) bcon(x, lb, ub), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0, @(x) bcon(x, lb, ub), []);
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, @(x) bcon(x, lb, ub));
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'l'
            % Calls with options
            [xx, ff, ee, oo] = feval(solv, objective, x0, @(x) lcon(x, A, b, Aeq, beq, lb, ub), options);
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 1e-7)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0, @(x) lcon(x, A, b, Aeq, beq, lb, ub), []) ;
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 1e-7)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, @(x) lcon(x, A, b, Aeq, beq, lb, ub) ) ;
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 1e-7)) 
                where; keyboard;;
            end
        case 'n'
            % Calls with options
            nonlcon_aug = @(x) allcon(x, A, b, Aeq, beq, lb, ub, nonlcon); % All constraints: cineq <= 0, ceq = 0
            [xx, ff, ee, oo] = feval(solv, objective, x0, nonlcon_aug, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 1e-4)) 
                where; keyboard;;
            end
            % Calls with empty options/without options
            [xx, ff, ee, oo] = feval(solv, objective, x0, nonlcon_aug, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 1e-4)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, nonlcon_aug);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 1e-4)) 
                where; keyboard;;
            end
        end
    case 'pdfo'
        switch type
        case 'u'
            [xx, ff, ee, oo] = feval(solv, objective, x0, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'b'
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], lb, ub, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], lb, ub);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, [], [], [], [], lb, ub, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'l'
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        case 'n'
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            oo
            if (~iseq(x, f, exitflag, output, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
            [xx, ff, ee, oo] = feval(solv, objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, []);
            oo
            if (~iseq(x1, f1, exitflag1, output1, xx, ff, ee, oo, 0)) 
                where; keyboard;;
            end
        end
    end
end
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
    where; keyboard;
    eq = 0;
end
if (prec == 0 && (exitflag ~= ee|| oo.funcCount ~= output.funcCount))
    eq = 0;
end

diff = max([abs(ff-f)/(1+abs(f)), norm(xx-x)/(1+norm(x)), abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation))]) 

return

function [cineq, ceq] = tricon(x)
    cineq = [];
    ceq = [];
return

function [cineq, ceq] = bcon(x, lb, ub)
    cineq = [lb(lb>-inf)-x(lb>-inf); x(ub<inf)-ub(ub<inf)];
    ceq = [];
return

function [cineq, ceq] = lcon(x, A, b, Aeq, beq, lb, ub)
    cineq = [A*x-b; lb(lb>-inf)-x(lb>-inf); x(ub<inf)-ub(ub<inf)];
    ceq = Aeq*x-beq;
return

function [cineq, ceq] = allcon(x, Aineq, bineq, Aeq, beq, bl, bu, nonlcon) % All constraints: cineq <= 0, ceq = 0
lcineq = [bl(bl>-inf)-x(bl>-inf); x(bu<inf)-bu(bu<inf)]; 
if(~isempty(Aineq))
    lcineq = [Aineq*x - bineq; lcineq];
end
if (isempty(Aeq))
    lceq = [];
else
    lceq = Aeq*x - beq;
end
if (isempty(nonlcon))
    nlcineq = [];
    nlceq = [];
else
    [nlcineq, nlceq] = nonlcon(x); % Nonliear constraints: nlcineq <= 0, nlceq = 0
end

cineq = [lcineq; nlcineq];
ceq = [lceq; nlceq];

function where()
    callstack = dbstack;
    invoker_file = callstack(2).file;
    invoker_name = callstack(2).name;
    invoker_line = callstack(2).line;
    fprintf('%s:%s:%d\n', invoker_file, invoker_name, invoker_line);
return
