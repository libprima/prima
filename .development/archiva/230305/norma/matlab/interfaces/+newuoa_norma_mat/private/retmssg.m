function retmssg(info, iprint, nf, f, x, solver)

    % Local variables    % Should be an integer of default kind

    if iprint == 0
        return;
    end

    if info == infos('ftarget_achieved')
        mssg = 'the target function value is achieved.';
    elseif info == infos('maxfun_reached')
        mssg = 'the objective function has been evaluated MAXFUN times.';
    elseif info == infos('small_tr_radius')
        mssg = 'the trust region radius reaches its lower bound.';
    elseif info == infos('trsubp_failed')
        mssg = 'a trust region step has failed to reduce the quadratic model.';
    elseif info == infos('nan_x')
        mssg = 'NaN occurs in x.';
    elseif info == infos('nan_inf_f')
        mssg = 'the objective function returns NaN or +INFINITY.';
    elseif info == infos('nan_model')
        mssg = 'NaN occurs in the models.';
    end

    if iprint >= 1
        if iprint >= 3
            fprintf('\n');
        end
        fprintf('Return from %s because %s\n', solver, deblank(mssg));
        fprintf('At the return from %sNumber of function evaluations = %d\n', solver, nf);
        fprintf('Least function value = %15.6fThe corresponding X is:%15.6f\n', f, x);
        fprintf('\n');
    end

    if iprint <= -1
        fout = strcat(solver, '_output.txt');
        OUTUNIT = fopen(deblank(fout), 'a');
        if OUTUNIT < 0
            fprintf('Fail to open file %s%\n', deblank(fout));
        else
            if iprint <= -3
                fprintf(OUTUNIT, '\n');
            end
            fprintf(OUTUNIT, 'Return from %s because %s\n', solver, deblank(mssg));
            fprintf(OUTUNIT, 'At the return from %sNumber of function evaluations = %d\n', solver, nf);
            fprintf(OUTUNIT, 'Least function value = %15.6fThe corresponding X is:%15.6f\n', f, x);
            fclose(OUTUNIT);
        end
        %print '(/1A /)', 'The output is printed to ' // trim(fout) // '.'
    end

end