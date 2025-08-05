function retmsg(info, iprint, nf, f, x, solver)

    % Local variables    % Should be an integer of default kind

    if iprint == 0
        return;
    end

    if info == infos('ftarget_achieved')
        msg = 'the target function value is achieved.';
    elseif info == infos('maxfun_reached')
        msg = 'the objective function has been evaluated MAXFUN times.';
    elseif info == infos('small_tr_radius')
        msg = 'the trust region radius reaches its lower bound.';
    elseif info == infos('trsubp_failed')
        msg = 'a trust region step has failed to reduce the quadratic model.';
    elseif info == infos('nan_x')
        msg = 'NaN occurs in x.';
    elseif info == infos('nan_inf_f')
        msg = 'the objective function returns NaN or +INFINITY.';
    elseif info == infos('nan_model')
        msg = 'NaN occurs in the models.';
    end

    if iprint >= 1
        if iprint >= 3
            fprintf('\n');
        end
        fprintf('Return from %s because %s\n', solver, deblank(msg));
        fprintf('At the return from %s    Number of function evaluations = %d\n', solver, nf);
        fprintf('Least function value = %.16f    The corresponding X is: %.16f\n', f, x);
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
            fprintf(OUTUNIT, 'Return from %s because %s\n', solver, deblank(msg));
            fprintf(OUTUNIT, 'At the return from %s    Number of function evaluations = %d\n', solver, nf);
            fprintf(OUTUNIT, 'Least function value = %.16f    The corresponding X is: %.16f\n', f, x);
            fclose(OUTUNIT);
        end
        %print '(/1A /)', 'The output is printed to ' // trim(fout) // '.'
    end

end
