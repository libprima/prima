function rhomsg(iprint, nf, f, rho, x, solver)

    % Local variables    % Should be an integer of default kind

    if iprint == 0
        return;
    end

    if iprint >= 2
        if iprint >= 3
            fprintf('\n');
        end
        fprintf('New RHO = %sNumber of function evaluations = %s\n', rho, nf);
        fprintf('Least function value = %sThe corresponding X is:%s\n', f, x);
    end

    if iprint <= -2
        fout = strcat(solver, '_output.txt');
        OUTUNIT = fopen(deblank(fout), 'a');
        if OUTUNIT < 0
            fprintf('Fail to open file %s%\n', deblank(fout));
        else
            if iprint <= -3
                fprintf(OUTUNIT, '\n');
            end
            fprintf(OUTUNIT, 'New RHO = %sNumber of function evaluations = %s\n', rho, nf);
            fprintf(OUTUNIT, 'Least function value = %sThe corresponding X is:%s\n', f, x);
            fclose(OUTUNIT);
        end
    end

end