function fmssg(iprint, nf, f, x, solver)

    % Local variables    % Should be an integer of default kind

    if iprint == 0
        return;
    end

    if iprint >= 3
        fprintf('Function number%sF = %sThe corresponding X is:%s\n', nf, f, x);
    end

    if iprint <= -3
        fout = strcat(solver, '_output.txt');
        fexist = exist(deblank(fout), 'file');
        if fexist
            fstat = 'old';
        else
            fstat = 'new';
        end
        OUTUNIT = fopen(deblank(fout), 'a');
        if OUTUNIT < 0
            fprintf('Fail to open file %s%\n', deblank(fout));
        else
            fprintf(OUTUNIT, 'Function number%sF = %sThe corresponding X is:%s\n', nf, f, x);
            fclose(OUTUNIT);
        end
    end

end