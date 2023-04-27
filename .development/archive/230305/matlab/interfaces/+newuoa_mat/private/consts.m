function const = consts(name)

    switch name
        % Some default values
        case 'rhobeg_dft'
            const = 1;
        case 'rhoend_dft'
            const = 1.0e-06;
        case 'ftarget_dft'
            const = -Inf;
        case 'iprint_dft'
            const = 0;
        case 'maxfun_dim_dft'
            const = 500;
        case {'funcmax', 'constrmax'}
            const = 10^40;

        % Maximal amount of memory (Byte) allowed for XHIST, FHIST, NLCHIST
        case 'maxmemory'
            const = 2.1e+09;
        otherwise
            error('Error (consts.m): unknown constant %s', name);
    end
end
