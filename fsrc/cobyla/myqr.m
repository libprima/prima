maxdim = 1000;
sigma = 1e6;
nr = 10;
subnr = 10;
for ir = 1 : nr + 3
    ir
    if ir == nr + 1
        m = 0;
    elseif ir == nr + 1
        m = 1;
    elseif ir == nr + 3
        m = maxdim;
    else
        m = floor(rand*maxdim);
    end
    for subir = 1 : subnr + 4
        subir
        if subir == subnr + 1
            n = 0;
        elseif subir == subnr + 2
            n = 1;
        elseif subir == subnr + 3
            n = m;
        elseif subnr == subnr + 4
            n = 2*m;
        else
            n = floor(2*rand*m);
        end
        r = floor(2*rand*min(m, n));
        A = sigma*randn(m, r)*diag(sigma.^randn(r, 1))*randn(r, n);
        economy = (rand > 0.5);
        debugging = (rand > 0.2);
        pivote = (rand > 0.5);
        if pivote
            [Q, R, P] = qrfac(A, economy, debugging);
        else
            [Q, R] = qrfac(A, economy, debugging);
        end
    end
end
