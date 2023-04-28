function trrad = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)

    % Generic module

    % Current trust region radius    % Norm of current trust region step    % Ratio threshold for contraction    % Ratio threshold for expansion    % Contraction factor    % Expansion factor    % Reduction ratio

    if ratio <= eta1
        trrad = gamma1 * dnorm;
    elseif ratio <= eta2
        trrad = max(0.5 * delta, dnorm);
    else
        trrad = max(0.5 * delta, gamma2 * dnorm);
    end

    % For noisy problems, the following may work better.
    %if (ratio <= eta1) then
    %trrad = gamma1*dnorm
    %else if (ratio <= eta2) then  % Ensure TRRAD >= DELTA
    %trrad = delta
    %else  % Ensure TRRAD > DELTA with a constant factor
    %trrad = max(delta*(1.0_RP + gamma2)/2.0_RP, gamma2*dnorm)
    %end if

end