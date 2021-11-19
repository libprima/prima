function [Q, Rdiag] =  qrexc(A, Q, Rdiag, i)

DEBUGGING = true;

m = size(A, 1);
n = size(A, 2);
EPS = eps;

% Postconditions;
if (DEBUGGING)
     assert(n >= 2, 'N >= 2');
     assert(size(Q, 1) == m && size(Q, 2) == m, 'SIZE(Q) == [m, m]');
    tol = max(1.0E-10, min(1.0E-1, 1.0E6 * EPS * n));
     assert(isorth(Q, tol), 'The columns of Q are orthonormal');
end

%====================!;
% Calculation starts !;
%====================!;

if (n < 2) % Should not happen.;
    return;
end

for k = i : n - 1
    i, k, n
    size(Rdiag)
    size(Q)
    size(A)
    hypt = sqrt(Rdiag(k + 1)^2 + (Q(:, k)'*A(:, k + 1))^2);
    G = planerot([Rdiag(k + 1); (Q(:, k)'*A(:, k + 1))]);
    Q(:, [k, k + 1]) = (Q(:, [k + 1, k])*G');
    Rdiag([k, k + 1]) = [hypt, (Rdiag(k + 1) / hypt) * Rdiag(k)];
end

%====================!;
%  Calculation ends  !;
%====================!;

% Postconditions;
if (DEBUGGING)
     assert(size(Q, 1) == m && size(Q, 2) == m, 'SIZE(Q) == [m, m]');
     assert(isorth(Q, tol), 'The columns of Q are orthonormal');
end

function is_orth = isorth(A, tol)

n = size(A, 2);

if (n > size(A, 1))
    is_orth = false;
else
    if (nargin == 2)
        is_orth = all(abs((A'*A) - eye(n)) <= max(tol, tol * max(abs(A))), 'all');
    else
        is_orth = all(abs((A'*A) - eye(n)) <= 0, 'all');
    end
end
