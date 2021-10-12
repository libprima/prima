nr = 100;

n = 500;

e = 0;

for ir = 1: nr
    A = randn(n,n);
    for i = 1 : n
        A(i, i+1:n) = 0;
        %A(i+1:n, i) = 0;
    end
    %B = tinv(A, 'ut');
    B = tinv(A, 'lt');
    e = max( e, norm(B'*A'-eye(n,n))/(max(abs(A), [], 'all')+max(abs(B), [], 'all')));
end
e


function B = tinv(A, matrix_type)
n = size(A, 1);

if (strcmpi(matrix_type, 'LT'))
    B = zeros(n,n);
    for i = 1:n
        B(i,i) = 1/A(i,i);
        B(i, 1:i-1)= -(A(i, 1:i - 1) / A(i, i))*B(1:i-1,1:i-1);
    end
end

if (strcmpi(matrix_type, 'UT'))
    B = zeros(n,n);
    for i = 1:n
        B(i,i) = 1/A(i,i);
        B(1:i-1, i)= -B(1:i-1,1:i-1)*(A(1:i - 1, i) / A(i, i));
    end
end

end
