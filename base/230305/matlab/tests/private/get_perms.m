function permutations = get_perms(m, n)
%GET_PERMS gets m random permutations of [1:n]. If factorial(n) > m, then the permutations are distinct
% from each other.
permutations = NaN(m, n);
if factorial(n) <= m
    permutations(1:factorial(n), :) = perms(1:n);
    for i = factorial(n) + 1 : m
        permutations(i, :) = permutations(mod(i-1, factorial(n)) + 1, :);
    end
else
    for i = 1 : m
        while(true)
            permutations(i, :) = randperm(n);
            if all(sum(abs(permutations(1:i-1, :) - permutations(i, :)), 2) > 0)
                break
            end
        end
    end
end
return
