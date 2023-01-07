function verisym(A, tol)
	% VERISYM verifies whether a matrix A is symmetric up to TOL.

	funname = 'VERISYM';

	if size(A, 1) ~= size(A, 2)
		error('Error: %s: A is not square.', deblank(funname));
	end
	if tol > 0
		if max(abs((A - transpose(A)))) > tol * max(max(abs(A)), 1)
			error('Error: %s: A is not symmetric up to TOL.', deblank(funname));
		end
	else
		if max(abs((A - transpose(A)))) > 0
			error('Error: %s: A is not symmetric up to TOL.', deblank(funname));
		end
	end
end
