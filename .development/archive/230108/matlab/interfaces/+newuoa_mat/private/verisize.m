function verisize(varargin)

    funname = 'VERISIZE';

    maxarg = 3; % Maximal number of inputs
    nvararg = length(varargin); % Number of inputs

    % Interpret the input.
    % Expected inputs: [x, n, m], yet input m may be omitted.
    if (nvararg < 2)
        error(sprintf('%s:TooFewInputs', funname), '%s: at least 2 input.', funname);
    elseif (nvararg == 2)
        x = varargin{1};
        n = varargin{2};
        m = 1;
    elseif (nvararg >= 3 && nvararg <= maxarg)
        x = varargin{1};
        n = varargin{2};
        m = varargin{3};
    else
        error(sprintf('%s:TooManyInputs', funname), '%s: at most %d inputs.', funname, maxarg);
    end

	if size(x, 1) ~= n
		error('Error: %s: SIZE(X, 1) /= N.', deblank(funname));
	end
	if size(x, 2) ~= m
		error('Error: %s: SIZE(X, 2) /= M.', deblank(funname));
	end
end