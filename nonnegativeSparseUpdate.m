function S = nonnegativeSparseUpdate(S0, AtA, AtY, lambda, parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves problem
% argmin_(S >= 0) 1 / 2 * ||Y - A * S||_2^2 + ||lambda .* S||_1
% using FISTA (Beck & Teboulle 2009)
%
% Inputs:
% - S0: starting point for S
% - AtA: matrix A' * A
% - AtY: matrix A' * Y
% - lambda: sparsity parameters
%
% Optional fields for argument "parameter":
%    - MaximumIteration: maximum number of iteration of the algorithm.
%    - RelativeDifferenceTolerance: relative difference tolerance.
%    - normConstrained: 1 to force the norm of the rows of S to be non-increasing, 0 otherwise (default: 0).
%
% Output:
% S: solution of the problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if nargin < 5
        parameters = [];
    end

    % precomputation
    L = max(eig(AtA));

    % initializations
    S = S0;
    prev_S = S;
    t = 1;

    % operators
    gradient = @(s) AtA * s - AtY;
    proximal = @(x, threshold) max(x-threshold, 0);
    
    if isfield(parameters, 'normConstrained')
        if parameters.normConstrained == 1  
            norm_limits = dimNorm(S0, 2);
            proximal = @(x, threshold) normProjection(max(x-threshold, 0), norm_limits);
        end
    end

    for k = 1 : parameters.MaximumIteration
        t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
        w = (t - 1) / t_next;
        t = t_next;

        R = (1 + w) * S - w * prev_S;
        prev_S = S;
        S = proximal(R - gradient(R) / L, lambda / L);

        if norm(al(prev_S - S), 'fro') / norm(al(S), 'fro') <  parameters.RelativeDifferenceTolerance
            break;
        end
    end

end

