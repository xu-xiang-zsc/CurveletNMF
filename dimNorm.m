function n = dimNorm(x, dim, norm_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provides the norm 1 or 2 or Inf along any given dimension of a matrix
%
% Inputs: 
% - x: the matrix.
% - dim: dimension of the computation.
% - type: type of the norm among 1, 2 or Inf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2
        dim = 1;
    end
    if nargin < 3
        norm_type = 2;
    end
    
    if norm_type == 2
        n = sqrt(sum(x .* x, dim));
    elseif norm_type == 1
        n = sum(abs(x), dim);
    elseif norm_type == Inf
        n = max(abs(x), [], dim); 
    else
        error('Authorized values for norm_type are: 1, 2 and Inf.\n');
    end
    
end




