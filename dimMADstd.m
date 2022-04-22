% sigmas = dimMADstd(X, dim)
% Provides an estimates of the standard deviations of "X" along dimension "dim"
% using the MAD (median absolute deviation) estimator.

function sigmas = dimMADstd(X, dim)
% compute the standard deviation of X based on the MAD estimator along dimension dim
    if nargin == 1
       dim = 1; 
    end
    
    medval = median(X, dim);
    sigmas = 1.4826 * median(abs(bsxfun(@minus, X, medval)), dim);

end



