% function X = normProjection(X,scales)
% Projection of X on the set of matrices with i-th column (resp. row) L2 norm
% inferior or equal to the i-th element of the row (resp. column) vector
% "scales".

function X = normProjection(X, scales)
    [m, n] = size(X);
    if size(scales, 1) == m
        n = dimNorm(X, 2);
        i = n > 0;
        X(i, :) = diag(min(n(i, :), scales(i, :)) ./ n(i, :)) * X(i, :);

    else
        n = dimNorm(X, 1);
        i = n > 0;
        X(:, i) = X(:, i) * diag(min(n(:, i), scales(:, i)) ./ n(:, i));
    end

end



