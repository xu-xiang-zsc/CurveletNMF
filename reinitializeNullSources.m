function data = reinitializeNullSources(data, Y)

% Straightforward reinitialization of a line of S and a column of A
% by picking one column in the residue.
%
% Inputs:
% - data: structure with fields A and S.
% - Y: the data matrix.
%
% Output:
% - data : same structure as input, after reinitialization of null sources
%   Non-null sources are left untouched.


    if sum(sum(data.S, 2) == 0) > 0 || sum(sum(isnan(data.S))) > 0
        fprintf(1, 'Reinitialization of %i null line(s) in S.\n', sum(sum(data.S, 2) == 0));

        indices = ((sum(data.S, 2) == 0) + sum(isnan(data.S), 2) > 0) > 0;
        [data.A(:, indices), data.S(indices, :)] =  m_fastExtractNMF(Y - data.A * data.S, sum(indices));
    end

end



function [A, S] = m_fastExtractNMF(residual, r)

    if r > 0
        [m, n] = size(residual);
        A = zeros(m, r);
        S = zeros(r, n);
        
        for i = 1 : r
            residual = residual .* (residual > 0);
            if sum(residual(:)) ~= 0
                % compute square norm of residual to select maximum one
                res2 = sum(residual .* residual);
                j = find(res2 == max(res2));
                if ~isempty(j)
                    j = j(1);
                    if res2(j) > 0
                        %normalize maximum residual
                        A(:, i) = residual(:, j) / sqrt(res2(j));
                        %compute scalar product with the rest od the residual and keep only
                        %positive coefficients
                        S(i, :) = A(:, i)' * residual;
                        S(i, :) = max(S(i, :), 0);
                        %compute new residual
                        residual = residual - A(:, i) * S(i, :);
                    end
                end
            end
        end
    else
        A = [];
        S = [];
    end
end

