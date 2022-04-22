function data = NMFnormalization(data,normalized_matrix_str, norm_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
% - data with fields A (m x r) and S (r x n).
% - normalized_matrix_str: a string 
% if normalized_matrix_str == 'A': columns of A are normalized and the
% weigths are set on S.
% if normalized_matrix_str == 'S': columns of S are normalized and the
% weigths are set on A.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if nargin < 2
        normalized_matrix_str = 'A';
    end
    if nargin < 3
        norm_type = 2;
    end

    % norms
    Ns = dimNorm(data.S', 1, norm_type);
    Na = dimNorm(data.A, 1, norm_type);
    N = Na .* Ns;
    act = N > 0;

    % normalize active sources and mixtures
    if strcmp(normalized_matrix_str, 'S')
        data.S(act, :) = bsxfun(@times, data.S(act, :), 1 ./ Ns(act)');
        data.A(:, act) = bsxfun(@times, data.A(:, act), N(act) ./ Na(act));
    elseif strcmp(normalized_matrix_str, 'A')
        data.S(act, :) = bsxfun(@times, data.S(act, :), (N(act) ./ Ns(act))');
        data.A(:, act) = bsxfun(@times, data.A(:, act), 1 ./ Na(act));
    else
        error('Field normalized_matrix must be either string ''A'' or ''S''.\n');        
    end

    %sets to zero inactive sources and mixtures
    data.S(~act, :) = 0 * data.S(~act, :);
    data.A(:, ~act) = 0 * data.A(:, ~act);

end
    




