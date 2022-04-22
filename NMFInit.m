function data = NMFInit(param, init_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the data structure with fields A and S.
% The initialization consists in alternating between a couple of unconstrained
% (fast) updates and constrained (more precise) updates.
%
% Inputs:
% - Y: data matrix.
% - rank: number of sources.
% - init_type: 1:random, 2:PCA, 3:OSP, 4:VCA
%
% Output: 
% - data with fields A and S such that Y \approx A * S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y = param.Y;
    r = param.rank;
    shape = param.S.sourcesShape;

    m = size(Y, 1);
    
    if(init_type == 1)      % random     
        data.A = rand(m, r);  
        
    elseif(init_type == 2)  % PCA
        R = Y * Y.';
        [V, D] = eig(R);
        [~, index] = sort(diag(D), 'descend');
        V_sort = V(:, index);
        A = V_sort(:, 1:r); 
        % normalize on mixing matrix A
        for k = 1 : r
            A(:, k) = A(:, k) / norm(A(:, k));
        end
        data.A = abs(A); % promise to non-negative of A   
        
    elseif(init_type == 3)  % OSP
        % OSP init
        HIM = Y';
        HIM = reshape(HIM,shape(1),shape(2),size(Y,1));
        [A, ~] = OSP(HIM, r);
        data.A = A;
        
    elseif (init_type == 4) % VCA
        [x1, ~, ~, ~] = dataProj(Y, r, 'proj_type', 'affine');
        [A_vca, ~, ~] = VCA(x1,'Endmembers',r,'verbose','off'); 
        data.A = A_vca;    
        
    else
        data.A = rand(m, r);
        
    end    
    
    
    % unconstrained updates and constrained updates
    options.MaximumIteration = 50;
    options.RelativeDifferenceTolerance = 0;
    
    for k = 1 : 2        
        % unconstrained updates
        i = sum(data.A) > 0;
        data.S(i, :) = max(data.A(:, i) \ Y, 0);
        i = sum(data.S, 2) > 0;
        data.A(:, i) = max(Y / data.S(i, :), 0);

        % constrained updates
        data.S = nonnegativeSparseUpdate(data.S, data.A'*data.A, data.A'*Y, 0, options);
        % reinitialize S if need be
        data = reinitializeNullSources(data, Y);
        data.A = nonnegativeSparseUpdate(data.A', data.S*data.S', data.S*Y', 0, options)';    
    end    
    
end

