function result = TSNMF_core(parameters)

% TSNMF algorithm --- sparse regularization in the Curvelet domain (version 1.0)
%
% Aims at solving :
% argmin_(A>=0, S>=0) 1/2*||Y-A*S||_2^2 + lambda*||S * W^T||_1
%
% Required fields of the parameters structure:
% - Y: data matrix.
% - rank: number of sources.
%
% Optional fields:
% - MaximumIteration: number of iterations of the main loop
% - phaseRatio: transition between decreasing thresholding phase
%   and refinement phase in percent of the iterations.
%
% Optional subfields for S:
% - S.MaximumIteration: number of iterations of the update of S
% - S.sourcesShape: shape of the 2D images ([width, height])
% - S.tau_MAD: constant coefficient for the final threshold computation (default: 1).
% - S.directSparsity: penalize the coarse scale or not (default: 1)
% - S.uniformFirstThreshold: use the same threshold for each source at the
% first iteration (default: 1). If the sources have very different dynamics, better set to 0.
% - S.scales: number of scales of Curvelet
% - S.angles: number of angles of Curvelet
% - coeff_len: total number of coefficients of Curvelet (depends on the
% S.scale and S.angles)
%
% Optional subfields for A:
% - A.MaximumIteration: number of iterations of the update of A
% - A.RelativeDifferenceTolerance: relative difference tolerance
%   for convergence of the Forward-Backward sub-iterations (default: 10^-9)
% - normConstrained: 1 to force the norm of the rows of S to be non-increasing, 0 otherwise (default: 0).
%
% Output: structure with fields:
% - A: the mixing matrix.
% - S: the sparse sources matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% ---------------------------init-------------------------- %%
    % check for required parameters
    if ~isfield(parameters, 'Y')
        error('matrix Y needed.\n')
    end
    if ~isfield(parameters, 'rank')
        error('number of sources needed.\n')
    end   
    
    % init A and S
    data = NMFInit(parameters, 3); % 2: PCA init     
    
    %---------------------------------------------------------------------%
    SAD_avg = SAD_Evaluate(data.A, parameters.reference.A, parameters.rank);
    fprintf(1, 'Init of SAD = %f\n',SAD_avg);
    %---------------------------------------------------------------------%    
    
    % prepare Curvelet operator
    parameters.S.dimensions = 2 : (length(parameters.S.sourcesShape) + 1);
    parameters.shape2line = @(S) reshape(S, [size(S, 1), size(S, 2) * size(S, 3)]);
    parameters.line2shape = @(S) reshape(S, [size(S, 1), parameters.S.sourcesShape]);       
    
    
    %% first threshold
    data = NMFnormalization(data, 'A');   
    data.lambda = [];
    
    % compute maximum on each scale
    S0 = data.A' * (data.A * data.S - parameters.Y);        
    % for each source, 2D-Curvelet is adopted
    for i = 1 : size(S0, 1)
        S0_T = reshape(S0(i,:), parameters.S.sourcesShape(1), parameters.S.sourcesShape(2));
        C_YW{i} = fdct_wrapping(S0_T, 1, 1, parameters.S.scales, parameters.S.angles); % Curvelet transform 
        [temp_lambda] = computeCurvScaleFun(C_YW{i});
        data.lambda  = [data.lambda; temp_lambda];
    end
    
    % uniformize among sources (this may not be helpful if the sources have very different dynamics)
    if parameters.S.uniformFirstThreshold
        data.lambda = bsxfun(@times, ones(parameters.rank, 1), max(data.lambda, [], 1));
    end

    data.lambda(:, 1) = (parameters.S.directSparsity ~= 0) * data.lambda(:, 1);
    
    data.S = 0 * data.S;

    %% other variables
    parameters.refinement_beginning = floor(parameters.phaseRatio * parameters.MaximumIteration);


    %% ---------------------- main looping ---------------------- %%
    
    Total_SADmin = zeros(parameters.MaximumIteration,parameters.rank+1);
    Total_Imin = zeros(parameters.MaximumIteration,parameters.rank);
    Total_rmse = zeros(parameters.MaximumIteration,1);
    %------------------------------------------------------------------   
        
    for iter = 1 : parameters.MaximumIteration
        data.iteration = iter;
        
        %% threshold
        Lambdas = genCurFullMatrix(data.lambda, parameters.S);
        
        if parameters.S.reweightedL1
            S_inv = data.A \ parameters.Y;
            data.S_inv = S_inv;
            for i = 1 : size(S_inv,1)
                temp = reshape(S_inv(i,:), parameters.S.sourcesShape(1), parameters.S.sourcesShape(2));
                S_w{i} = fdct_wrapping(temp, 1, 1, parameters.S.scale); % Curvelet transform 
                sig_sources(i,:) = computeCurvMAD(S_w{i});
            end
            Sigmas = genCurFullMatrix( sig_sources, parameters.S);

            for k = 1 : length(Lambdas)
                for j = 1 : length(Lambdas{k})
                    for t = 1 : length(Lambdas{k}{j})
                        tp = (Lambdas{k}{j}{t}) ./ (1 + (abs(S_w{k}{j}{t}) ./ Sigmas{k}{j}{t}).^2);
                        Lambdas{k}{j}{t} = tp;
                % Lambdas = Lambdas ./ (1 + (abs(S_w) ./ Sigmas).^2);
                    end
                end
            end
        end
        
        data.Lambdas = Lambdas;

        % update S 
        AtA = data.A' * data.A;
        AtY = parameters.line2shape(data.A' * parameters.Y);
        S = parameters.line2shape(data.S);
        S = CP_Update_S(S, AtA, AtY, data.Lambdas, parameters.S);
        data.S = parameters.shape2line(S);
        clear AtA AtY S    
        
        % reinitialize lines of S if need be
        data = reinitializeNullSources(data, parameters.Y);
        
        
        % update A
        % A column norms must be non-increasing in order to converge to a stationary point
        % this is only helpful in the refinement steps (since the cost function is then fixed)
        
        if data.iteration >= parameters.refinement_beginning
            parameters.A.normConstrained = 1;
        end
        data = NMFnormalization(data, 'S');
        SSt = data.S * data.S';
        SYt = data.S * parameters.Y';
        data.A = nonnegativeSparseUpdate(data.A', SSt, SYt, 0, parameters.A)';
        
        data = NMFnormalization(data, 'A');            
        
        
        % Update thresholds
        data = m_thresholdManagement(data, parameters); 
        
        %-----------------------------------------------------------------%
        
        
        [SAD_avg,SAD_min,Imin] = SAD_Evaluate(data.A, parameters.reference.A, parameters.rank);
%         fprintf(1, 'Iter# %i/%i SAD = %f,[%i %i %i %i],lambda=%f, sg=%f, ',...
%         data.iteration, parameters.MaximumIteration,SAD_avg, ...
%         Imin(1), Imin(2), Imin(3), Imin(4), data.lambda(1,2), parameters.S.tau_MAD * sigma_grad(1,2));    
        fprintf(1, '%i  %f  %f  %f  %f  %f  ', data.iteration, SAD_avg, SAD_min(1), SAD_min(2), SAD_min(3), SAD_min(4));    
        Total_SADmin(iter, 1:end-1) = SAD_min(1:end);
        Total_SADmin(iter, end) = SAD_avg;
        Total_Imin(iter, 1:end) = Imin(1:end);
        [f_S] = sunsal(data.A, parameters.Y, 'POSITIVITY','yes', 'ADDONE', 'yes');
        Total_rmse(iter,1) = rmse_f(parameters.Y, data.A*f_S);
        fprintf(1, '  %f\n', Total_rmse(iter,1));    
        %-----------------------------------------------------------------%
        
        if(data.iteration == 137)
            break;
        end
        
    end
    
    %% ---------------------terminate------------------------------- %%
    result.A = data.A;
    result.S = data.S;
    result.rmse = Total_rmse;
    result.SADmin = Total_SADmin;
    result.Imin = Total_Imin;
    
end


%% --------------------- Parameter update ----------------------- %%
function data = m_thresholdManagement(data, parameters)

    sigma_update_ratio = 0.02;

    % online estimation of the noise on the gradient, based on the MAD estimator of the residual
    if data.iteration == 1 % not necessary to update all the residues at once
        ind = 1 : size(data.A, 1);
        
        data.sigma_res = zeros(size(data.A, 1), parameters.coeff_len); % this value is related with the size of curvelet coefficients
        for k = ind
            res = data.A(k, :) * data.S - parameters.Y(k, :);
            res = reshape(res, parameters.S.sourcesShape(1), parameters.S.sourcesShape(2));   
            resW = fdct_wrapping(res, 1, 1, parameters.S.scales,parameters.S.angles); % Curvelet transform 
            data.sigma_res(k,:) = computeCurvMAD(resW);
        end
    else
        ind = randperm(size(data.A, 1));
        ind = ind(1 : floor(size(data.A, 1) * sigma_update_ratio + 0.99));
        
        indr = data.A(ind, :) * data.S - parameters.Y(ind, :);
        indr = reshape(indr, [numel(ind), parameters.S.sourcesShape]);  
        indrWt_down = [];
        for z = 1 : size(indr,1)
            temp = reshape(indr(z,:,:), parameters.S.sourcesShape(1), parameters.S.sourcesShape(2));
            indrWt{z} = fdct_wrapping(temp,1,1, parameters.S.scales,parameters.S.angles); % Curvelet transform    
            indrWt_down = [indrWt_down; computeCurvMAD(indrWt{z})];
        end
        data.sigma_res(ind,:) = indrWt_down;
    end

    sigma_grad = data.sigma_res(1 : size(data.S, 1), :);
    for k = 1 : size(data.S, 1)
        sigma_grad(k, :) = sqrt(sum(bsxfun(@times, data.A(:, k).^2, data.sigma_res.^2)));
    end
    

    %% decreasing lambda

    if data.iteration <= parameters.refinement_beginning

        %linear decrease to tau_MAD * sigma_grad when reaching the refinement steps
        data.lambda = max(parameters.S.tau_MAD * sigma_grad, data.lambda ...
            - 1 / (parameters.refinement_beginning + 1 - data.iteration)...
            * (data.lambda - parameters.S.tau_MAD * sigma_grad));
        
        data.lambda(:, 1) = (parameters.S.directSparsity ~= 0) * data.lambda(:, 1);
    end   

end