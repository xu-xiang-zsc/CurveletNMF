function S = CP_Update_S(S0, AtA, AtY, lambda, parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves argmin_{S>=0} ||Y-A*S||_2^2 + ||lambda .* (S W^T)||_1 
% using Chambolle-Pock algorithm
%
% Inputs:
% - S0: starting point for S
% - AtA: matrix A' * A
% - AtY: matrix A' * Y
% - lambda: sparsity parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    if nargin < 4
        err('''lambda'' is required.\n');
    end
    if nargin < 5
        parameters = [];
    end

    parameters.dimensions = find([1, parameters.sourcesShape] > 1);
    s2l = @(S) reshape(S, [size(S, 1), numel(S) / size(S, 1)]);
    l2s = @(S) reshape(S, [size(S, 1), parameters.sourcesShape]);

   
    x = S0;
    
    for i = 1 : size(x, 1)
        temp = reshape(x(i,:,:), parameters.sourcesShape(1), parameters.sourcesShape(2));
        X_cur{i} = fdct_wrapping(temp, 1, 1, parameters.scales, parameters.angles); % Curvelet transform 
    end
    
    y_L1 = X_cur;
    y_pos = min(0, x);
    x1 = x;
    
    sigma = 0.675; 
    tau = 0.675;
    theta = 1; 

    HI1 = inv(tau * AtA + eye(size(AtA, 1)));
    tauAY = tau * AtY;


    %% main looping 
        
    for k = 1 : parameters.MaximumIteration

        sigX = l2s(sigma * x);
        for i = 1 : size(sigX,1)
            tp = reshape(sigX(i,:,:), parameters.sourcesShape(1), parameters.sourcesShape(2));
            sig_cur{i} = fdct_wrapping(tp, 1, 1, parameters.scales, parameters.angles); % Curvelet transform 
        end
        
        cfs = y_L1;
        for s = 1 : length(y_L1)
          for w = 1 : length(y_L1{s})
              for v = 1 : length(y_L1{s}{w}) 
                cfs{s}{w}{v} = (y_L1{s}{w}{v} + sig_cur{s}{w}{v});
              end
          end
        end
        
        y_L1 = curv_proj(cfs, lambda);

        y_pos = min(0, y_pos + sigma * x);

        Y_wt = zeros(length(y_L1), parameters.sourcesShape(1), parameters.sourcesShape(2));
        for i = 1 : length(y_L1)
            Y_wt(i,:,:) =  ifdct_wrapping(y_L1{i}, 1, parameters.sourcesShape(1), parameters.sourcesShape(2));
        end

        
        x2 = HI1 * s2l(x1 - tau * (Y_wt + y_pos) + tauAY);
        x = l2s((1 + theta) * x2 - theta * s2l(x1));
        
        %-----------------------------------------------------------------%
%         val = norm(al(x1-x),'fro')/norm(al(x),'fro');
%         fprintf(1, 'S # %i, res = %f\n', k, val);   
        %-----------------------------------------------------------------%
        
        x1 = l2s(x2);   
        
    end

    S = max(0, x);

end


function y = curv_proj(x, lambda)    
    y = x;
    for s = 1 : length(x)
          for w = 1 : length(x{s})
              for v = 1 : length(x{s}{w})
                  y_temp = y{s}{w}{v};
                  lambda_temp = lambda{s}{w}{v};
                  y_temp(y_temp > lambda_temp) = lambda_temp(y_temp > lambda_temp);
                  y_temp(y_temp < -lambda_temp) = -lambda_temp(y_temp < -lambda_temp);
                  y{s}{w}{v} = y_temp;
              end
          end
    end
end