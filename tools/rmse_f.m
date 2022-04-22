% RMSE calculator
% original and reconstruct  are --bands*N-- format

function rmses = rmse_f(original, reconstruct)

    [bands, N] = size(original);

    fdd = zeros(N, 2);
    rmses = 0.0;

    for i=1:N
        sum1 = 0.0;
    
        A = isnan(reconstruct(:,i));
        if A(1) == 0
            for j=1:bands
                sum1 = sum1 + (original(j,i)-reconstruct(j,i))*(original(j,i)-reconstruct(j,i));        
            end
            sum1 = sum1 / bands;
            sum1 = sqrt(sum1);
        end
    
    
        rmses = rmses + sum1;
        fdd(i,1) = i;
        fdd(i,2) = rmses;
        
%         fprintf('i = %d, rmses = %4f  \n', i, rmses);
        
    end

    rmses = rmses / N;

