function [avg_SAD,SADmin,Imin] = SAD_Evaluate(estA, refA, r) 

    SAD_tb = zeros(r, r); 
    for i = 1 : r
        for j = 1 : r
            SAD_tb(i,j) = SAD(estA(:, i), refA(:, j));
        end
    end

    [SADmin, Imin] = min(SAD_tb, [], 2);
    avg_SAD = mean(SADmin);
    
end