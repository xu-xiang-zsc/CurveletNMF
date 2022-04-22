function Coeff = genCurFullMatrix(lambda, params)

% generates a Threshold Matrix repeating the values in lambda for size "shape"
%
% Inputs:
% - lambda: values to repeat.
% - shape: sizes of the image data.
%
% Output:
% - Coeff: Cell of all bands of threshold matrix


    X = zeros(params.sourcesShape(1), params.sourcesShape(2));
    C = fdct_wrapping(X, 1, 1, params.scales, params.angles); % just want to obtain the size of Curvelet coefficients
    
    for i = 1 : size(lambda, 1)
        
        count = 1;
        for s = 1 : length(C)
          for w = 1 : length(C{s})
%             C{s}{w}(1:end) = lambda(i,count);
            C{s}{w}(1:end) = lambda(count);
            count = count + 1;
          end
        end
        
        Coeff{i} = C;
    end

end

