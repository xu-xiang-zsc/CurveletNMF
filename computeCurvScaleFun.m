% calculate the max value of each scale and each direction for Curvelet
% Coefficients

function [cfs] = computeCurvScaleFun(C)

    cfs = [];
    for s = 1:length(C)
      for w = 1:length(C{s})
          max_v = max(abs(C{s}{w}(:)));
          cfs = [cfs, max_v];
      end
    end
    
end