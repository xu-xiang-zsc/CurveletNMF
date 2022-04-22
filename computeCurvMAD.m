function cfs = computeCurvMAD(C)

    cfs = [];
    for s = 1 : length(C)
      for w = 1 : length(C{s})
        cfs = [cfs, dimMADstd(C{s}{w}(:))];
      end
    end
    
end