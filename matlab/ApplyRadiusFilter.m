function sensitivity = ApplyRadiusFilter(sensitivity, params)
    [H, Hs] = RadiusFilter(params);
    s = H * (sensitivity(:) ./ Hs);
    sensitivity = reshape(s, params.nely, params.nelx);
end