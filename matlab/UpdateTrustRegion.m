function d = UpdateTrustRegion(c0, c, cost, d, params)
    omega = min((c0 -c) ./ (c0 - cost));
    if omega < 1 && omega >= 0
        factor = min(d * 0.7);
        factor = max(factor, 1e-3);
    elseif omega < 0
        factor = min(d * 0.5);
        factor = max(factor, 1e-3);
    else
        factor = min(d * 1.5);
        factor = min(factor, 1.0);
    end

    d = factor;

    if params.verbose
        fprintf("    Trust region factor: %.4f, omega: %.4f\n", d, omega);
    end
end