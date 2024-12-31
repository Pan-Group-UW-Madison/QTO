function d = UpdateTrustRegion(c0, c, cost, d, params)
    if (isempty(params.fixedD) == false && params.fixedD)
        return;
    end

    factors = zeros(length(c0), 1);
    for i = 1:length(c0)
        omega = (c0(i) - c) / (c0(i) - cost);
        if omega < 1 && omega >= 0
            factors(i) = min(d * 0.7);
            factors(i) = max(factors(i), 1e-3);
        elseif omega < 0
            factors(i) = min(d * 0.5);
            factors(i) = max(factors(i), 1e-3);
        else
            factors(i) = min(d * 1.5);
            if strcmp(params.objective, 'minimum compliance')
                factors(i) = min(factors(i), 1.0);
            else
                factors(i) = min(factors(i), 0.6);
            end
        end
    end

    d = min(factors);

    if params.verbose
        fprintf("    c: %.4f\n", c);
        for i = 1:length(cost)
            fprintf("    C0: %.4f, Cost: %.4f\n", c0(i), cost(i));
        end
        fprintf("    Trust region factor: %.4f, omega: %.4f\n", d, omega);
    end
end