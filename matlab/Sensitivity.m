function sensitivity = Sensitivity(x, params, StiffnessMat)
    if strcmp(params.objective, 'minimum compliance')
        sensitivityBase = Compliance(params, StiffnessMat);
    else
        sensitivityBase = Mechanism(params, StiffnessMat);
    end

    if params.dim == 2
        xPhys = zeros(params.nely, params.nelx);
        for i = 1:params.NumMaterial
            xPhys = xPhys + params.E(i) * x(:, :, i);
        end
        xPhys = xPhys + params.Emin;

        xMask = zeros(params.nely, params.nelx);
        for i = 1:params.NumMaterial
            xMask(x(:, :, i) > params.Emin) = 1;
        end

        sensitivity = zeros(params.nely, params.nelx, params.NumMaterial);
        for i = 1:params.NumMaterial
            xRho = x(:, :, i);
            xFilter = zeros(params.nely, params.nelx);
            xFilter(xRho > 0) = params.E(i) - params.Emin;
            xFilter(xRho == 0) = xPhys(xRho == 0) * (params.E(i) - params.Emin);
            xFilter(xMask == 0) = xPhys(xMask == 0);

            s = sensitivityBase .* xFilter;
            if strcmp(params.filter, 'radius')
                s = params.H * (s(:) ./ params.Hs);
                sensitivity(:, :, i) = reshape(s, params.nely, params.nelx);
            elseif strcmp(params.filter, 'pde')
                s = params.TF' * (params.LF' \ (params.LF \ (params.TF * s(:))));
                sensitivity(:, :, i) = reshape(s, params.nely, params.nelx);
            end

            if strcmp(params.BC, 'cantilever')
                halfNely = floor(params.nely / 2);
                upperSensitivity = sensitivity(1:halfNely, :, i);
                flippedLowerSensitivity = flip(sensitivity(halfNely+1:end, :, i), 1);
                sensitivity(1:halfNely, :, i) = 0.5 * (upperSensitivity + flippedLowerSensitivity);
                sensitivity(halfNely+1:end, :, i) = sensitivity(halfNely:-1:1, :, i);
            end
        end
    end
end

function sensitivity = Compliance(params, StiffnessMat)
    if params.dim == 2
        sensitivity = reshape(sum((params.U(StiffnessMat.edofMat) * StiffnessMat.KE) .* params.U(StiffnessMat.edofMat), 2), params.nely, params.nelx);

        if strcmp(params.BC, 'cantilever')
            halfNely = floor(params.nely / 2);
            upperSensitivity = sensitivity(1:halfNely, :);
            flippedLowerSensitivity = flip(sensitivity(halfNely+1:end, :), 1);
            sensitivity(1:halfNely, :) = 0.5 * (upperSensitivity + flippedLowerSensitivity);
            sensitivity(halfNely+1:end, :) = sensitivity(halfNely:-1:1, :);
        end
    end
end

function sensitivity = Mechanism(params, StiffnessMat)
    if params.dim == 2
        U1 = params.U(:, 1);
        U2 = params.U(:, 2);
        sensitivity = -reshape(sum((U1(StiffnessMat.edofMat) * StiffnessMat.KE) .* U2(StiffnessMat.edofMat), 2), params.nely, params.nelx);
    end
end