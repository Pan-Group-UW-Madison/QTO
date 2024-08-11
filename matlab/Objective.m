function [obj, params] = Objective(x, params, StiffnessMat)
    if strcmp(params.objective, 'minimum compliance')
        [obj, params] = Compliance(x, params, StiffnessMat);
    else
        [obj, params] = Mechanism(x, params, StiffnessMat);
    end
end

function [obj, params] = Compliance(x, params, StiffnessMat)
    if params.dim == 2
        xPhys = zeros(params.nely, params.nelx);
        for i = 1:params.NumMaterial
            xPhys = xPhys + params.E(i) * x(:, :, i);
        end
        xPhys = xPhys + params.Emin;

        sK = reshape(StiffnessMat.KE(:) * xPhys(:)', 64 * params.nelx * params.nely, 1);
        K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

        params.U(params.freedofs) = K(params.freedofs, params.freedofs) \ params.F(params.freedofs);

        obj = params.U' * params.F;
    end
end

function [obj, params] = Mechanism(x, params, StiffnessMat)
    if params.dim == 2
        xPhys = zeros(params.nely, params.nelx);
        for i = 1:params.NumMaterial
            xPhys = xPhys + params.E(i) * x(:, :, i);
        end
        xPhys = xPhys + params.Emin;

        sK = reshape(StiffnessMat.KE(:) * xPhys(:)', 64 * params.nelx * params.nely, 1);
        K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

        K(params.din, params.din) = K(params.din,params.din)+params.k1;
        K(params.dout, params.dout) = K(params.dout,params.dout)+params.k2;

        params.U(params.freedofs, :) = K(params.freedofs, params.freedofs) \ params.F(params.freedofs, :);

        obj = params.U(params.dout, 1);
    end
end