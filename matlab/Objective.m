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
            xPhys = xPhys + (params.E(i) - params.Emin) * x(:, :, i);
        end
        xPhys = xPhys + params.Emin;

        if strcmp(params.BC, 'Zhou and Rozvany')
            xPhys(params.fixedIndex) = 1;
        end

        sK = reshape(StiffnessMat.KE(:) * xPhys(:)', 64 * params.nelx * params.nely, 1);
        K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

        params.U(params.freedofs) = K(params.freedofs, params.freedofs) \ params.F(params.freedofs);

        obj = params.U' * params.F;
    else
        xPhys = zeros(params.nelz, params.nely, params.nelx);
        for i = 1:params.NumMaterial
            xPhys = xPhys + (params.E(i) - params.Emin) * x(:, :, :, i);
        end
        xPhys = xPhys + params.Emin;

        sK = reshape(StiffnessMat.KE(:) * xPhys(:)', 64 * params.nelx * params.nely * params.nelz, 1);
        K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

        params.U(params.freedofs) = K(params.freedofs, params.freedofs) \ params.F(params.freedofs);

        obj = params.U' * params.F;
    end
end

function [obj, params] = Mechanism(x, params, StiffnessMat)
    if params.dim == 2
        if params.useSuperResolution
            superX = zeros(params.nely * params.superResolutionRatio, params.nelx * params.superResolutionRatio, params.NumMaterial);
            for i = 1:params.NumMaterial
                superX(:, :, i) = kron(x(:, :, i), ones(params.superResolutionRatio, params.superResolutionRatio));
            end
            superX = superX + params.Emin;
    
            sK = reshape(StiffnessMat.KE(:) * superX(:)', 64 * params.nelx * params.nely * params.superResolutionRatio * params.superResolutionRatio, 1);
            K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

            K(params.din, params.din) = K(params.din,params.din)+params.k1;
            K(params.dout, params.dout) = K(params.dout,params.dout)+params.k2;

            params.superU(params.freedofs, :) = K(params.freedofs, params.freedofs) \ params.superF(params.freedofs, :);

            obj = params.superU(params.dout, 1);

            superX = zeros(params.nely * params.superResolutionRatio, params.nelx * params.superResolutionRatio, params.NumMaterial);
            for i = 1:params.NumMaterial
                superX(:, :, i) = kron(x(:, :, i), ones(params.superResolutionRatio, params.superResolutionRatio));
            end
            superX = superX + max(1e-4, params.Emin);

            superX = params.HM * (superX(:) ./ params.HsM);
            superX = reshape(superX, params.nely * params.superResolutionRatio, params.nelx * params.superResolutionRatio);

            sK = reshape(StiffnessMat.KE(:) * superX(:)', 64 * params.nelx * params.nely * params.superResolutionRatio * params.superResolutionRatio, 1);
            K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

            K(params.din, params.din) = K(params.din,params.din)+params.k1;
            K(params.dout, params.dout) = K(params.dout,params.dout)+max(5e-3, params.k2);

            params.superU(params.freedofs, :) = K(params.freedofs, params.freedofs) \ params.superF(params.freedofs, :);
        else
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
end