clear;

system("clear");

params.dim = 2;

params.nelx = 8;
params.nely = 4;

params.Density = 0.3;

params.rmin = 2;

params.epsilon = 5e-3;

% params.BC = 'inverter';
% params.objective = 'compliant mechanism';

params.BC = 'cantilever';
params.objective = 'minimum compliance';

params.NumMaterial = 1;

params.nu = 0.3;
params.E = 1.0;
params.density = 1.0;

params.k1 = 0.1;
params.k2 = 0.1;

params.Emin = 1e-2;

params.mass = params.Density;

addpath(genpath('../'));

params = BCs(params);

params.alldofs = 1:2 * (params.nely + 1) * (params.nelx + 1);
params.freedofs = setdiff(params.alldofs, params.fixeddofs);

x = params.Density * ones(params.nely, params.nelx);

StiffnessMat = ElementMat(params);

xPhys = zeros(params.nely, params.nelx);
for i = 1:params.NumMaterial
    xPhys = xPhys + params.E(i) * x(:, :, i);
end
xPhys = xPhys + params.Emin;

if strcmp(params.objective, 'minimum compliance')
    sK = reshape(StiffnessMat.KE(:) * xPhys(:)', 64 * params.nelx * params.nely, 1);
    K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

    params.U(params.freedofs) = K(params.freedofs, params.freedofs) \ params.F(params.freedofs);
else
    sK = reshape(StiffnessMat.KE(:) * xPhys(:)', 64 * params.nelx * params.nely, 1);
K = sparse(StiffnessMat.iK, StiffnessMat.jK, sK); K = (K + K') / 2;

K(params.din, params.din) = K(params.din,params.din)+params.k1;
K(params.dout, params.dout) = K(params.dout,params.dout)+params.k2;

params.U(params.freedofs, :) = K(params.freedofs, params.freedofs) \ params.F(params.freedofs, :);
end

obj = params.U' * params.F;

alldofs = params.alldofs;
fixeddofs = params.fixeddofs;
freedofs = params.freedofs;

F = params.F;
U = params.U;

save('SingleMaterialCostFunc.mat', 'alldofs', 'fixeddofs', 'freedofs', 'K', 'F', 'U');