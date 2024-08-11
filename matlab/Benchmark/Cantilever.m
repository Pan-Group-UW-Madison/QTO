clear;

system("clear");

addpath(genpath('../'));

params.dim = 2;

params.nelx = 600;
params.nely = 400;
params.Density = 0.3;

params.ClusterNelx = 15;
params.ClusterNely = 10;

params.rmin = 4;
params.epsilon = 1e-3;

% params.MilpSolver = 'milp';
% params.MilpSolver = 'lp';
params.MilpSolver = 'multilevel-lp';
params.PreStage = 1;

params.tolerance = 1e-9;

params.BC = 'cantilever';
params.objective = 'minimum compliance';
params.filter = 'radius';
% params.filter = 'pde';

params.NumMaterial = 5;

params.nu = 0.3;
if params.NumMaterial == 1
    params.E = 1.0;
    params.density = 1.0;
elseif params.NumMaterial == 4
    params.E = [0.43, 0.7, 0.94, 1.0];
    params.density = [0.3, 0.5, 0.8, 1.0];
else
    params.E = [0.43, 0.7, 0.85, 0.94, 1.0];
    params.density = [0.3, 0.5, 0.65, 0.8, 1.0];
end

if strcmp(params.BC, 'inverter')
    params.xSymmetric = false;
    params.ySymmetric = true;
else
    params.xSymmetric = false;
    params.ySymmetric = false;
end

if strcmp(params.BC, 'inverter')
    params.k1 = 0.1;
    params.k2 = 0.1;
else
    if params.NumMaterial == 1
        params.Density0 = 0.6;
        params.N = 10;
        params.Emin = 1e-4;
    else
        params.Density0 = 0.3;
        params.N = 1;
        params.Emin = 1e-9;
    end
end

params.maxFem = 100;

params.d0 = 0.3;
params.verbose = true;
params.visualizeStep = true;
params.visualizeLevel = false;

result = MultiCutsTopOpt(params);

% print result
fprintf('Number of Fem: %d\n', result.numFem);
fprintf('Time of Fem: %.4fs\n', result.timeFem);
fprintf('Time of Opt: %.4fs\n', result.timeOpt);

fprintf('Objective: %.4f\n', result.obj);

for i = 1:length(result.timeFemStage)
    fprintf('    Stage %3d: Time of Fem: %8.4fs, Time of Opt: %8.4fs\n', i, result.timeFemStage(i),...
        result.timeOptStage(i));
end

for i = 1:params.NumMaterial
    fprintf('    Material %d: %.4f\n', i, sum(sum(result.x(:, :, i))) / (params.nelx * params.nely));
end

% visualize result
x = result.x;
Visualize(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '.png']);