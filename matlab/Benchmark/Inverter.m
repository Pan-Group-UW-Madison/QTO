clear;

system("clear");

addpath(genpath('../'));

params.dim = 2;

params.nelx = 600;
params.nely = 300;
params.Density = 0.3;

params.ClusterNelx = 50;
params.ClusterNely = 25;

params.rmin = 4;
params.epsilon = 5e-3;

% params.MilpSolver = 'milp';
% params.MilpSolver = 'lp';
params.MilpSolver = 'multilevel-lp';
params.PreStage = 1;

params.tolerance = 1e-9;

params.BC = 'inverter';
params.objective = 'compliant mechanism';
params.filter = 'radius';

params.NumMaterial = 5;

params.nu = 0.3;
if params.NumMaterial == 1
    params.E = 1.0;
    params.density = 1.0;
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
    params.Density0 = 0.3;
    params.N = 1;
    params.Emin = 1e-9;
end

params.maxFem = 100;

params.d0 = 0.2;
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
    fprintf('    Stage %3d: Time of Fem: %.4fs, Time of Opt: %.4fs\n', i, result.timeFemStage(i),...
        result.timeOptStage(i));
end

% visualize result
x = result.x;
Visualize(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '.png']);