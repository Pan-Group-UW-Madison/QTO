clear;

system("clear");

addpath(genpath('../'));

params.dim = 2;

params.nelx = 200;
params.nely = 100;
params.Density = 0.3;

params.ClusterNelx = 50;
params.ClusterNely = 25;

params.rmin = 2;
params.epsilon = 5e-3;

params.MilpSolver = 'milp';
% params.MilpSolver = 'lp';
% params.MilpSolver = 'multilevel-lp';
params.PreStage = 1;

params.singleCut = false;
params.fixedD = false;
params.storeIterResult = false;
params.compareMilp = false;
params.useQuantum = false;
params.useSuperResolution = false;
params.superResolutionRatio = 2;

params.tolerance = 1e-9;

params.BC = 'inverter';
params.objective = 'compliant mechanism';
params.filter = 'radius';

params.NumMaterial = 1;

params.nu = 0.3;
if params.NumMaterial == 1
    params.E = 1.0;
    params.density = 1.0;
else
    params.E = [0.43, 0.7, 0.85, 0.94, 1.0];
    params.density = [0.3, 0.5, 0.65, 0.8, 1.0];
    % params.E = [0.5, 1.0];
    % params.density = [0.3, 1.0];
    % params.E = [0.1, 1.0];
    % params.density = [1.0, 1.0];
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
    params.Emin = 1e-4;
else
    params.Density0 = 0.3;
    params.N = 1;
    params.Emin = 1e-9;
end

params.maxFem = 200;

params.d0 = 0.3;
params.verbose = true;
params.visualizeStep = false;
params.visualizeLevel = false;

result = MultiCutsTopOpt(params);

% print result
fprintf('Number of Fem: %d\n', result.numFem);
fprintf('Time of Fem: %.4fs\n', result.timeFem);
fprintf('Time of Opt: %.4fs\n', result.timeOpt);

params.useSuperResolution = false;
params.Emin = 1e-9;
StiffnessMat = ElementMat(params);

params = BCs(params);

params.alldofs = 1:2 * (params.nely + 1) * (params.nelx + 1);
params.freedofs = setdiff(params.alldofs, params.fixeddofs);

[params.H, params.Hs] = RadiusFilter(params);
inverterParams = params;
inverterParams.rmin = 1;
[params.HM, params.HsM] = RadiusFilter(inverterParams);

[result.obj, ~] = Objective(result.x, params, StiffnessMat);

fprintf('Objective: %.4f\n', result.obj);

for i = 1:length(result.timeFemStage)
    fprintf('    Stage %3d: Time of Fem: %.4fs, Time of Opt: %.4fs', i, result.timeFemStage(i),...
        result.timeOptStage(i));
    if mod(i, 2) == 0 || i == length(result.timeFemStage)
        fprintf('\n');
    else
        fprintf(', ');
    end
end

% visualize result
x = result.x;
Visualize(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '.png']);