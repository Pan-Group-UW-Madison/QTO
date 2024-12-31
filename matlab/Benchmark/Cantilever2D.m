clear;

system("clear");

addpath(genpath('../'));

params.dim = 2;

params.nelx = 200;
params.nely = 100;
params.Density = 0.3;

params.ClusterNelx = 15;
params.ClusterNely = 10;

params.rmin = 4;
params.epsilon = 5e-3;

params.MilpSolver = 'milp';
% params.MilpSolver = 'lp';
% params.MilpSolver = 'multilevel-lp';
% params.MilpSolver = 'lagrangian';
% params.MilpSolver = 'dantzig-wolfe';
params.PreStage = 1;

params.singleCut = false;
params.fixedD = false;
params.storeIterResult = false;
params.compareMilp = true;
params.useQuantum = false;

params.tolerance = 1e-9;

params.BC = 'cantilever';
params.objective = 'minimum compliance';
params.filter = 'radius';
% params.filter = 'pde';

params.NumMaterial = 2;

params.nu = 0.3;
if params.NumMaterial == 1
    params.E = 1.0;
    params.density = 1.0;
elseif params.NumMaterial == 4
    params.E = [0.43, 0.7, 0.94, 1.0];
    params.density = [0.3, 0.5, 0.8, 1.0];
else
    % params.E = [0.43, 0.7, 0.85, 0.94, 1.0];
    % params.density = [0.3, 0.5, 0.65, 0.8, 1.0];
    % params.E = [0.43, 0.85, 1.0, 0.94, 0.7];
    % params.density = [0.3, 0.65, 1.0, 0.8, 0.5];
    % params.E = [0.4, 0.7, 0.85, 0.9, 1.0];
    % params.density = [0.3, 0.5, 0.65, 0.8, 1.0];
    % params.E = [0.4, 0.85, 1.0, 0.9, 0.7];
    % params.density = [0.3, 0.65, 1.0, 0.8, 0.5];

    params.E = [0.6, 1.0];
    params.density = [0.4, 1.0];
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
    fprintf('  Stage %3d: Time of Fem: %8.4fs, Time of Opt: %8.4fs\t', i, result.timeFemStage(i),...
        result.timeOptStage(i));
    
    if mod(i, 2) == 0 || i == length(result.timeFemStage)
        fprintf('\n');
    end
end

for i = 1:params.NumMaterial
    fprintf('  Material %d: %.4f\t', i, sum(sum(result.x(:, :, i))) / (params.nelx * params.nely));
    
    if mod(i, 3) == 0 || i == params.NumMaterial
        fprintf('\n');
    end
end

% visualize result
% xTemp = result.x;
% x = zeros(params.nely, params.nelx, params.NumMaterial);
% reorder = [1, 3, 5, 4, 2];
% for i = 1:params.NumMaterial
%     x(:, :, reorder(i)) = xTemp(:, :, i);
% end

x = result.x;
if strcmp(params.MilpSolver, 'dantzig-wolfe') && ~isempty(params.useQuantum) && params.useQuantum
    params.MilpSolver = 'dantzig-wolfe-quantum';
end
VisualizeBinary(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '.png']);