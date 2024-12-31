clear;

system("clear");
system("export OMP_NUM_THREADS=40");

addpath(genpath('../'));

params.dim = 2;

params.nelx = 120;
params.nely = 40;
params.Density = 0.3;

params.ClusterNelx = 15;
params.ClusterNely = 10;

params.rmin = 2;
params.epsilon = 5e-3;

% params.MilpSolver = 'milp';
% params.MilpSolver = 'lp';
% params.MilpSolver = 'multilevel-lp';
% params.MilpSolver = 'lagrangian';
params.MilpSolver = 'dantzig-wolfe';
params.PreStage = 1;

params.tolerance = 1e-9;

params.BC = 'mbb';
params.objective = 'minimum compliance';
params.filter = 'radius';
% params.filter = 'pde';

params.singleCut = false;
params.fixedD = false;
params.storeIterResult = false;
params.compareMilp = true;

params.NumMaterial = 1;

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
    params.E = [0.55, 1.0];
    params.density = [0.5, 1.0];

    % params.E = [0.6, 1.0];
    % params.density = [0.4, 1.0];
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
    params.Emin = 1e-9;
else
    params.Density0 = 0.3;
    params.N = 1;
    params.Emin = 1e-9;
end

params.maxFem = 200;

params.d0 = 0.1;
params.verbose = true;
params.visualizeStep = false;
params.visualizeLevel = false;

result = MultiCutsTopOpt(params);

% print result
fprintf('Number of Fem: %d\n', result.numFem);
fprintf('Time of Fem: %.4fs\n', result.timeFem);
fprintf('Time of Opt: %.4fs\n', result.timeOpt);

fprintf('Objective: %.4f\n', result.obj);

for i = 1:length(result.timeFemStage)
    fprintf('    Stage %3d: Time of Fem: %8.4fs, Time of Opt: %8.4fs\t', i, result.timeFemStage(i),...
        result.timeOptStage(i));
    
    if mod(i, 2) == 0 || i == length(result.timeFemStage)
        fprintf('\n');
    end
end

for i = 1:params.NumMaterial
    fprintf('    Material %d: %.4f\n', i, sum(sum(result.x(:, :, i))) / (params.nelx * params.nely));
end

mass = 0;
for i = 1:params.NumMaterial
    mass = mass + sum(sum(result.x(:, :, i))) / (params.nelx * params.nely) * params.density(i);
end
fprintf('    Mass: %.4f\n', mass);

% visualize result
x = result.x;
if params.singleCut == false
    Visualize(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '_' num2str(params.N) '.png']);
    % VisualizeBinary(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '_' num2str(params.N) '.png']);
else
    Visualize(x, params, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '_' num2str(params.N) '_single_cut.png']);
end

if params.storeIterResult
    clf;
    hold on;

    set(gcf, 'position', [100, 200, 512, 256]);

    iteOffset = 0;
    for i = 1:length(result.objResult)
        plot((1:length(result.objResult{i}))+iteOffset, result.objResult{i}, 'r', 'LineWidth', 2);
        plot((1:length(result.costResult{i}))+iteOffset, result.costResult{i}, 'b', 'LineWidth', 2);
        plot((1:length(result.upperBoundResult{i}))+iteOffset, result.upperBoundResult{i}, 'k', 'LineWidth', 2);
        iteOffset = iteOffset + length(result.objResult{i});
    end

    xlabel('Iteration', 'Interpreter', 'latex');
    ylabel('Objective', 'Interpreter', 'latex');

    legend('Objective', 'Cost', 'Upper Bound', 'Location', 'Best', 'Interpreter', 'latex');

    set(gca, 'YScale', 'log');

    hold off;

    saveas(gcf, ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '_history.png']);
end