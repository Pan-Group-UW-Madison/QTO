system("clear");

params.dim = 2;

nelxList = [120, 240, 360, 480, 600];
nelyList = [80, 160, 240, 320, 400];

params.Density = 0.6;

params.rmin = 2;
params.epsilon = 5e-3;

params.MilpSolver = 'milp';
params.PreStage = 0;

params.tolerance = 1e-6;

params.BC = 'cantilever';
params.objective = 'minimum compliance';

params.NumMaterial = 1;

params.nu = 0.3;
params.E = 1.0;
params.density = 1.0;

params.xSymmetric = false;
params.ySymmetric = false;
params.mass = params.Density;

addpath(genpath('../'));

timeFem = zeros(length(nelxList), 1);
timeOpt = zeros(length(nelxList), 1);

for i = 1:length(nelxList)
    params.nelx = nelxList(i);
    params.nely = nelyList(i);

    params.d = params.nelx * params.nely;

    params.Emin = 1e-2;

    x = params.Density * ones(params.nely, params.nelx);

    StiffnessMat = ElementMat(params);

    params = BCs(params);

    params.alldofs = 1:2 * (params.nely + 1) * (params.nelx + 1);
    params.freedofs = setdiff(params.alldofs, params.fixeddofs);

    [params.H, params.Hs] = RadiusFilter(params);
    
    t1 = tic;
    [obj, params] = Objective(x, params, StiffnessMat);
    sensitivity = Sensitivity(x, params, StiffnessMat);
    timeFem(i) = toc(t1);

    t1 = tic;
    for j = 1:10
        [multiCutsResult, history] = MultiCuts(x(:), obj, sensitivity(:), params, []);
    end
    timeOpt(i) = toc(t1) / 10;

    % visualization
    x = reshape(multiCutsResult.x, params.nely, params.nelx);
    set(gcf, 'position', [100, 200, 240, floor(240*params.nely/params.nelx)]);
    set(gca, 'Position', [0, 0, 1, 1])
    myColorMap = jet(256);
    myColorMap(end,:) = 1;
    colormap(myColorMap);
    imagesc(1 - x);
    clim([0 1]);
    axis equal;
    axis off;
    drawnow;
    filename = ['Result/' num2str(params.nelx) 'x' num2str(params.nely) '_' params.MilpSolver '_' params.BC '_' num2str(params.NumMaterial) '.png'];
    saveas(gcf, filename);
end

fprintf('Single Material Unit Test of Minimum Compliance using MILP Solver\n');
for i = 1:length(nelxList)
    fprintf('    nelx: %4d, nely: %4d, timeFem: %.4fs, timeOpt: %.4fs\n', nelxList(i), nelyList(i), timeFem(i), timeOpt(i));
end