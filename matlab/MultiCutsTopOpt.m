function result = MultiCutsTopOpt(params)
    StiffnessMat = ElementMat(params);

    params = BCs(params);

    params.alldofs = 1:2 * (params.nely + 1) * (params.nelx + 1);
    params.freedofs = setdiff(params.alldofs, params.fixeddofs);

    if strcmp(params.filter, 'radius')
        [params.H, params.Hs] = RadiusFilter(params);
    elseif strcmp(params.filter, 'pde')
        [params.LF, params.TF] = PDEFilter(params);
    end

    if strcmp(params.BC, 'inverter')
        EminList = [1e-2, 1e-9];
        MassList = [params.Density, params.Density];
    else
        DensityMin = params.Density;
        DensityMax = params.Density0;

        N = params.N;
        A = -(N-1) / log(DensityMin / DensityMax);
        if isnan(A)
            A = 1;
        end
        MassList = [DensityMax * exp(-(0:N-1) / A) DensityMin];
        MassList = round(MassList, 3);
        EminList = 1e-2 * ones(length(MassList), 1);
        EminList(1) = 1e-2;
        EminList(end) = params.Emin;
    end

    % initialize density
    if params.NumMaterial == 1
        x = MassList(1) * ones(params.nely, params.nelx);
    else
        x = zeros(params.nely, params.nelx, params.NumMaterial);
        densitySum = sum(params.density);
        for i = 1:params.NumMaterial
            x(:, :, i) = params.density(i) / densitySum * ones(params.nely, params.nelx);
        end
    end

    result.numFem = 0;
    result.timeFem = 0;
    result.timeOpt = 0;

    params.d = 1.0;

    stage = 1;

    result.timeFemStage = zeros(length(MassList), 1);
    result.timeOptStage = zeros(length(MassList), 1);

    while stage <= length(MassList)
        fprintf("Stage %d\n", stage);

        params.Emin = EminList(stage);
        params.mass = MassList(stage);

        collection.x = [];
        collection.weight = [];
        collection.obj = [];
        collection.cost = [];

        history = [];

        % jump start
        if stage == 1 || stage == length(MassList)
            t1 = tic;
            [obj, params] = Objective(x, params, StiffnessMat);
            sensitivity = Sensitivity(x, params, StiffnessMat);
            t2 = toc(t1);

            result.numFem = result.numFem + 1;
            result.timeFem = result.timeFem + t2;
            result.timeFemStage(stage) = result.timeFemStage(stage) + t2;
        end

        collection = OptReshape(x, sensitivity, obj, 0, collection);

        t1 = tic;
        [multiCutsResult, ~] = MultiCuts(collection.x, collection.obj, collection.weight, params, history);
        t2 = toc(t1);

        result.timeOpt = result.timeOpt + t2;
        result.timeOptStage(stage) = result.timeOptStage(stage) + t2;

        x = multiCutsResult.x;
        x = reshape(x, params.nely, params.nelx, params.NumMaterial);

        cTarget = obj;

        collection.x = [];
        collection.weight = [];
        collection.obj = [];
        collection.cost = [];

        % normal start
        t1 = tic;
        [obj, params] = Objective(x, params, StiffnessMat);
        sensitivity = Sensitivity(x, params, StiffnessMat);
        t2 = toc(t1);

        result.numFem = result.numFem + 1;
        result.timeFem = result.timeFem + t2;
        result.timeFemStage(stage) = result.timeFemStage(stage) + t2;

        if stage == 1
            params.d = params.d0;
        else
            params.d = UpdateTrustRegion(cTarget, obj, multiCutsResult.obj, params.d, params);
        end
        fprintf("  Initial trust region: %.4f\n", params.d);

        upperBound = obj;
        optimalX = x;

        collection = OptReshape(x, sensitivity, obj, multiCutsResult.obj, collection);

        if params.verbose
            fprintf("  FEM: %3d, Objective: %8.4f, Upper Bound: %8.4f\n", result.numFem, obj, upperBound);
        end

        if params.visualizeStep
            Visualize(x, params, ['Step/step_' num2str(result.numFem) '.png']);
        end

        % main loop
        while (result.numFem < params.maxFem)
            t1 = tic;
            [multiCutsResult, history] = MultiCuts(collection.x, collection.obj, collection.weight, params, history);
            t2 = toc(t1);
    
            result.timeOpt = result.timeOpt + t2;
            result.timeOptStage(stage) = result.timeOptStage(stage) + t2;
    
            x = multiCutsResult.x;
            x = reshape(x, params.nely, params.nelx, params.NumMaterial);

            t1 = tic;
            [obj, params] = Objective(x, params, StiffnessMat);
            sensitivity = Sensitivity(x, params, StiffnessMat);
            t2 = toc(t1);
    
            result.numFem = result.numFem + 1;
            result.timeFem = result.timeFem + t2;
            result.timeFemStage(stage) = result.timeFemStage(stage) + t2;
    
            collection = OptReshape(x, sensitivity, obj, multiCutsResult.obj, collection);
    
            fprintf("  Lambda: %d\n", multiCutsResult.lambda);
            newD = UpdateTrustRegion(collection.obj(multiCutsResult.lambda), obj, collection.cost(multiCutsResult.lambda + 1), params.d(multiCutsResult.lambda), params);
            params.d = [params.d newD];

            % convergence of the cuts
            condition1 = abs(multiCutsResult.obj - upperBound) / abs(upperBound) < params.epsilon && abs(obj - multiCutsResult.obj) / abs(upperBound) < params.epsilon;
            condition2 = obj > upperBound && abs(obj - upperBound) / abs(upperBound) > 5 * params.epsilon && abs(obj - upperBound) / abs(upperBound) < 10 * params.epsilon;

            condition = condition1 || condition2;

            if params.verbose
                fprintf("  FEM: %3d, Objective: %8.4f, Upper Bound: %8.4f, Cost: %8.4f, difference ratio: %8.4f, %8.4f\n", result.numFem, obj, upperBound, multiCutsResult.obj, ...
                    abs(multiCutsResult.obj - upperBound) / abs(upperBound), abs(obj - multiCutsResult.obj) / abs(upperBound));
            end

            if params.visualizeStep
                Visualize(x, params, ['Step/step_' num2str(result.numFem) '.png']);
            end

            if condition
                if obj < upperBound
                    upperBound = obj;
                    optimalX = x;
                end

                x = optimalX;
                sensitivity = optimalSensitivity;
                obj = upperBound;

                break;
            end

            if obj < upperBound
                upperBound = obj;
                optimalX = x;
                optimalSensitivity = sensitivity;
            end
        end

        stage = stage + 1;

        % clear
        clear collection history;
        if stage < length(MassList)
            params.d = max(MassList(stage - 1) - MassList(stage), params.d(end));
        else
            params.d = params.d0;
        end
    end

    if params.visualizeStep
        Visualize(x, params, ['Step/step_' num2str(result.numFem+1) '.png']);
    end

    result.x = optimalX;
    result.obj = upperBound;
end

function collection = OptReshape(x, sensitivity, obj, cost, collection)
    collection.x = [collection.x, x(:)];
    collection.weight = [collection.weight, sensitivity(:)];
    collection.obj = [collection.obj, obj];
    collection.cost = [collection.cost, cost];
end