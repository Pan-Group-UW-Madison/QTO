function result = MultiCutsTopOpt(params)
    if strcmp(params.BC, 'inverter') && params.useSuperResolution == true
        superParams = params;
        superParams.nelx = params.nelx * params.superResolutionRatio;
        superParams.nely = params.nely * params.superResolutionRatio;

        StiffnessMat = ElementMat(superParams);

        superParams = BCs(superParams);
        params.superU = superParams.U;
        params.superF = superParams.F;
        params.din = superParams.din;
        params.dout = superParams.dout;

        params.alldofs = 1:2 * (superParams.nely + 1) * (superParams.nelx + 1);
        params.fixeddofs = superParams.fixeddofs;
        params.freedofs = setdiff(params.alldofs, params.fixeddofs);
    else
        StiffnessMat = ElementMat(params);

        params = BCs(params);

        params.alldofs = 1:2 * (params.nely + 1) * (params.nelx + 1);
        params.freedofs = setdiff(params.alldofs, params.fixeddofs);
    end

    if strcmp(params.filter, 'radius')
        [params.H, params.Hs] = RadiusFilter(params);
        if strcmp(params.BC, 'inverter')
            inverterParams = params;
            inverterParams.rmin = params.rmin;
            if params.useSuperResolution
                inverterParams.nelx = params.nelx * params.superResolutionRatio;
                inverterParams.nely = params.nely * params.superResolutionRatio;
            end
            [params.HM, params.HsM] = RadiusFilter(inverterParams);
        end
    elseif strcmp(params.filter, 'pde')
        [params.LF, params.TF] = PDEFilter(params);
    end

    if strcmp(params.BC, 'inverter')
        % EminList = [1e-2, 3e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6];
        % MassList = params.Density * ones(1, 9);
        EminList = [1e-2, 1e-4];
        MassList = params.Density * ones(1, 2);
    else
        DensityMin = params.Density;
        DensityMax = params.Density0;

        N = params.N;
        if N == 1
            A = 1;
        else
            A = -(N-1) / log(DensityMin / DensityMax);
        end
        MassList = [DensityMax * exp(-(0:N-1) / A) DensityMin];
        MassList = round(MassList, 3);
        EminList = 1e-2 * ones(length(MassList), 1);
        EminList(end) = params.Emin;
    end

    % initialize density
    if params.NumMaterial == 1
        x = MassList(1) * ones(params.nely, params.nelx);
    else
        x = zeros(params.nely, params.nelx, params.NumMaterial);
        densitySum = sum(params.density);
        for i = 1:params.NumMaterial
            x(:, :, i) = params.density(i) / densitySum * MassList(1) * ones(params.nely, params.nelx);
        end
    end

    result.numFem = 0;
    result.timeFem = 0;
    result.timeOpt = 0;

    params.d = 1.0;

    stage = 1;

    result.timeFemStage = zeros(length(MassList), 1);
    result.timeOptStage = zeros(length(MassList), 1);

    if params.storeIterResult
        result.objResult = cell(length(MassList), 1);
        result.costResult = cell(length(MassList), 1);
        result.upperBoundResult = cell(length(MassList), 1);
    end

    while stage <= length(MassList)
        if params.verbose
            fprintf("Stage %d\n", stage);
        end

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
            sensitivity = Sensitivity(x, params, StiffnessMat, stage == 1);
            t2 = toc(t1);

            result.numFem = result.numFem + 1;
            result.timeFem = result.timeFem + t2;
            result.timeFemStage(stage) = result.timeFemStage(stage) + t2;

            if params.verbose
                fprintf("  FEM: %3d, Objective: %8.4f\n", result.numFem, obj);
            end

            collection = OptReshape(x, sensitivity, obj, 0, collection);

            t1 = tic;
            if strcmp(params.BC, 'inverter') && strcmp(params.MilpSolver, 'milp') && stage == 1
                params.MilpSolver = 'lp';
                [multiCutsResult, history] = MultiCuts(collection.x, collection.obj, collection.weight, params, history);
                params.MilpSolver = 'milp';
            else
                [multiCutsResult, history] = MultiCuts(collection.x, collection.obj, collection.weight, params, history);
            end
            t2 = toc(t1);

            result.timeOpt = result.timeOpt + t2;
            result.timeOptStage(stage) = result.timeOptStage(stage) + t2;

            x = multiCutsResult.x;
            x = reshape(x, params.nely, params.nelx, params.NumMaterial);

            t1 = tic;
            [obj, params] = Objective(x, params, StiffnessMat);
            sensitivity = Sensitivity(x, params, StiffnessMat, false);
            t2 = toc(t1);

            result.numFem = result.numFem + 1;
            result.timeFem = result.timeFem + t2;
            result.timeFemStage(stage) = result.timeFemStage(stage) + t2;

            if stage == 1
                params.d = [params.d params.d0];
                collection = OptReshape(x, sensitivity, obj, multiCutsResult.obj, collection);

                upperBound = obj;
                optimalX = x;
                optimalSensitivity = sensitivity;
                optimalD = params.d0;
            else
                collection.x = [];
                collection.weight = [];
                collection.obj = [];
                collection.cost = [];

                history = [];

                collection = OptReshape(x, sensitivity, obj, multiCutsResult.obj, collection);

                upperBound = obj;
            end
        else
            collection = OptReshape(x, sensitivity, obj, multiCutsResult.obj, collection);

            upperBound = obj;
        end

        if params.verbose
            fprintf("  FEM: %3d, Objective: %8.4f, Upper Bound: %8.4f\n", result.numFem, obj, upperBound);
        end

        if params.visualizeStep
            Visualize(x, params, ['Step/step_' num2str(result.numFem) '.png']);
        end

        stackedIteration = 0;
        innerLoop = 0;

        % main loop
        while (result.numFem < params.maxFem)
            innerLoop = innerLoop + 1;
            t1 = tic;
            if (innerLoop == 1 && stage > 1 && stage < length(MassList) && strcmp(params.MilpSolver, 'milp'))
                params.MilpSolver = 'lp';
                [multiCutsResult, history] = MultiCuts(collection.x, collection.obj, collection.weight, params, history);
                params.MilpSolver = 'milp';
            else
                [multiCutsResult, history] = MultiCuts(collection.x, collection.obj, collection.weight, params, history);
            end
            t2 = toc(t1);
    
            result.timeOpt = result.timeOpt + t2;
            result.timeOptStage(stage) = result.timeOptStage(stage) + t2;
    
            x = multiCutsResult.x;
            x = reshape(x, params.nely, params.nelx, params.NumMaterial);

            t1 = tic;
            [obj, params] = Objective(x, params, StiffnessMat);
            sensitivity = Sensitivity(x, params, StiffnessMat, false);
            t2 = toc(t1);
    
            result.numFem = result.numFem + 1;
            result.timeFem = result.timeFem + t2;
            result.timeFemStage(stage) = result.timeFemStage(stage) + t2;
    
            collection = OptReshape(x, sensitivity, obj, multiCutsResult.obj, collection);
    
            if params.verbose
                fprintf("  Lambda");
                for i = 1:length(multiCutsResult.lambda)
                    fprintf(" %d ", multiCutsResult.lambda(i));
                end
                fprintf("\n");
            end
            newD = UpdateTrustRegion(collection.obj(multiCutsResult.lambda), obj, multiCutsResult.obj, params.d(multiCutsResult.lambda), params);
            params.d = [params.d newD];

            % convergence of the cuts
            condition1 = abs(obj - upperBound) / abs(upperBound) < params.epsilon;
            condition2 = abs(obj - multiCutsResult.obj) / abs(upperBound) < params.epsilon;
            condition3 = obj > upperBound && multiCutsResult.obj > upperBound && stage > 1 && abs(obj - upperBound) / abs(upperBound) < 1e-2 && abs(obj - multiCutsResult.obj) / abs(upperBound) < 1e-2;

            if strcmp(params.objective, 'minimum compliance')
                condition = (condition1 && condition2) || condition3;
            else
                condition = condition1;
            end

            if params.verbose
                fprintf("  FEM: %3d, Objective: %8.4f, Upper Bound: %8.4f, Cost: %8.4f, difference ratio: %8.4f, %8.4f\n", result.numFem, obj, upperBound, multiCutsResult.obj, ...
                    abs(obj - upperBound) / abs(upperBound), abs(obj - multiCutsResult.obj) / abs(upperBound));
            end

            if params.visualizeStep
                Visualize(x, params, ['Step/step_' num2str(result.numFem) '.png']);
            end

            if multiCutsResult.obj > upperBound
                stackedIteration = stackedIteration + 1;
            else
                stackedIteration = 0;
            end

            % if condition || (stackedIteration >= 3 && strcmp(params.BC, 'inverter') && stage > 1)
            if condition || (stackedIteration >= 5)
                if obj < upperBound
                    upperBound = obj;
                    optimalX = x;
                    optimalSensitivity = sensitivity;
                    optimalD = newD;
                end

                x = optimalX;
                sensitivity = optimalSensitivity;
                obj = upperBound;

                break;
            end

            if params.storeIterResult
                result.objResult{stage} = [result.objResult{stage} obj];
                result.costResult{stage} = [result.costResult{stage} multiCutsResult.obj];
                result.upperBoundResult{stage} = [result.upperBoundResult{stage} upperBound];
            end

            if obj < upperBound
                stackedIteration = 0;
                upperBound = obj;
                optimalX = x;
                optimalSensitivity = sensitivity;
                optimalD = newD;
            end
        end

        stage = stage + 1;

        % clear
        clear collection history;
        if stage < length(MassList) - 1
            params.d = max(MassList(stage - 1) - MassList(stage) + 1e-3, optimalD);
        else
            params.d = optimalD;
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