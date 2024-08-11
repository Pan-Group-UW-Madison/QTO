function result = MultiCutsSubLp(x, obj, weight, d, params)
    % use multilevel solution as an initial guess
    params.ClusterScale = floor(params.nelx / params.ClusterNelx);
    clusteringSizeList = linspace(1, params.ClusterScale, params.ClusterScale);
    idx = zeros(length(clusteringSizeList), 1);

    minNelx = floor(params.nelx / params.ClusterScale);
    minNely = floor(params.nely / params.ClusterScale);

    for i = 1:length(clusteringSizeList)
        nelx = minNelx * clusteringSizeList(i);
        nely = minNely * clusteringSizeList(i);
        if (mod(params.nelx, nelx) == 0) && (mod(params.nely, nely) == 0)
            idx(i) = 1;
        end
    end

    clusteringSizeList = clusteringSizeList(logical(idx));

    if length(clusteringSizeList) > 5
        clusteringSizeList = [clusteringSizeList(1:2) clusteringSizeList(end-3:end)];
    end

    y = cell(length(clusteringSizeList), 1);
    xSub = cell(length(clusteringSizeList), 1);
    weightSub = cell(length(clusteringSizeList), 1);

    paramsSub = params;
    paramsSub.MilpSolver = 'lp';
    
    % first level
    paramsSub.nelx = minNelx * clusteringSizeList(1);
    paramsSub.nely = minNely * clusteringSizeList(1);
    paramsSub.ClusterScale = floor(params.nelx / paramsSub.nelx);
    [xSub{1}, weightSub{1}] = ClusterSub(x, weight, paramsSub);

    result = MultiCutsSub(xSub{1}, obj, weightSub{1}, d, paramsSub);
    y{1} = result.x;

    if params.visualizeLevel
        if isfield(params, 'step') == false
            params.step = 1;
        end

        clf;
        Visualize(reshape(y{1}, paramsSub.nely, paramsSub.nelx, params.NumMaterial), paramsSub, ['Step/step_' num2str(params.step) '_small_original.png']);
    end

    % second level
    paramsSub.nelx = minNelx * clusteringSizeList(2);
    paramsSub.nely = minNely * clusteringSizeList(2);
    paramsSub.ClusterScale = floor(params.nelx / paramsSub.nelx);
    [xSub{2}, weightSub{2}] = ClusterSub(x, weight, paramsSub);

    result = MultiCutsSub(xSub{2}, obj, weightSub{2}, d, paramsSub);
    y{2} = result.x;

    if params.visualizeLevel
        clf;
        Visualize(reshape(y{2}, paramsSub.nely, paramsSub.nelx, params.NumMaterial), paramsSub, ['Step/step_' num2str(params.step) '_medium_original.png']);
    end

    n = size(x, 2);
    for i = 3:length(clusteringSizeList)
        targetParams = params;
        targetParams.nelx = minNelx * clusteringSizeList(i-2);
        targetParams.nely = minNely * clusteringSizeList(i-2);
        targetParams.ClusterScale = clusteringSizeList(i) / clusteringSizeList(i-2);
        y1 = Expand(y{i-2}, targetParams);
        targetParams.nelx = minNelx * clusteringSizeList(i-1);
        targetParams.nely = minNely * clusteringSizeList(i-1);
        targetParams.ClusterScale = clusteringSizeList(i) / clusteringSizeList(i-1);
        y2 = Expand(y{i-1}, targetParams);

        numEle = minNely*clusteringSizeList(i) * minNelx*clusteringSizeList(i);

        y1 = reshape(y1, minNely*clusteringSizeList(i) * minNelx*clusteringSizeList(i), params.NumMaterial);
        y2 = reshape(y2, minNely*clusteringSizeList(i) * minNelx*clusteringSizeList(i), params.NumMaterial);

        if params.visualizeLevel
            clf;
            targetParams.nelx = minNelx * clusteringSizeList(i);
            targetParams.nely = minNely * clusteringSizeList(i);
            Visualize(reshape(y1, targetParams.nely, targetParams.nelx, params.NumMaterial), targetParams, ['Step/step_' num2str(params.step) '_small.png']);
            Visualize(reshape(y2, targetParams.nely, targetParams.nelx, params.NumMaterial), targetParams, ['Step/step_' num2str(params.step) '_medium.png']);

            params.step = params.step + 1;
        end

        if params.NumMaterial == 1
            y1 = y1(:);
            y2 = y2(:);

            idx1 = find(abs(y1 - y2) > 1e-6);
            idx2 = find(abs(y1 - y2) <= 1e-6);

            idx1Compressed = idx1;
            idx2Compressed = idx2;
        else
            idx1Compressed = zeros(minNely*clusteringSizeList(i) * minNelx*clusteringSizeList(i), 1);
            for j = 1:params.NumMaterial
                idx1Compressed(abs(y1(:, j) - y2(:, j)) > 1e-6) = 1;
            end
            idx2Compressed = find(idx1Compressed == 0);
            idx1Compressed = find(idx1Compressed);

            idx1 = zeros(length(idx1Compressed) * params.NumMaterial, 1);
            idx2 = zeros(length(idx2Compressed) * params.NumMaterial, 1);

            for j = 1:params.NumMaterial
                idx1((j-1)*length(idx1Compressed)+1:j*length(idx1Compressed)) = (j-1)*numEle + idx1Compressed;
                idx2((j-1)*length(idx2Compressed)+1:j*length(idx2Compressed)) = (j-1)*numEle + idx2Compressed;
            end

            y1 = y1(:);
        end

        if n > 1
            paramsSub.nelx = minNelx * clusteringSizeList(i);
            paramsSub.nely = minNely * clusteringSizeList(i);
            paramsSub.ClusterScale = floor(params.nelx / paramsSub.nelx);
            [x0, weight0] = ClusterSub(x, weight, paramsSub);

            l = length(idx1) / params.NumMaterial;

            numDesignVar = l * params.NumMaterial;
            numVar = l * params.NumMaterial + 1;

            etaOffset = l * params.NumMaterial + 1;
            numSensitivity = n;
            numTrustRegion = n;

            f = zeros(1, numVar);
            f(etaOffset) = 1;

            if params.NumMaterial == 1
                nnz = numSensitivity*numVar + numTrustRegion*numDesignVar + numDesignVar;
                rhs = zeros(numSensitivity + numTrustRegion + 1, 1);
                rhs(numSensitivity+numTrustRegion+1) = params.mass;
            else
                nnz = numSensitivity*numVar + numTrustRegion*numDesignVar + numDesignVar*2;
                rhs = zeros(numSensitivity + numTrustRegion + l + 1, 1);
                rhs(numSensitivity+numTrustRegion+1) = params.mass;
                rhs(numSensitivity+numTrustRegion+2:numSensitivity+numTrustRegion+1+l) = 1;
            end

            row = zeros(nnz, 1);
            col = zeros(nnz, 1);
            val = zeros(nnz, 1);

            % sensitivity
            for ii = 1:numSensitivity
                % sensitivity
                row((ii-1)*numVar+1:(ii-1)*numVar+numVar) = ii;
                col((ii-1)*numVar+1:(ii-1)*numVar+numVar) = 1:numVar;
                val((ii-1)*numVar+1:(ii-1)*numVar+numDesignVar) = -weight0(idx1, ii);
                val(ii*numVar) = -1;

                rhs(ii) = -obj(ii) - weight0(:, ii)' * x0(:, ii) + weight0(idx2, ii)' * y1(idx2);
            end

            % trust region
            for ii = 1:numTrustRegion
                row(numSensitivity*numVar+(ii-1)*numDesignVar+1:numSensitivity*numVar+ii*numDesignVar) = numSensitivity+ii;
                col(numSensitivity*numVar+(ii-1)*numDesignVar+1:numSensitivity*numVar+ii*numDesignVar) = 1:numDesignVar;

                moveLimit = zeros(numEle, 1);
                for jj = 1:params.NumMaterial
                    moveLimit = moveLimit + x0((jj-1)*numEle+1:jj*numEle, ii);
                end
                compressedX0 = moveLimit;
                moveLimit = 1 - 2 * moveLimit;
                moveLimit = moveLimit';

                d0 = d(ii) * numEle - sum(compressedX0.^2);

                offset = numSensitivity*numVar + (ii-1)*numDesignVar;
                for jj = 1:params.NumMaterial
                    val(offset+(jj-1)*l+1:offset+jj*l) = moveLimit(idx1Compressed);

                    d0 = d0 - moveLimit(idx2Compressed) * y1((jj-1)*numEle+idx2Compressed);
                end

                rhs(numSensitivity+ii) = d0;
            end

            % volume/mass constraint
            offset = numSensitivity*numVar + numTrustRegion*numDesignVar;
            row(offset+1:offset+numDesignVar) = ones(1, numDesignVar) + numSensitivity + numTrustRegion;
            for ii = 1:params.NumMaterial
                col(offset+(ii-1)*l+1:offset+ii*l) = (ii-1)*l+(1:l);
                val(offset+(ii-1)*l+1:offset+ii*l) = params.density(ii) * ones(1, l) / numEle;

                rhs(numSensitivity+numTrustRegion+1) = rhs(numSensitivity+numTrustRegion+1) - params.density(ii) * sum(y1((ii-1)*numEle+idx2Compressed)) / numEle;
            end

            % material usage
            if params.NumMaterial > 1
                offset = numSensitivity*numVar + numTrustRegion*numDesignVar + numDesignVar;
                for jj = 1:params.NumMaterial
                    row(offset+jj:params.NumMaterial:offset+numDesignVar) = (1:l) + numSensitivity + numTrustRegion + 1;
                    col(offset+jj:params.NumMaterial:offset+numDesignVar) = (jj-1)*l+(1:l);
                    val(offset+jj:params.NumMaterial:offset+numDesignVar) = 1;
                end
            end

            model.obj = f;
            model.A = sparse(row, col, val);
            model.rhs = rhs;
            model.sense = '<';

            model.vtype(1:numVar) = 'C';
    
            model.lb = zeros(numVar, 1);
            model.ub = ones(numVar, 1);
            model.lb(etaOffset) = -inf;
            model.ub(etaOffset) = inf;
    
            model.modelsense = 'min';
    
            gurobiParams.outputflag = 0;
            gurobiParams.MIPGap = params.tolerance;
            gurobiParams.Feasibilitytol = params.tolerance;
            gurobiParams.OptimalityTol = params.tolerance;
            gurobiParams.PreSolve = params.PreStage;
            gurobiParams.Heuristics = params.PreStage;
            gurobiParams.TimeLimit = 300;
    
            while true
                gurobiResult = gurobi(model, gurobiParams);

                if strcmp(gurobiResult.status, 'OPTIMAL')
                    break;
                else
                    fprintf("    Warning: infeasible\n");
                    model.rhs(numSensitivity+1:numSensitivity+numTrustRegion) = model.rhs(numSensitivity+1:numSensitivity+numTrustRegion) + 0.02 * numEle;
                end
            end

            y{i} = y1;
            y{i}(idx1) = round(gurobiResult.x(1:end-1));
            eta = gurobiResult.x(end);
        else
            paramsSub.nelx = minNelx * clusteringSizeList(i);
            paramsSub.nely = minNely * clusteringSizeList(i);
            paramsSub.ClusterScale = floor(params.nelx / paramsSub.nelx);
            [x0, weight0] = ClusterSub(x, weight, paramsSub);

            l = length(idx1) / params.NumMaterial;

            numVar = l * params.NumMaterial;

            model.obj = -weight0(idx1);
            moveLimit = zeros(numEle, 1);
            for jj = 1:params.NumMaterial
                moveLimit = moveLimit + x0((jj-1)*numEle+1:jj*numEle, 1);
            end
            compressedX0 = moveLimit;
            moveLimit = 1 - 2 * moveLimit;
            moveLimit = moveLimit';

            d0 = d(1) * numEle - sum(compressedX0.^2);

            if params.NumMaterial == 1
                nnz = 2*params.NumMaterial*l;
                rhs = [params.mass; d0];
            else
                nnz = 3*params.NumMaterial*l;
                rhs = [params.mass; d0; ones(l, 1)];
            end

            row = zeros(nnz, 1);
            col = zeros(nnz, 1);
            val = zeros(nnz, 1);

            % volume/mass constraint
            row(1:l*params.NumMaterial) = ones(1, l*params.NumMaterial);
            col(1:l*params.NumMaterial) = 1:l*params.NumMaterial;
            for ii = 1:params.NumMaterial
                val((ii-1)*l+1:ii*l) = params.density(ii) * ones(1, l) / numEle;

                rhs(1) = rhs(1) - params.density(ii) * sum(y1((ii-1)*numEle+idx2Compressed, 1)) / numEle;
            end

            % trust region
            for ii = 1:params.NumMaterial
                row(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = 2*ones(1, l);
                col(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = (ii-1)*l+1:ii*l;
                val(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = moveLimit(idx1Compressed);

                d0 = d0 - moveLimit(idx2Compressed) * y1((ii-1)*numEle+idx2Compressed, 1);
            end
            rhs(2) = d0;

            % material usage
            if params.NumMaterial > 1
                offset = 2*params.NumMaterial*l;
                for jj = 1:params.NumMaterial
                    row(offset+jj:params.NumMaterial:offset+params.NumMaterial*l) = (1:l) + 2;
                    col(offset+jj:params.NumMaterial:offset+params.NumMaterial*l) = (jj-1)*l+(1:l);
                    val(offset+jj:params.NumMaterial:offset+params.NumMaterial*l) = 1;
                end
            end

            model.A = sparse(row, col, val);
            model.rhs = rhs;
            model.sense = '<';
            model.vtype = 'C';
            model.lb = zeros(numVar, 1);
            model.ub = ones(numVar, 1);
            model.modelsense = 'min';

            gurobiParams.outputflag = 0;
            gurobiParams.MIPGap = params.tolerance;
            gurobiParams.Feasibilitytol = params.tolerance;
            gurobiParams.OptimalityTol = params.tolerance;
            gurobiParams.PreSolve = params.PreStage;
            gurobiParams.Heuristics = params.PreStage;
            gurobiParams.TimeLimit = 300;

            while true
                gurobiResult = gurobi(model, gurobiParams);

                if strcmp(gurobiResult.status, 'OPTIMAL')
                    break;
                else
                    fprintf("    Warning: infeasible\n");
                    model.rhs(2) = model.rhs(2) + 0.02 * numEle;
                end
            end

            y{i} = y1;
            y{i}(idx1) = round(gurobiResult.x);
        end
    end

    if n > 1
        result.x = y{end};
        result.obj = eta;
    else
        result.x = y{end};
        result.obj = obj(1) - weight' * (result.x - x);
    end
end