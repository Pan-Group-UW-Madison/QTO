function result = MultiCutsSubLagrangian(x, obj, weight, d, params)
    n = size(x, 2);
    if n == 1
        numCon = 2;
    else
        numCon = 2 * n + 1;
    end

    l0 = size(x, 1) / params.NumMaterial;

    oldLagrangian = zeros(numCon, 1);

    difference = 0.5;

    clusterSize = floor(params.nelx * params.nely / 400);

    freeVar = 1:l0;

    x0 = x(:, 1);

    loop = 0;
    while difference > 1e-3 || clusterSize > 1
        loop = loop + 1;

        [xSub, weightSub] = Cluster(x, weight, params, clusterSize);

        lagrangian = DualSub(xSub, obj, weightSub, d, params, freeVar);

        if loop == 1
            oldLagrangian = 0.5 * lagrangian;
        end

        result1 = PrimalSub(x, x0, obj, weight, d, params, lagrangian, freeVar);
        result2 = PrimalSub(x, x0, obj, weight, d, params, oldLagrangian, freeVar);

        x1 = result1.x;
        x2 = result2.x;

        x1 = reshape(x1, params.nely, params.nelx, params.NumMaterial);
        x2 = reshape(x2, params.nely, params.nelx, params.NumMaterial);

        diff = zeros(params.nely, params.nelx);

        for i = 1:params.NumMaterial
            diff = diff + abs(x1(:, :, i) - x2(:, :, i));
        end

        % diff = min(diff(:), 1);

        % freeVar = find(diff > 1e-3);

        % fprintf(" difference: %5d\n", length(freeVar));

        difference = norm(lagrangian - oldLagrangian) / norm(lagrangian);

        oldLagrangian = 0.5 * (lagrangian + oldLagrangian);

        clusterSize = max(floor(clusterSize / 2), 1);
    end

    result = PrimalSub(x, x0, obj, weight, d, params, lagrangian, freeVar);
end

function [xSub, weightSub] = Cluster(x, weight, params, clusterRatio)
    if clusterRatio == 1
        xSub = x;
        weightSub = weight;
        return;
    end

    M = size(x, 2);

    l = size(x, 1) / params.NumMaterial;

    fullClusteredSize = floor(l / clusterRatio);
    outClusteredSize = l - fullClusteredSize * clusterRatio;

    xSub = zeros((fullClusteredSize + outClusteredSize) * params.NumMaterial, M);
    weightSub = zeros((fullClusteredSize + outClusteredSize) * params.NumMaterial, M);

    for m = 1:M
        xReshaped = reshape(x(:, m), params.nely, params.nelx, params.NumMaterial);
        weightReshaped = reshape(weight(:, m), params.nely, params.nelx, params.NumMaterial);

        xM = zeros(fullClusteredSize + outClusteredSize, params.NumMaterial);
        weightM = zeros(fullClusteredSize + outClusteredSize, params.NumMaterial);
        for n = 1:params.NumMaterial
            xFlatten = xReshaped(:, :, n);
            xFlatten = xFlatten(:);
            weightFlatten = weightReshaped(:, :, n);
            weightFlatten = weightFlatten(:);

            xSubReshaped = reshape(xFlatten(1:fullClusteredSize*clusterRatio), clusterRatio, fullClusteredSize);
            xM(1:fullClusteredSize, n) = sum(xSubReshaped) / clusterRatio;
            xM(fullClusteredSize+1:end, n) = sum(xFlatten(clusterRatio*fullClusteredSize+1:end)) / outClusteredSize;

            weightSubReshaped = reshape(weightFlatten(1:fullClusteredSize*clusterRatio), clusterRatio, fullClusteredSize);
            weightM(1:fullClusteredSize, n) = sum(weightSubReshaped);
            weightM(fullClusteredSize+1:end, n) = sum(weightFlatten(clusterRatio*fullClusteredSize+1:end));
        end

        xSub(:, m) = xM(:);
        weightSub(:, m) = weightM(:);
    end
end

function result = PrimalSub(x, y, obj, weight, d, params, lagrangian, freeVar)
    n = size(x, 2);

    l0 = size(x, 1) / params.NumMaterial;
    fixedVar = setdiff(1:l0, freeVar);
    
    numEle = l0;

    if params.NumMaterial == 1
        idx1 = freeVar;
        idx2 = fixedVar;

        idx1Compressed = idx1;
        idx2Compressed = idx2;
    else
        idx1Compressed = freeVar;
        idx2Compressed = fixedVar;

        idx1 = zeros(length(idx1Compressed) * params.NumMaterial, 1);
        idx2 = zeros(length(idx2Compressed) * params.NumMaterial, 1);

        for j = 1:params.NumMaterial
            idx1((j-1)*length(idx1Compressed)+1:j*length(idx1Compressed)) = (j-1)*numEle + idx1Compressed;
            idx2((j-1)*length(idx2Compressed)+1:j*length(idx2Compressed)) = (j-1)*numEle + idx2Compressed;
        end
    end

    if n > 1
        l = numEle;

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
            val((ii-1)*numVar+1:(ii-1)*numVar+numDesignVar) = -weight(idx1, ii);
            val(ii*numVar) = -1;

            rhs(ii) = -obj(ii) - weight(:, ii)' * x(:, ii) + weight(idx2, ii)' * y(idx2);
        end

        % trust region
        for ii = 1:numTrustRegion
            row(numSensitivity*numVar+(ii-1)*numDesignVar+1:numSensitivity*numVar+ii*numDesignVar) = numSensitivity+ii;
            col(numSensitivity*numVar+(ii-1)*numDesignVar+1:numSensitivity*numVar+ii*numDesignVar) = 1:numDesignVar;

            moveLimit = zeros(numEle, 1);
            for jj = 1:params.NumMaterial
                moveLimit = moveLimit + x((jj-1)*numEle+1:jj*numEle, ii);
            end
            compressedX0 = moveLimit;
            moveLimit = 1 - 2 * moveLimit;
            moveLimit = moveLimit';

            d0 = d(ii) * numEle - sum(compressedX0.^2);

            offset = numSensitivity*numVar + (ii-1)*numDesignVar;
            for jj = 1:params.NumMaterial
                val(offset+(jj-1)*l+1:offset+jj*l) = moveLimit(idx1Compressed);

                d0 = d0 - moveLimit(idx2Compressed) * y((jj-1)*numEle+idx2Compressed);
            end

            rhs(numSensitivity+ii) = d0;
        end

        % volume/mass constraint
        offset = numSensitivity*numVar + numTrustRegion*numDesignVar;
        row(offset+1:offset+numDesignVar) = ones(1, numDesignVar) + numSensitivity + numTrustRegion;
        for ii = 1:params.NumMaterial
            col(offset+(ii-1)*l+1:offset+ii*l) = (ii-1)*l+(1:l);
            val(offset+(ii-1)*l+1:offset+ii*l) = params.density(ii) * ones(1, l) / numEle;

            rhs(numSensitivity+numTrustRegion+1) = rhs(numSensitivity+numTrustRegion+1) - params.density(ii) * sum(y((ii-1)*numEle+idx2Compressed)) / numEle;
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

        result.x = y;
        result.x(idx1) = gurobiResult.x(1:numDesignVar);
        result.obj = gurobiResult.x(etaOffset);
    else
        l = length(idx1) / params.NumMaterial;

        numVar = l * params.NumMaterial;

        model.obj = -weight(idx1);
        moveLimit = zeros(numEle, 1);
        for jj = 1:params.NumMaterial
            moveLimit = moveLimit + x((jj-1)*numEle+1:jj*numEle, 1);
        end
        moveLimit = 1 - 2 * moveLimit;
        moveLimit = moveLimit';

        nnz = 2*params.NumMaterial*l;

        row = zeros(nnz, 1);
        col = zeros(nnz, 1);
        val = zeros(nnz, 1);

        % volume/mass constraint
        row(1:l*params.NumMaterial) = ones(1, l*params.NumMaterial);
        col(1:l*params.NumMaterial) = 1:l*params.NumMaterial;
        for ii = 1:params.NumMaterial
            val((ii-1)*l+1:ii*l) = params.density(ii) * ones(1, l) / numEle;
        end

        % trust region
        for ii = 1:params.NumMaterial
            row(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = 2*ones(1, l);
            col(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = (ii-1)*l+1:ii*l;
            val(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = moveLimit(idx1Compressed);
        end

        A1 = sparse(row, col, val);
        model.obj = model.obj' + lagrangian' * A1;

        % material usage
        if params.NumMaterial > 1
            row = zeros(params.NumMaterial*l, 1);
            col = zeros(params.NumMaterial*l, 1);
            val = zeros(params.NumMaterial*l, 1);

            for jj = 1:params.NumMaterial
                row(jj:params.NumMaterial:params.NumMaterial*l) = (1:l);
                col(jj:params.NumMaterial:params.NumMaterial*l) = (jj-1)*l+(1:l);
                val(jj:params.NumMaterial:params.NumMaterial*l) = 1;
            end

            numCon = l;

            rhs = ones(l, 1);
        else
            row = [];
            col = [];
            val = [];

            numCon = 0;

            rhs = [];
        end

        model.A = sparse(row, col, val, numCon, numVar);
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

        result.x = y;
        result.x(idx1) = gurobiResult.x;
        result.obj = obj(1) - weight' * (result.x - x);
    end
end

function lagrangian = DualSub(x, obj, weight, d, params, freeVar)
    n = size(x, 2);

    if n == 1
        lagrangian = ones(2, 1);
    else
        lagrangian = ones(2 * n + 1, 1);
    end
    
    l0 = size(x, 1) / params.NumMaterial;
    numEle = l0;

    if n > 1
        l = floor(size(x, 1) / params.NumMaterial);
        numVar = l * params.NumMaterial + 1;
        numDesignVar = l * params.NumMaterial;
        if params.NumMaterial == 1
            numAdditionalCon = 2 * n;
        else
            numAdditionalCon = 2 * n + numEle;
        end

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
            rhs = ones(numSensitivity + numTrustRegion + 1 + l, 1);
            rhs(numSensitivity+numTrustRegion+1) = params.mass;
        end

        row = zeros(nnz, 1);
        col = zeros(nnz, 1);
        val = zeros(nnz, 1);

        % sensitivity
        for ii = 1:numSensitivity
            % sensitivity
            row((ii-1)*numVar+1:(ii-1)*numVar+numDesignVar) = ii;
            col((ii-1)*numVar+1:(ii-1)*numVar+numDesignVar) = 1:numDesignVar;
            val((ii-1)*numVar+1:(ii-1)*numVar+numDesignVar) = -weight(:, ii);
            
            row(ii*numVar) = ii;
            col(ii*numVar) = etaOffset;
            val(ii*numVar) = -1;

            rhs(ii) = -obj(ii) - weight(:, ii)' * x(:, ii);
        end

        % trust region
        for ii = 1:numTrustRegion
            moveLimit = zeros(l, 1);
            for jj = 1:params.NumMaterial
                moveLimit = moveLimit + x((jj-1)*l+1:jj*l, ii);
            end
            compressedX0 = moveLimit;
            moveLimit = 1 - 2 * moveLimit;
            moveLimit = moveLimit';

            row(numSensitivity*numVar+(ii-1)*numDesignVar+1:numSensitivity*numVar+ii*numDesignVar) = numSensitivity+ii;
            col(numSensitivity*numVar+(ii-1)*numDesignVar+1:numSensitivity*numVar+ii*numDesignVar) = 1:numDesignVar;

            offset = numSensitivity*numVar + (ii-1)*numDesignVar;
            for jj = 1:params.NumMaterial
                val(offset+(jj-1)*l+1:offset+jj*l) = moveLimit;
            end

            rhs(numSensitivity+ii) = d(ii) * l - sum(compressedX0.^2);
        end

        % volume/mass constraint
        offset = numSensitivity*numVar + numTrustRegion*numDesignVar;
        row(offset+1:offset+numDesignVar) = ones(1, numDesignVar) + numSensitivity + numTrustRegion;
        for ii = 1:params.NumMaterial
            col(offset+(ii-1)*l+1:offset+ii*l) = (ii-1)*l+1:ii*l;
            val(offset+(ii-1)*l+1:offset+ii*l) = params.density(ii) * ones(1, l) / l;
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

        A = [sparse(row, col, val)' sparse(1:numVar, 1:numVar, ones(numVar, 1))];

        model.obj = [rhs; ones(numVar, 1)];
        model.A = -A;
        model.rhs = f;
        model.sense = '<';
        model.vtype = 'C';

        model.lb = zeros(numVar+numAdditionalCon, 1);
        model.ub = inf(numVar+numAdditionalCon, 1);

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
                lagrangian = gurobiResult.x(1:2*n + 1);
                break;
            else
                fprintf("    Warning: infeasible\n");
                break;
            end
        end
    else
        l = floor(size(x, 1) / params.NumMaterial);
        numVar = l * params.NumMaterial;
        if params.NumMaterial == 1
            numAdditionalCon = 2;
        else
            numAdditionalCon = 2 + numEle;
        end

        moveLimit = zeros(l, 1);
        for ii = 1:params.NumMaterial
            moveLimit = moveLimit + x((ii-1)*l+1:ii*l, 1);
        end
        compressedX0 = moveLimit;
        moveLimit = 1 - 2 * moveLimit;
        moveLimit = moveLimit';

        d0 = d(1) * l - sum(compressedX0.^2);

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

        % mass constraint
        for ii = 1:params.NumMaterial
            row((ii-1)*l+1:ii*l) = ones(1, l);
            col((ii-1)*l+1:ii*l) = (ii-1)*l+1:ii*l;
            val((ii-1)*l+1:ii*l) = params.density(ii) * ones(1, l) / l;
        end

        % trust region
        row(params.NumMaterial*l+1:2*params.NumMaterial*l) = 2*ones(1, l*params.NumMaterial);
        col(params.NumMaterial*l+1:2*params.NumMaterial*l) = 1:numVar;
        for ii = 1:params.NumMaterial
            val(params.NumMaterial*l+(ii-1)*l+1:params.NumMaterial*l+ii*l) = moveLimit;
        end

        % material usage
        if params.NumMaterial > 1
            offset = 2*params.NumMaterial*l;
            for jj = 1:params.NumMaterial
                row(offset+jj:params.NumMaterial:offset+params.NumMaterial*l) = (1:l) + 2;
                col(offset+jj:params.NumMaterial:offset+params.NumMaterial*l) = (jj-1)*l+(1:l);
                val(offset+jj:params.NumMaterial:offset+params.NumMaterial*l) = 1;
            end
        end

        A = [sparse(row, col, val)' sparse(1:numVar, 1:numVar, ones(numVar, 1))];

        model.obj = [rhs; ones(numVar, 1)];
        model.A = -A;
        model.rhs = -weight;
        model.sense = '<';
        model.vtype = 'C';
        model.lb = zeros(numVar+numAdditionalCon, 1);
        model.ub = inf(numVar+numAdditionalCon, 1);

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
                lagrangian = gurobiResult.x(1:2);
                break;
            else
                fprintf("    Warning: infeasible\n");
                break;
            end
        end
    end
end