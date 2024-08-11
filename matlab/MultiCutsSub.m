

function result = MultiCutsSub(x, obj, weight, d, params)
    if strcmp(params.MilpSolver, 'multilevel-admm')
        result = MultiCutsSubAdmm(x, obj, weight, d, params);

        return;
    end

    if strcmp(params.MilpSolver, 'multilevel-lp')
        result = MultiCutsSubLp(x, obj, weight, d, params);

        return;
    end

    n = size(x, 2);
    if n > 1
        l = floor(size(x, 1) / params.NumMaterial);
        numVar = l * params.NumMaterial + 1;
        numDesignVar = l * params.NumMaterial;

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
        model.A = sparse(row, col, val);
        model.rhs = rhs;
        model.sense = '<';
        if strcmp(params.MilpSolver, 'milp')
            model.vtype(1:numDesignVar) = 'B';
            model.vtype(etaOffset) = 'C';
        else
            model.vtype(1:numVar) = 'C';
        end

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
                model.rhs(numSensitivity+1:numSensitivity+numTrustRegion) = model.rhs(numSensitivity+1:numSensitivity+numTrustRegion) + 0.02 * l;
            end
        end

        result.x = gurobiResult.x(1:numDesignVar);
        result.obj = gurobiResult.x(etaOffset);
    else
        l = floor(size(x, 1) / params.NumMaterial);
        numVar = l * params.NumMaterial;

        model.obj = -weight;
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

        model.A = sparse(row, col, val);
        model.rhs = rhs;
        model.sense = '<';
        if strcmp(params.MilpSolver, 'milp')
            model.vtype = 'B';
        else
            model.vtype = 'C';
        end
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
                model.rhs(2) = model.rhs(2) + 0.02 * l;
            end
        end

        result.x = gurobiResult.x;
        result.obj = obj(1) - weight' * (result.x - x);
    end
end