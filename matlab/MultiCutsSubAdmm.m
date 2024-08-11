function result = MultiCutsSubAdmm(x, obj, weight, params)
    % use multilevel solution as an initial guess
    y = InitialGuess(x, obj, weight, params);

    n = size(x, 2);
    if n > 1
    else
        l = floor(size(x, 1) / params.NumMaterial);
        numVar = size(x, 1);

        % mass constraint
        massConstraint = ones(1, size(x, 1)) / l;
        for i = 1:params.NumMaterial
            massConstraint((i-1)*l+1:i*l) = massConstraint((i-1)*l+1:i*l) * params.density(i);
        end
        [row1, col1, val1, numNewVar1, numCon1, d1] = Convert(massConstraint, params.Density, 2);
        % move limit constraint
        moveLimit = 1 - 2 * x(:, 1);
        moveLimitConstraint = moveLimit';
        moveLimitRhs = params.d(1);
        for i = 1:params.NumMaterial
            moveLimitRhs = moveLimitRhs + moveLimitConstraint((i-1)*l+1:i*l) * x((i-1)*l+1:i*l, 1);
        end
        [row2, col2, val2, numNewVar2, numCon2, d2] = Convert(moveLimitConstraint, moveLimitRhs, 2);

        A1 = sparse(row1, col1, val1, numCon1, numVar+numNewVar1);
        A2 = sparse(row2, col2, val2, numCon2, numVar+numNewVar2);

        A1EqBin = A1(1:end-1, 1:numVar);
        A1EqCon = A1(1:end-1, numVar+1:end);
        A1IneqEq = A1(end, numVar+1:end);

        A2EqBin = A2(1:end-1, 1:numVar);
        A2EqCon = A2(1:end-1, numVar+1:end);
        A2IneqEq = A2(end, numVar+1:end);

        % material usage
        for i = 1:params.NumMaterial
            row(i:params.NumMaterial:params.NumMaterial*l) = (1:l);
            col(i:params.NumMaterial:params.NumMaterial*l) = (i-1)*l+(1:l);
            val(i:params.NumMaterial:params.NumMaterial*l) = 1;
        end
        C1 = sparse(row, col, val, l, numVar);
        C2 = sparse(1:l, 1:l, ones(l, 1), l, l);

        numBinaryVar = size(x, 1) + l;
        numContinuousVar = numNewVar1 + numNewVar2;

        binaryVar = zeros(numBinaryVar, 1);
        binaryVar(1:size(x, 1)) = y(:);
        binaryVar(size(x, 1)+1:end) = ones(l, 1) - C1 * y(:);

        numBinCon1 = numCon1 - 1;
        numBinCon2 = numCon2 - 1;

        lambda = ones(numBinCon1+numBinCon2+l, 1);
        rho = 25;

        z1 = InitializeZ(A1EqBin, A1EqCon, binaryVar(1:numVar), numNewVar1);
        z2 = InitializeZ(A2EqBin, A2EqCon, binaryVar(1:numVar), numNewVar2);
        
        Q1 = [C1' * C1, C1'; C1, C2];
        Q1(1:numVar, 1:numVar) = Q1(1:numVar, 1:numVar) + (A1EqBin' * A1EqBin + A2EqBin' * A2EqBin);
        Q1 = 0.5 * rho * Q1;

        Q2 = [A1EqCon' * A1EqCon, A1EqCon' * A2EqCon; A2EqCon' * A1EqCon, A2EqCon' * A2EqCon];
        Q2 = 0.5 * rho * Q2;

        A = [A1IneqEq, sparse(1, numNewVar2); sparse(1, numNewVar1), A2IneqEq];

        for i = 1:3
            % first block
            f = [-weight(:)' zeros(1, l)];

            f(1:numVar) = f(1:numVar) + lambda(1:numBinCon1)' * A1EqBin + lambda(numBinCon1+1:numBinCon1+numBinCon2)' * A2EqBin;
            f(1:numVar) = f(1:numVar) + rho * (z1' * A1EqCon' * A1EqBin + z2' * A2EqCon' * A2EqBin);
            f(1:numVar) = f(1:numVar) + lambda(numBinCon1+numBinCon2+1:end)' * C1;
            f(1:numVar) = f(1:numVar) - rho * ones(1, l) * C1;
            f(numVar+1:end) = f(numVar+1:end) + lambda(numBinCon1+numBinCon2+1:end)' * C2;
            f(numVar+1:end) = f(numVar+1:end) - rho * ones(1, l);

            binModel.obj = f;
            binModel.Q = Q1;

            binModel.vtype(1:numVar) = 'B';
            binModel.vtype(numVar+1:numBinaryVar) = 'C';
            binModel.lb = zeros(numBinaryVar, 1);
            binModel.ub = ones(numBinaryVar, 1);
            binModel.modelsense = 'min';
            binModel.A = sparse(0, numBinaryVar);
            binModel.start = binaryVar;

            binParams.outputflag = 0;
            binParams.MIPGap = 1e-9;
            binParams.Feasibilitytol = 1e-9;
            binParams.OptimalityTol = 1e-9;
            binParams.PreSolve = 1;
            binParams.Heuristics = 1;
            binParams.TimeLimit = 100;

            binResult = gurobi(binModel, binParams);

            binaryVar = binResult.x;

            % second block
            f = zeros(1, numContinuousVar);
            f(1:numNewVar1) = lambda(1:numBinCon1)' * A1EqCon;
            f(numNewVar1+1:end) = lambda(numBinCon1+1:numBinCon1+numBinCon2)' * A2EqCon;
            f(1:numNewVar1) = f(1:numNewVar1) + rho * binaryVar(1:numVar)' * A1EqBin' * A1EqCon;
            f(numNewVar1+1:end) = f(numNewVar1+1:end) + rho * binaryVar(1:numVar)' * A2EqBin' * A2EqCon;

            conModel.obj = f;
            conModel.Q = Q2;

            conModel.vtype = 'C';
            binModel.modelsense = 'min';
            conModel.A = A;
            conModel.rhs = [d1(end); d2(end)];
            conModel.sense = '<';
            conModel.start = [z1; z2];

            conParams.outputflag = 0;
            conParams.MIPGap = 1e-9;
            conParams.Feasibilitytol = 1e-9;
            conParams.OptimalityTol = 1e-9;
            conParams.PreSolve = 1;
            conParams.Heuristics = 1;
            conParams.TimeLimit = 100;

            conResult = gurobi(conModel, conParams);

            z1 = conResult.x(1:numNewVar1);
            z2 = conResult.x(numNewVar1+1:end);

            % update lambda
            r1 = A1EqBin * binaryVar(1:numVar) + A1EqCon * z1;
            r2 = A2EqBin * binaryVar(1:numVar) + A2EqCon * z2;
            r3 = C1 * binaryVar(1:numVar) + C2 * binaryVar(numVar+1:end) - ones(l, 1);

            lambda = lambda + rho * [r1; r2; r3];
        end
    end

    result.x = binaryVar(1:numVar);
    % result.x = y(:);
    result.obj = obj(1) - weight' * (result.x - x);
end

function x = InitialGuess(x, obj, weight, params)
    paramsSub = params;
    paramsSub.nelx = params.ClusterNelx;
    paramsSub.nely = params.ClusterNely;
    paramsSub.ClusterScale = floor(params.nelx / params.ClusterNelx);
    paramsSub.MilpSolver = 'lp';

    [xSub, weightSub] = ClusterSub(x, weight, paramsSub);

    result = MultiCutsSub(xSub(:), obj, weightSub(:), paramsSub);
    x = result.x;

    x = Expand(x, paramsSub);
end

function [row, col, val, numNewVar, numCon, d] = Convert(weight, rhs, clustering)
    L = 0;
    numNewVar = 0;
    l = size(weight, 2);

    for i = 1:ceil(log(size(weight, 2)) / log(clustering))
        L = L + 1;
        l = ceil(l / 2);
        numNewVar = numNewVar + l;

        if l == clustering
            break;
        end
    end

    numVar = size(weight, 2) + numNewVar;

    offsetByLevel = zeros(L+1, 1);
    offsetByLevel(1) = 0;
    offsetByLevel(2) = size(weight, 2);
    conOffsetByLevel = zeros(L+1, 1);
    conOffsetByLevel(1) = 0;
    conOffsetByLevel(2) = ceil(size(weight, 2) / 2);
    l = size(weight, 2);

    for i = 1:L-1
        l = ceil(l / 2);
        offsetByLevel(i + 2) = offsetByLevel(i + 1) + l;
        conOffsetByLevel(i + 2) = conOffsetByLevel(i + 1) + ceil(l/2);
    end

    l1 = size(weight, 2);
    l2 = ceil(l1 / 2);
    nnz = l1 + l2;
    numCon = l2;
    for i = 1:L-1
        l1 = ceil(l1 / 2);
        l2 = ceil(l2 / 2);
        nnz = nnz + l1 + l2;
        numCon = numCon + l2;
    end
    nnz = nnz + 2;
    numCon = numCon + 1;

    row = zeros(nnz, 1);
    col = zeros(nnz, 1);
    val = zeros(nnz, 1);

    l = size(weight, 2);
    l1 = ceil(l / 2);
    row(1:3:3*l1) = 1:l1;
    row(2:3:3*l1) = 1:l1;
    row(3:3:3*l1) = 1:l1;
    col(1:3:3*l1) = 1:2:2*l1;
    col(2:3:3*l1) = 2:2:2*l1;
    col(3:3:3*l1) = offsetByLevel(2) + (1:l1);
    val(1:3:3*l1) = weight(1:2:2*l1);
    val(2:3:3*l1) = weight(2:2:2*l1);
    val(3:3:3*l1) = -1;

    offset = 3*ceil(l / 2);
    l1 = ceil(l / 2);
    for i = 2:L
        if mod(l1, 2) == 0
            l2 = ceil(l1 / 2);

            row(offset+1:3:offset+3*l2) = conOffsetByLevel(i) + (1:l2);
            row(offset+2:3:offset+3*l2) = conOffsetByLevel(i) + (1:l2);
            row(offset+3:3:offset+3*l2) = conOffsetByLevel(i) + (1:l2);
            col(offset+1:3:offset+3*l2) = offsetByLevel(i) + (1:2:2*l2);
            col(offset+2:3:offset+3*l2) = offsetByLevel(i) + (2:2:2*l2);
            col(offset+3:3:offset+3*l2) = offsetByLevel(i+1) + (1:l2);
            val(offset+1:3:offset+3*l2) = 1;
            val(offset+2:3:offset+3*l2) = 1;
            val(offset+3:3:offset+3*l2) = -1;
        else
            l2 = ceil(l1 / 2) - 1;

            row(offset+1:3:offset+3*l2) = conOffsetByLevel(i) + (1:l2);
            row(offset+2:3:offset+3*l2) = conOffsetByLevel(i) + (1:l2);
            row(offset+3:3:offset+3*l2) = conOffsetByLevel(i) + (1:l2);
            col(offset+1:3:offset+3*l2) = offsetByLevel(i) + (1:2:2*l2);
            col(offset+2:3:offset+3*l2) = offsetByLevel(i) + (2:2:2*l2);
            col(offset+3:3:offset+3*l2) = offsetByLevel(i+1) + (1:l2);
            val(offset+1:3:offset+3*l2) = 1;
            val(offset+2:3:offset+3*l2) = 1;
            val(offset+3:3:offset+3*l2) = -1;

            j1 = ceil(l1 / 2) - 1;
            j2 = ceil(l1 / 2);
            row(offset+3*j1+1:offset+3*j1+2) = conOffsetByLevel(i) + j2;
            col(offset+3*j1+1) = offsetByLevel(i) + l1;
            col(offset+3*j1+2) = offsetByLevel(i+1) + j2;
            val(offset+3*j1+1:offset+3*j1+2) = [1, -1];
        end

        offset = offset + l1 + ceil(l1 / 2);
        l1 = ceil(l1 / 2);
    end

    row(offset+1:offset+2) = numCon;
    col(offset+1:offset+2) = numVar-1:numVar;
    val(offset+1:offset+2) = [1, 1];

    d = zeros(numCon, 1);
    d(end) = rhs;
end

function z = InitializeZ(A1, A2, y, numNewVar)
    clustering = 2;

    z = zeros(numNewVar, 1);

    L = 0;
    l0 = size(y, 1);
    l = l0;

    for i = 1:ceil(log(l0) / log(clustering))
        L = L + 1;
        l = ceil(l / 2);

        if l == clustering
            break;
        end
    end

    conOffsetByLevel = zeros(L+1, 1);
    conOffsetByLevel(1) = 0;
    conOffsetByLevel(2) = ceil(l0 / 2);

    l = l0;
    for i = 1:L-1
        l = ceil(l / 2);
        conOffsetByLevel(i + 2) = conOffsetByLevel(i + 1) + ceil(l/2);
    end

    l1 = ceil(l0 / 2);
    z(1:ceil(l0 / 2)) = A1(1:l1, :) * y;
    for i=2:L
        l2 = conOffsetByLevel(i+1) - conOffsetByLevel(i);
        newRow = conOffsetByLevel(i) + (1:l2);
        oldRow = conOffsetByLevel(i-1) + (1:l1);
        z(newRow) = A2(newRow, oldRow) * z(oldRow);
        l1 = ceil(l1 / 2);
    end
end