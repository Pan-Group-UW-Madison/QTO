function [result, history] = MultiCuts(x, obj, weight, params, history)
    n = size(x, 2);

    maxM = min(n, 4); % maximum cardinality of the cut
    maxR = min(maxM, n);

    if n > 1
        multiCutResult = MultiCutsSub(x(:, end), obj(end), weight(:, end), params.d(end), params);
        history = AddSolution(history, multiCutResult.x, multiCutResult.obj, n);

        % if the last cut is the best
        if multiCutResult.obj < min(history.objFuncHistory{1}(1:end-1))
            result.x = multiCutResult.x;
            result.obj = multiCutResult.obj;
            result.lambda = n;

            history = FlagSolution(history, n);

            return;
        end

        minObj = inf(maxR, 1);
        xForMinObj = ones(n, size(x, 1));
        lambdaForMinObj = zeros(maxM, maxM);

        minObj(1) = multiCutResult.obj;
        xForMinObj(1, :) = multiCutResult.x;
        lambdaForMinObj(1, 1) = n;

        for p = 2:maxM
            candidates = nchoosek(1:n, p);
            estimatedLowerBound = zeros(size(candidates, 1), 1);
            for j = 1:size(candidates, 1)
                estimatedLowerBound(j) = max(history.objFuncHistory{1}(candidates(j, :)));

                for mm = 2:p-1
                    largerCandidates = nchoosek(candidates(j, :), mm);
                    for ll = 1:size(largerCandidates, 1)
                        [isExist, ~, objFunc] = CheckSolution(history, largerCandidates(ll, :));
                        if isExist
                            estimatedLowerBound(j) = max(estimatedLowerBound(j), objFunc);
                        end
                    end
                end
            end

            minObj(p) = min(estimatedLowerBound);
        end

        minVal1 = minObj(1);
        idx1 = 1;
        minVal2 = min(minObj(2:end));

        if minVal1 < minVal2
            lambda = lambdaForMinObj(idx1, 1:idx1);

            history = FlagSolution(history, lambda);

            result.x = xForMinObj(idx1, :);
            result.obj = minVal1;
            result.lambda = lambda;

            return;
        end

        minObj(2:end) = inf;

        % branch and bound for the optimal cut
        for m = 2:maxM
            candidates = nchoosek(1:n, m);
            % remove any existing solution
            existingSolution = zeros(size(candidates));
            solutionCounter = 0;
            for j = 1:size(candidates, 1)
                useFlag = FlagCheck(history, candidates(j, :));
                if useFlag
                    solutionCounter = solutionCounter + 1;
                    existingSolution(solutionCounter, :) = candidates(j, :);
                end
            end
            existingSolution = existingSolution(1:solutionCounter, :);
            if ~isempty(existingSolution)
                candidates = setdiff(candidates, existingSolution, 'rows');
            end
            estimatedLowerBound = zeros(size(candidates, 1), 1);
            for j = 1:size(candidates, 1)
                estimatedLowerBound(j) = max(history.objFuncHistory{1}(candidates(j, :)));

                for mm = 2:m-1
                    largerCandidates = nchoosek(candidates(j, :), mm);
                    for ll = 1:size(largerCandidates, 1)
                        [isExist, ~, objFunc] = CheckSolution(history, largerCandidates(ll, :));
                        if isExist
                            estimatedLowerBound(j) = max(estimatedLowerBound(j), objFunc);
                        end
                    end
                end
            end
            [estimatedLowerBound, idx] = sort(estimatedLowerBound);
            candidates = candidates(idx, :);

            objCandidates = inf(size(candidates, 1), 1);

            for j = 1:size(candidates, 1)
                [isExist, xCandidate, objFunc] = CheckSolution(history, candidates(j, :));
                if isExist
                    objCandidates(j) = objFunc;
                else
                    multiCutResult = MultiCutsSub(x(:, candidates(j, :)), obj(candidates(j, :)), weight(:, candidates(j, :)), params.d(candidates(j, :)), params);
                    xCandidate = multiCutResult.x;
                    objCandidates(j) = multiCutResult.obj;
                    history = AddSolution(history, multiCutResult.x, multiCutResult.obj, candidates(j, :));
                end
                if j == 1 || objCandidates(j) < min(objCandidates(1:j-1))
                    xForMinObj(m, :) = xCandidate;
                    minObj(m) = objCandidates(j);
                end

                if j == size(candidates, 1) || min(objCandidates(1:j)) < min(estimatedLowerBound(j+1:end)) || min(minObj(1:m-1)) < min(estimatedLowerBound(j+1:end))
                    [val, idx] = min(objCandidates(1:j));
                    minObj(m) = val;
                    lambdaForMinObj(m, 1:m) = candidates(idx, :);
                    break;
                end
            end

            for p = m+1:maxM
                candidates = nchoosek(1:n, p);
                estimatedLowerBound = zeros(size(candidates, 1), 1);
                for j = 1:size(candidates, 1)
                    estimatedLowerBound(j) = max(history.objFuncHistory{1}(candidates(j, :)));

                    for mm = 2:p-1
                        largerCandidates = nchoosek(candidates(j, :), mm);
                        for ll = 1:size(largerCandidates, 1)
                            [isExist, ~, objFunc] = CheckSolution(history, largerCandidates(ll, :));
                            if isExist
                                estimatedLowerBound(j) = max(estimatedLowerBound(j), objFunc);
                            end
                        end
                    end
                end

                minObj(p) = min(estimatedLowerBound);
            end

            [minVal1, idx1] = min(minObj(1:m));
            if m ~= maxR
                minVal2 = min(minObj(m+1:end));
            else
                minVal2 = inf;
            end

            if minVal1 < minVal2
                lambda = lambdaForMinObj(idx1, 1:idx1);

                history = FlagSolution(history, lambda);

                result.x = xForMinObj(idx1, :);
                result.obj = minVal1;
                result.lambda = lambda;

                break;
            end
        end
    else
        multiCutResult = MultiCutsSub(x, obj, weight, params.d, params);
        history = AddSolution(history, multiCutResult.x, multiCutResult.obj, 1);

        result.x = multiCutResult.x;
        result.obj = multiCutResult.obj;
        result.lambda = 1;
    end
end

function history = AddSolution(history, x, objFunc, lambda)
    if isempty(history)
        history.solutionHistory = cell(6, 1);
        history.xHistory = cell(6, 1);
        history.objFuncHistory = cell(6, 1);
        history.useFlagHistory = cell(6, 1);
    end

    n = size(lambda, 2);
    history.solutionHistory{n} = [history.solutionHistory{n}; lambda];
    history.xHistory{n} = [history.xHistory{n}; x'];
    history.objFuncHistory{n} = [history.objFuncHistory{n}; objFunc];
    history.useFlagHistory{n} = [history.useFlagHistory{n}; false];
end

function [isExist, x, objFunc] = CheckSolution(history, lambda)
    n = length(lambda);
    partialSolution = history.solutionHistory{n};

    x = [];
    objFunc = 1;

    isExist = false;
    for i = 1:size(partialSolution, 1)
        if norm(partialSolution(i, :) - lambda) < 1e-6
            isExist = true;
            x = history.xHistory{n}(i, :);
            objFunc = history.objFuncHistory{n}(i);
            return;
        end
    end
end

function history = FlagSolution(history, lambda)
    n = length(lambda);
    partialSolution = history.solutionHistory{n};

    for i = 1:size(partialSolution, 1)
        if norm(partialSolution(i, :) - lambda) < 1e-6
            history.useFlagHistory{n}(i) = true;
            return;
        end
    end
end

function isUsed = FlagCheck(history, lambda)
    n = length(lambda);
    partialSolution = history.solutionHistory{n};

    isUsed = false;
    for i = 1:size(partialSolution, 1)
        if norm(partialSolution(i, :) - lambda) < 1e-6
            isUsed = history.useFlagHistory{n}(i);
            return;
        end
    end
end