function [xResult, weightResult] = ClusterSub(x, weight, params)
    M = size(x, 2);
    xResult = zeros(params.nely * params.nelx * params.NumMaterial, M);
    weightResult = zeros(params.nely * params.nelx * params.NumMaterial, M);

    for m = 1:M
        if params.dim == 2
            if params.NumMaterial == 1
                xSub = zeros(params.nely, params.nelx);
                weightSub = zeros(params.nely, params.nelx);
                
                xReshaped = reshape(x(:, m), params.nely*params.ClusterScale, params.nelx*params.ClusterScale);
                weightReshaped = reshape(weight(:, m), params.nely*params.ClusterScale, params.nelx*params.ClusterScale);

                for j = 1:params.ClusterScale
                    for i = 1:params.ClusterScale
                        xSub = xSub + xReshaped(i:params.ClusterScale:params.nely*params.ClusterScale, j:params.ClusterScale:params.nelx*params.ClusterScale);
                        weightSub = weightSub + weightReshaped(i:params.ClusterScale:params.nely*params.ClusterScale, j:params.ClusterScale:params.nelx*params.ClusterScale);
                    end
                end
                xSub = xSub / params.ClusterScale / params.ClusterScale;

                if strcmp(params.BC, 'cantilever')
                    halfNely = floor(params.nely / 2);

                    upperX = xSub(1:halfNely, :);
                    flippedLowerX = flip(xSub(halfNely+1:end, :), 1);
                    xSub(1:halfNely, :) = 0.5 * (upperX + flippedLowerX);
                    xSub(halfNely+1:end, :) = xSub(halfNely:-1:1, :);

                    upperWeight = weightSub(1:halfNely, :);
                    flippedLowerWeight = flip(weightSub(halfNely+1:end, :), 1);
                    weightSub(1:halfNely, :) = 0.5 * (upperWeight + flippedLowerWeight);
                    weightSub(halfNely+1:end, :) = weightSub(halfNely:-1:1, :);
                end

                xResult(:, m) = xSub(:);
                weightResult(:, m) = weightSub(:);
            else
                xSub = zeros(params.nely, params.nelx, params.NumMaterial);
                weightSub = zeros(params.nely, params.nelx, params.NumMaterial);

                xReshaped = reshape(x(:, m), params.nely*params.ClusterScale, params.nelx*params.ClusterScale, params.NumMaterial);
                weightReshaped = reshape(weight(:, m), params.nely*params.ClusterScale, params.nelx*params.ClusterScale, params.NumMaterial);

                for n = 1:params.NumMaterial
                    for j = 1:params.ClusterScale
                        for i = 1:params.ClusterScale
                            xSub(:, :, n) = xSub(:, :, n) + xReshaped(i:params.ClusterScale:params.nely*params.ClusterScale, j:params.ClusterScale:params.nelx*params.ClusterScale, n);
                            weightSub(:, :, n) = weightSub(:, :, n) + weightReshaped(i:params.ClusterScale:params.nely*params.ClusterScale, j:params.ClusterScale:params.nelx*params.ClusterScale, n);
                        end
                    end
                    xSub(:, :, n) = xSub(:, :, n) / (params.ClusterScale * params.ClusterScale);

                    if strcmp(params.BC, 'cantilever')
                        halfNely = floor(params.nely / 2);
    
                        upperX = xSub(1:halfNely, :, n);
                        flippedLowerX = flip(xSub(halfNely+1:end, :, n), 1);
                        xSub(1:halfNely, :, n) = 0.5 * (upperX + flippedLowerX);
                        xSub(halfNely+1:end, :, n) = xSub(halfNely:-1:1, :, n);
    
                        upperWeight = weightSub(1:halfNely, :, n);
                        flippedLowerWeight = flip(weightSub(halfNely+1:end, :, n), 1);
                        weightSub(1:halfNely, :, n) = 0.5 * (upperWeight + flippedLowerWeight);
                        weightSub(halfNely+1:end, :, n) = weightSub(halfNely:-1:1, :, n);
                    end
                end

                xResult(:, m) = xSub(:);
                weightResult(:, m) = weightSub(:);
            end
        else
        end
    end
end