function x = Expand(x, params)
    if params.dim == 2
        if params.NumMaterial == 1
            clusterScale = params.ClusterScale;
            xReshaped = reshape(x, params.nely, params.nelx);
            [lowY, lowX] = meshgrid(1:params.nelx, 1:params.nely);
            [highY, highX] = meshgrid(linspace(1, params.nelx, params.nelx*clusterScale), linspace(1, params.nely, params.nely*clusterScale));
            x = interp2(lowY, lowX, xReshaped, highY, highX, 'cubic');
        else
            clusterScale = params.ClusterScale;
            xReshaped = reshape(x, params.nely, params.nelx, params.NumMaterial);
            x = zeros(params.nely*clusterScale, params.nelx*clusterScale, params.NumMaterial);
            [lowY, lowX] = meshgrid(1:params.nelx, 1:params.nely);
            [highY, highX] = meshgrid(linspace(1, params.nelx, params.nelx*clusterScale), linspace(1, params.nely, params.nely*clusterScale));
            for n = 1:params.NumMaterial
                x(:, :, n) = interp2(lowY, lowX, xReshaped(:, :, n), highY, highX, 'cubic');
            end
        end
    else
    end
end