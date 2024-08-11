function Visualize(y, params, filename)
    if params.NumMaterial == 1
        x = y;
        x = reshape(x, params.nely, params.nelx);
    else
        x = zeros(params.nely, params.nelx);
        for i = 1:params.NumMaterial
            x = x + i*y(:, :, i);
        end
        x = x / params.NumMaterial;
    end
    
    % visualize result
    clf;
    if params.ySymmetric
        x = [flip(x, 1); x];
    end
    set(gcf, 'position', [100, 200, 256, 256*size(x, 1)/size(x, 2)]);
    set(gca, 'Position', [0, 0, 1, 1])
    myColorMap = jet(256);
    myColorMap(end,:) = 1;
    colormap(myColorMap);
    imagesc(1 - x);
    clim([0 1]);
    axis equal;
    axis off;
    drawnow;
    saveas(gcf, filename);
end