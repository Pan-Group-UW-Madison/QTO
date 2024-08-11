function [H, Hs] = RadiusFilter(params)
    iH = ones(params.nelx * params.nely * (2 * (ceil(params.rmin) - 1) + 1)^2, 1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;

    for i1 = 1:params.nelx
        for j1 = 1:params.nely
            e1 = (i1 - 1) * params.nely + j1;
            for i2 = max(i1 - (ceil(params.rmin) - 1), 1):min(i1 + (ceil(params.rmin) - 1), params.nelx)
                for j2 = max(j1 - (ceil(params.rmin) - 1), 1):min(j1 + (ceil(params.rmin) - 1), params.nely)
                    e2 = (i2 - 1) * params.nely + j2;
                    k = k + 1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0, params.rmin - sqrt((i1 - i2)^2 + (j1 - j2)^2));

                    if strcmp(params.BC, 'cantilever') == false
                        if i1 < params.rmin && i2 >= 2 * i1
                            sH(k) = sH(k) * 2;
                        end
                    end
                end
            end
        end
    end

    H = sparse(iH, jH, sH);
    Hs = sum(H, 2);
end