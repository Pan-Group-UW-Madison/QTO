function StiffnessMat = ElementMat(params)
    if params.dim == 2
        A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
        A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
        B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
        B12 = [2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
        StiffnessMat.KE = 1 / (1 - params.nu ^ 2) / 24 * ([A11 A12; A12' A11] + params.nu * [B11 B12; B12' B11]);

        nodenrs = reshape(1 : (1 + params.nelx) * (1 + params.nely), 1 + params.nely, 1 + params.nelx);
        edofVec = reshape(2 * nodenrs(1 : end - 1, 1 : end - 1) + 1, params.nelx * params.nely, 1);
        StiffnessMat.edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * params.nely + [2 3 0 1] -2 -1], params.nelx * params.nely, 1);
        StiffnessMat.iK = reshape(kron(StiffnessMat.edofMat, ones(8, 1))', 64 * params.nelx * params.nely, 1);
        StiffnessMat.jK = reshape(kron(StiffnessMat.edofMat, ones(1, 8))', 64 * params.nelx * params.nely, 1);
    end
end