function [LF, TF] = PDEFilter(params)
    nodenrs = reshape(1:(1+params.nelx)*(1+params.nely),1+params.nely,1+params.nelx);

    Rmin = params.rmin/2/sqrt(3);
    
    KEF = Rmin^2 * [4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4] / 6 + ...
                   [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4] / 36;

    edofVecF = reshape(nodenrs(1:end-1,1:end-1), params.nelx*params.nely, 1);
    edofMatF = repmat(edofVecF,1,4)+repmat([0 params.nely+(1:2) 1], params.nelx*params.nely, 1);

    iKF = reshape(kron(edofMatF,ones(4,1))', 16*params.nelx*params.nely, 1);
    jKF = reshape(kron(edofMatF,ones(1,4))', 16*params.nelx*params.nely, 1);
    sKF = reshape(KEF(:)*ones(1,params.nelx*params.nely), 16*params.nelx*params.nely, 1);
    KF = sparse(iKF, jKF, sKF);

    LF = chol(KF, 'lower');

    iTF = reshape(edofMatF, 4*params.nelx*params.nely, 1);
    jTF = reshape(repmat(1:params.nelx*params.nely,4,1)', 4*params.nelx*params.nely, 1);
    sTF = repmat(1/4, 4*params.nelx*params.nely, 1);
    TF = sparse(iTF, jTF, sTF);
end