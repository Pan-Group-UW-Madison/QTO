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
    else
        A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
        -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
        k = 1/144*A'*[1; params.nu];

        K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
            k(2) k(1) k(2) k(4) k(6) k(7);
            k(2) k(2) k(1) k(4) k(7) k(6);
            k(3) k(4) k(4) k(1) k(8) k(8);
            k(5) k(6) k(7) k(8) k(1) k(2);
            k(5) k(7) k(6) k(8) k(2) k(1)];
        K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
            k(8)  k(9)  k(12) k(5)  k(3)  k(5);
            k(10) k(10) k(13) k(7)  k(4)  k(6);
            k(6)  k(5)  k(11) k(9)  k(2)  k(10);
            k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
            k(11) k(4)  k(6)  k(12) k(10) k(13)];
        K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
            k(7)  k(6)  k(4)  k(10) k(13) k(10);
            k(5)  k(5)  k(3)  k(8)  k(12) k(9);
            k(9)  k(10) k(2)  k(6)  k(11) k(5);
            k(12) k(13) k(10) k(11) k(6)  k(4);
            k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
        K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
            k(11) k(14) k(11) k(12) k(9)  k(8);
            k(11) k(11) k(14) k(12) k(8)  k(9);
            k(13) k(12) k(12) k(14) k(7)  k(7);
            k(10) k(9)  k(8)  k(7)  k(14) k(11);
            k(10) k(8)  k(9)  k(7)  k(11) k(14)];
        K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
            k(2) k(1)  k(8)  k(4) k(6)  k(11);
            k(8) k(8)  k(1)  k(5) k(11) k(6);
            k(3) k(4)  k(5)  k(1) k(8)  k(2);
            k(5) k(6)  k(11) k(8) k(1)  k(8);
            k(4) k(11) k(6)  k(2) k(8)  k(1)];
        K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
            k(11) k(14) k(7)  k(12) k(9)  k(2);
            k(7)  k(7)  k(14) k(10) k(2)  k(9);
            k(13) k(12) k(10) k(14) k(7)  k(11);
            k(10) k(9)  k(2)  k(7)  k(14) k(7);
            k(12) k(2)  k(9)  k(11) k(7)  k(14)];
        StiffnessMat.KE = 1/((params.nu+1)*(1-2*params.nu))*...
            [K1  K2  K3  K4;
            K2'  K5  K6  K3';
            K3' K6  K5' K2';
            K4  K3  K2  K1'];

        nEle = params.nelx*params.nely*params.nelz;

        nodegrd = reshape(1:(params.nely+1)*(params.nelx+1), params.nely+1, params.nelx+1);
        nodeids = reshape(nodegrd(1:end-1,1:end-1),params.nely*params.nelx,1);
        nodeidz = 0:(params.nely+1)*(params.nelx+1):(params.nelz-1)*(params.nely+1)*(params.nelx+1);
        nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
        edofVec = 3*nodeids(:)+1;
        StiffnessMat.edofMat = repmat(edofVec, 1, 24)+ ...
            repmat([0 1 2 3*params.nely + [3 4 5 0 1 2] -3 -2 -1 ...
            3*(params.nely+1)*(params.nelx+1)+[0 1 2 3*params.nely + [3 4 5 0 1 2] -3 -2 -1]],nEle,1);
        StiffnessMat.iK = reshape(kron(StiffnessMat.edofMat,ones(24,1))',24*24*nEle,1);
        StiffnessMat.jK = reshape(kron(StiffnessMat.edofMat,ones(1,24))',24*24*nEle,1);
    end
end