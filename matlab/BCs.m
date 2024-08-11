function params = BCs(params)
    if params.dim == 2
        if strcmp(params.BC, 'cantilever')
            params.fixeddofs = 1:2 * (params.nely + 1);
            params.F = sparse(2*(params.nelx+1)*(params.nely+1)-2*floor(params.nely/2), 1, -1, 2 * (params.nelx + 1) * (params.nely + 1), 1);
            params.U = zeros(2 * (params.nelx + 1) * (params.nely + 1), 1);
        elseif strcmp(params.BC, 'mbb')
            params.fixeddofs = union(1:2:2 * (params.nely + 1), 2 * (params.nelx + 1) * (params.nely + 1));
            params.F = sparse(2, 1, -1, 2 * (params.nely + 1) * (params.nelx + 1), 1);
            params.U = zeros(2 * (params.nelx + 1) * (params.nely + 1), 1);
        elseif strcmp(params.BC, 'bridge')
            params.fixeddofs = union(1:2:2 * (params.nely + 1), 2 * (params.nelx + 1) * (params.nely + 1));
            params.F = sparse(2 * (params.nely + 1), 1, -0.5, 2 * (params.nely + 1) * (params.nelx + 1), 1);
            params.U = zeros(2 * (params.nelx + 1) * (params.nely + 1), 1);
        elseif strcmp(params.BC, 'two_point_brige')
            params.fixeddofs = union(1:2:2 * (params.nely + 1), 2 * (params.nelx + 1) * (params.nely + 1));
            params.F = sparse([2 * (params.nely + 1); 2 * (params.nely + 1) * (params.nelx / 2)], [1; 1], [-1; -1], 2 * (params.nely + 1) * (params.nelx + 1), 1);
            params.U = zeros(2 * (params.nelx + 1) * (params.nely + 1), 1);
        elseif strcmp(params.BC, 'inverter')
            params.fixeddofs = union(2:2*(params.nely+1):2*(params.nelx+1)*(params.nely+1),2*(params.nely+1):-1:2*(params.nely+1)-1);
            params.F = zeros(2*(params.nely+1)*(params.nelx+1), 2);
            params.din = 1;
            params.dout = 2*params.nelx*(params.nely+1)+1;
            params.F(params.din, 1) = 1;
            params.F(params.dout, 2) = -1;
            params.U = zeros(2 * (params.nelx + 1) * (params.nely + 1), 2);
        end
    end
end