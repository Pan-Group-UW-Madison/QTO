from mpi4py import MPI
import numpy as np

from .optimizer import SubOptimizer

class MMAOptimizer(SubOptimizer):
    def __init__(self, problem):
        super().__init__(problem)
    
    def solve(self):
        pass
    
    def solve_subproblem(m, n, opt_iter, xval, xmin, xmax, xold1, xold2, df0dx, fval,
                  dfdx, low, upp, a0=1, a=None, c=None, d=None, move=0.05,
                  asyinit=0.5, asydecr=0.7, asyincr=1.2,
                  low_bnd=0.002, up_bnd=1.0, albefa=0.1, feps=1e-6):
        """Solution update scheme with the method of moving asymptotes (MMA).
        The algorithm is available in https://doi.org/10.1002/nme.1620240207.

        Minimize:
            f_0(x) + a_0*z + sum(c_i*y_i + 0.5*d_i*(y_i)^2)
        Subjected to:
            f_i(x) - a_i*z - y_i <= 0,    i = 1, 2, ..., m
            xmin_j <= x_j <= xmax_j,      j = 1, 2, ..., n
            y_i >= 0,                     i = 1, 2, ..., m
            z >= 0

        Args:
            m: The number of general constraints.
            n: The number of the variables, x_j.
            opt_iter: Iteration counter.
            xval: Current values of the variables, x_j.
            xmin, xmax: Lower and upper bounds of the variables, x_j.
            xold1, xold2: The values of x_j at one and two iterations ago.
            df0dx: The derivatives of the objective function, f_0(x),
                with respect to the variables, x_j, calculated at xval.
            fval: The values of the constraint functions, f_i(x), calculated at xval.
            dfdx: An (m, n) array with the derivatives of the constraint functions,
                f_i(x), with respect to the variables, x_j, calculated at xval.
            low, upp: Lower and upper asymptotes from the previous iteration.
            a0, a, c, d: Coefficients in the objective function.
            move: Move limit of the variables, x_j.
            asyinit: Initial rate of the asymptotes.
            asydecr: Decreasing rate of the asymptotes when the variables are oscillating.
            asyincr: Increasing rate of the asymptotes when the variables are monotonically updated.
            low_bnd, up_bnd: Lower and upper bounds for determining the asymptotes.
            albefa: A parameter for determining the bounds of the variables.
            feps: A parameter for approximating the objective function.

        Returns:
            x_new: Optimal values of the variables, x_j, in the current subproblem.
            change: Maximum change of the variables, x_j.
            low, upp: Lower and upper asymptotes calculated and used in the current subproblem.
        """

        # Initialize the coefficients of the objective function
        if a is None:
            a = np.zeros(m)
        if c is None:
            c = np.full(m, 1000)
        if d is None:
            d = np.zeros(m)
        comm = MPI.COMM_WORLD
        n_global = comm.allreduce(n, op=MPI.SUM)
        epsimin = np.sqrt(m+n_global)*1e-9

        # Evaluate the asymptotes (low and upp)
        xrange = xmax - xmin
        if opt_iter < 2.5:
            low = xval - asyinit*xrange
            upp = xval + asyinit*xrange
        else:
            idx = (xval-xold1) * (xold1-xold2)  # Check the oscillation
            gamma = np.ones(n)
            gamma[idx > 0] = asyincr  # Relax the asymptotes if monotonic
            gamma[idx < 0] = asydecr  # Tighten the asymptotes if oscillating
            low = xval - gamma*(xold1-low)
            upp = xval + gamma*(upp-xold1)
            lowmin = xval - up_bnd*xrange
            lowmax = xval - low_bnd*xrange
            uppmin = xval + low_bnd*xrange
            uppmax = xval + up_bnd*xrange
            low = np.minimum(np.maximum(low, lowmin), lowmax)
            upp = np.maximum(np.minimum(upp, uppmax), uppmin)

        # Evaluate the bounds of the variables (alpha and beta)
        alpha = np.maximum.reduce([xmin, low+albefa*(xval-low), xval-move*xrange])
        beta = np.minimum.reduce([xmax, upp-albefa*(upp-xval), xval+move*xrange])

        # Evaluate the coefficients of the approximating objective fuction (p0 and q0)
        idx = df0dx > 0
        p0, q0 = np.zeros(n), np.zeros(n)
        p0[idx], p0[~idx] = 1.001*df0dx[idx], -0.001*df0dx[~idx]
        q0[idx], q0[~idx] = 0.001*df0dx[idx], -1.001*df0dx[~idx]
        ul1 = upp - low
        p0, q0 = (p0+feps/ul1)*(upp-xval)**2, (q0+feps/ul1)*(xval-low)**2

        # Evaluate the coefficients of the approximating constraint fuctions (P_mat and Q_mat)
        P_mat, Q_mat = np.zeros((m, n)), np.zeros((m, n))
        dfdx = dfdx.reshape(m, n)
        idx = dfdx > 0
        P_mat[idx], Q_mat[~idx] = dfdx[idx], -dfdx[~idx]
        P_mat, Q_mat = P_mat*(upp-xval)**2, Q_mat*(xval-low)**2

        # Evaluate the b vector for constraint functions
        b = P_mat@(1.0/(upp-xval)) + Q_mat@(1.0/(xval-low))
        b = comm.allreduce(b, op=MPI.SUM) - fval

        # Solve the subproblem with a primal-dual interior-point approach
        x_new = solve_subproblem(m, epsimin, low, upp, alpha, beta, p0, q0,
                                P_mat, Q_mat, a0, a, b, c, d)
        change = comm.allreduce(np.max(np.abs(x_new-xval), initial=0), op=MPI.MAX)
        return x_new, change, low, upp