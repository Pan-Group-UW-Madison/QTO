from mpi4py import MPI
import numpy as np
import ufl
from dolfinx import la
from dolfinx.fem import Function, form
from dolfinx.fem.petsc import create_matrix, assemble_matrix
from petsc4py import PETSc

from ..common import bcolors

class Optimizer:
    def __init__(self, descriptor, problem):
        self.comm = MPI.COMM_WORLD
        self.descriptor = descriptor
        self.problem = problem
        
    def summary(self):
        if self.comm.rank == 0:
            print("Optimization summary:")
            print("  Optimizer: " + bcolors.OKBLUE + self.optimizer_name + bcolors.ENDC)
            print("  Number of iterations: " + bcolors.OKBLUE + f"{self.num_iter}" + bcolors.ENDC)
            print("  Total analysis time: " + bcolors.OKBLUE + f"{self.analysis_time:.4f}" + bcolors.ENDC + " s")
            print("  Total optimization time: " + bcolors.OKBLUE + f"{self.optimization_time:.4f}" + bcolors.ENDC + " s")
            print("  Total running time: " + bcolors.OKBLUE + f"{self.running_time:.4f}" + bcolors.ENDC + " s", flush=True)

class SubOptimizer:
    def __init__(self, problem):
        self.comm = MPI.COMM_WORLD
        self.problem = problem
        
class DensityFilter():
    def __init__(self, problem, R, petsc_options={}):
        """Construct a PDE filter."""
        # Initialization
        rho = problem.rho_field
        rho_tilde = problem.rho_phys_field
        S0, S = rho.function_space, rho_tilde.function_space
        u0, u = ufl.TrialFunction(S0), ufl.TrialFunction(S)
        v, self.af = ufl.TestFunction(S), Function(S)

        self.rho, self.rho_tilde = rho, rho_tilde
        self.rho_tilde_wrap = la.create_petsc_vector_wrap(self.rho_tilde.x)
        self.af_wrap = la.create_petsc_vector_wrap(self.af.x)
        self.vec_s0, self.vec_s = rho.vector.copy(), rho_tilde.vector.copy()

        # Construct Kf and T matrices based on the Helmholtz PDE
        dx = ufl.Measure("dx", metadata={"quadrature_degree": 2})
        Kf_expr = (R**2*ufl.dot(ufl.grad(u), ufl.grad(v)) + u*v)*dx
        T_expr = u0*v*dx
        Kf_form, T_form = form(Kf_expr), form(T_expr)
        Kf_mat, self.T_mat = create_matrix(Kf_form), create_matrix(T_form)

        # Construct a filtering solver
        comm = MPI.COMM_WORLD
        self.solver = PETSc.KSP().create(comm)
        self.solver.setOperators(Kf_mat)
        prefix = f"filter_solver_{id(self)}"
        self.solver.setOptionsPrefix(prefix)

        # Apply PETSc options
        opts = PETSc.Options()
        opts.prefixPush(prefix)
        for key, value in petsc_options.items():
            opts[key] = value
        opts.prefixPop()
        self.solver.setFromOptions()
        Kf_mat.setOptionsPrefix(prefix)
        Kf_mat.setFromOptions()

        # Assemble Kf and T matrices
        assemble_matrix(Kf_mat, Kf_form)
        Kf_mat.assemble()
        assemble_matrix(self.T_mat, T_form)
        self.T_mat.assemble()
        self.T_mat_transpose = self.T_mat.copy()
        self.T_mat_transpose.transpose()

    def forward(self):
        """Compute the filtered variables."""
        self.T_mat.mult(self.rho.vector, self.vec_s)
        self.solver.solve(self.vec_s, self.rho_tilde_wrap)
        self.rho_tilde.x.scatter_forward()
        return self.rho_tilde
        pass

    def backward(self, sf_vectors):
        """Recover the sensitivities."""
        values = []
        for sf in sf_vectors:
            if sf is not None:
                self.solver.solve(sf, self.af_wrap)
                self.af.x.scatter_forward()
                self.T_mat_transpose.mult(self.af.vector, self.vec_s0)
                values.append(self.vec_s0.array.copy())
            else:
                values.append(None)
        return values

class Heaviside():
    def __init__(self, problem):
        self.rho_phys = problem.rho_phys_field

    def forward(self, beta, eta=0.5):
        denominator = np.tanh(beta*eta) + np.tanh(beta*(1-eta))
        self.drho = beta*(1-np.tanh(beta*(self.rho_phys.vector-eta))**2) / denominator
        self.rho_phys.vector.array = (
            np.tanh(beta*eta)+np.tanh(beta*(self.rho_phys.vector-eta))) / denominator
        self.rho_phys.x.scatter_forward()

    def backward(self, vectors):
        for vector in vectors:
            if vector is not None:
                vector.array *= self.drho