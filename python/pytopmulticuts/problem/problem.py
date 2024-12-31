import numpy as np
import ufl
from dolfinx.mesh import locate_entities_boundary, meshtags
from dolfinx.fem import VectorFunctionSpace, FunctionSpace, Function, Constant, dirichletbc, locate_dofs_topological

from petsc4py import PETSc
import dolfinx.io
from dolfinx.fem import form, Function
from dolfinx import la
from dolfinx.fem.petsc import (create_vector, create_matrix,
                               assemble_vector, assemble_matrix, set_bc)

from mpi4py import MPI

from ..common import bcolors

class Problem:
    def __init__(self, descriptor):
        if descriptor["mesh"] is None:
            if MPI.COMM_WORLD.rank == 0:
                raise ValueError("Mesh is not provided.")
            exit(1)
        else:
            self.mesh = descriptor["mesh"]
        
        if descriptor["prefix"] is None:
            self.prefix = "result/"
        else:
            self.prefix = descriptor["prefix"]
            
        if descriptor["problem_name"] is None:
            self.problem_name = "problem"
        else:
            self.problem_name = descriptor["problem_name"]
        
        self.comm = self.mesh.comm
        
        self.V = None
        self.S0 = None
        self.S = None
        
        self.objective = descriptor["objective"]
        
    def set_solver(self, petsc_options):
        self.lhs_mat = create_matrix(self.lhs_form)
        self.rhs_vec = create_vector(self.rhs_form)
        self.u_wrap = la.create_petsc_vector_wrap(self.u_field.x)
        
        self.solver = PETSc.KSP().create(self.u_field.function_space.mesh.comm)
        self.solver.setOperators(self.lhs_mat)
        prefix = f"linear_solver_{id(self)}"
        self.solver.setOptionsPrefix(prefix)
        
        # Apply PETSc options
        opts = PETSc.Options()
        opts.prefixPush(prefix)
        for key, value in petsc_options.items():
            opts[key] = value
        opts.prefixPop()
        self.solver.setFromOptions()
        for var in [self.lhs_mat, self.rhs_vec, self.l_vec]:
            if var is not None:
                var.setOptionsPrefix(prefix)
                var.setFromOptions()
        
    def solve_prime(self):
        """Solve K*x=F."""
        self.lhs_mat.zeroEntries()
        assemble_matrix(self.lhs_mat, self.lhs_form, bcs=self.bcs)
        self.lhs_mat.assemble()
        if self.spring_vec is not None:
            self.lhs_mat.setDiagonal(self.lhs_mat.getDiagonal()+self.spring_vec)
        self.solver.solve(self.rhs_vec, self.u_wrap)
        self.u_field.x.scatter_forward()
        
    def solve_adjoint(self):
        """Solve K*lambda=-L."""
        self.solver.solve(-self.l_vec, self.lam_wrap)
        self.lam.x.scatter_forward()
        
    def __del__(self):
        self.solver.destroy()
        self.lhs_mat.destroy()
        self.rhs_vec.destroy()
        self.u_wrap.destroy()
        # self.lam_wrap.destroy()
        if self.spring_vec is not None:
            self.spring_vec.destroy()
            self.l_vec.destroy()
        
class LinearElasticity(Problem):
    def __init__(self, descriptor):
        super().__init__(descriptor)
        
        if descriptor["problem_name"] is None:
            self.problem_name = "linear_elasticity"
            
        self.V = VectorFunctionSpace(self.mesh, ("CG", 1))
        self.S0 = FunctionSpace(self.mesh, ("DG", 0))
        self.u, self.v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        self.u_field = Function(self.V)
        self.rho_field = Function(self.S0)
        
        # if descriptor["interpolation"] == "continuous":
        self.S = FunctionSpace(self.mesh, ("CG", 1))
        self.rho_phys_field = Function(self.S)
        
        if isinstance(descriptor["young's modulus"], (int, float)):
            self.E_list = np.array([descriptor["young's modulus"]], dtype=np.float64)
        elif isinstance(descriptor["young's modulus"], (list, tuple)):
            self.E_list = np.array(descriptor["young's modulus"], dtype=np.float64)
        elif isinstance(descriptor["young's modulus"], np.ndarray):
            self.E_list = descriptor["young's modulus"]
        else:
            raise ValueError("Young's modulus is not in the correct format.")
            exit(1)
        if np.size(self.E_list, 0) > 1:
            if isinstance(descriptor["density"], (list, tuple)):
                self.rho_list = np.array(descriptor["density"], dtype=np.float64)
            elif isinstance(descriptor["density"], np.ndarray):
                self.rho_list = descriptor["density"]
            elif descriptor["density"] is None:
                raise ValueError("Density is not provided.")
                exit(1)
            else:
                raise ValueError("Density is not in the correct format.")
                exit(1)
        self.nu = descriptor["poisson's ratio"]
        
        E0 = self.E_list[-1]
        nu = self.nu
        
        if descriptor["interpolation"] == "discrete":
            self.interpolation = "discrete"
            self.eps = 1e-2
            E = (self.eps + (1-self.eps)*self.rho_phys_field) * E0
        else:
            self.interpolation = "continuous"
            p, eps = 3, 1e-6
            E = (eps + (1-eps)*self.rho_phys_field**p) * E0
        _lambda, mu = E*nu/(1+nu)/(1-2*nu), E/(2*(1+nu))
        
        # Kinematics
        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u):  # 3D or plane strain
            return 2*mu*epsilon(u) + _lambda*ufl.tr(epsilon(u))*ufl.Identity(len(u))
        
        self.dim = self.mesh.topology.dim
        self.disp_facets = locate_entities_boundary(self.mesh, self.dim-1, descriptor["disp_bc"])
        self.bcs = [dirichletbc(Constant(self.mesh, np.full(self.dim, 0.0)), locate_dofs_topological(self.V, self.dim-1, self.disp_facets), self.V)]
        
        tractions, facets, markers = [], [], []
        for marker, (traction, traction_bc) in enumerate(descriptor["traction_bcs"]):
            tractions.append(Constant(self.mesh, np.array(traction, dtype=float)))
            current_facets = locate_entities_boundary(self.mesh, self.dim-1, traction_bc)
            facets.extend(current_facets)
            markers.extend([marker,]*len(current_facets))
        facets = np.array(facets, dtype=np.int32)
        markers = np.array(markers, dtype=np.int32)
        _, unique_indices = np.unique(facets, return_index=True)
        facets, markers = facets[unique_indices], markers[unique_indices]
        sorted_indices = np.argsort(facets)
        facet_tags = meshtags(self.mesh, self.dim-1, facets[sorted_indices], markers[sorted_indices])
        
        metadata = {"quadrature_degree": descriptor["quadrature_degree"]}
        self.dx = ufl.Measure("dx", metadata=metadata)
        self.ds = ufl.Measure("ds", domain=self.mesh, metadata=metadata, subdomain_data=facet_tags)
        b = Constant(self.mesh, np.array(descriptor["body_force"], dtype=float))
        
        lhs = ufl.inner(sigma(self.u), epsilon(self.v))*self.dx
        rhs = ufl.dot(b, self.v)*self.dx
        for marker, t in enumerate(tractions):
            rhs += ufl.dot(t, self.v)*self.ds(marker)
        if descriptor["objective"] == "compliance":
            self.spring_vec = self.l_vec = None
        else:
            self.spring_vec, self.l_vec = create_mechanism_vectors(
                V, opt["in_spring"], opt["out_spring"])
        self.lhs_form = form(lhs)
        self.rhs_form = form(rhs)
            
        super().set_solver(descriptor["petsc_options"])
        
        assemble_vector(self.rhs_vec, self.rhs_form)
        self.rhs_vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        set_bc(self.rhs_vec, self.bcs)        
        
        # Define optimization-related variables
        self.f_int = ufl.inner(sigma(self.u_field), epsilon(self.v))*self.dx
        self.compliance = ufl.inner(sigma(self.u_field), epsilon(self.u_field))*self.dx
        if self.interpolation == "discrete":
            self.volume = self.rho_field*self.dx
        else:
            self.volume = self.rho_phys_field*self.dx
        self.total_volume = Constant(self.mesh, 1.0)*self.dx
        
    def summary(self):
        if self.comm.rank == 0:
            print(bcolors.WARNING + "Problem name: ", self.problem_name + bcolors.ENDC)
            print("  Number of ranks: " + bcolors.OKBLUE, self.comm.size, bcolors.ENDC)
            if self.mesh.topology.dim == 2:
                print("  Number of cells: " + bcolors.OKBLUE, self.mesh.topology.index_map(2).size_global, bcolors.ENDC)
                print("  Number of vertices: " + bcolors.OKBLUE, self.mesh.topology.index_map(0).size_global, bcolors.ENDC)
                print("  Number of dofs: " + bcolors.OKBLUE, 2*self.V.dofmap.index_map.size_global, bcolors.ENDC)
            elif self.mesh.topology.dim == 3:
                print("  Number of cells: " + bcolors.OKBLUE, self.mesh.topology.index_map(3).size_global, bcolors.ENDC)
                print("  Number of vertices: " + bcolors.OKBLUE, self.mesh.topology.index_map(0).size_global, bcolors.ENDC)
                print("  Number of dofs: " + bcolors.OKBLUE, 3*self.V.dofmap.index_map.size_global, bcolors.ENDC)
            print("  Number of materials: " + bcolors.OKBLUE, np.size(self.E_list, 0), bcolors.ENDC, flush=True)
    
    def save_results(self):
        xdmf = dolfinx.io.XDMFFile(self.mesh.comm, self.prefix+self.problem_name+".xdmf", "w")
        xdmf.write_mesh(self.mesh)
        if self.interpolation == "discrete":
            self.rho_field.name = "density"
            xdmf.write_function(self.rho_field)
        else:
            self.rho_phys_field.name = "density"
            xdmf.write_function(self.rho_phys_field)
        xdmf.close()