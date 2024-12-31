from .optimizer import Optimizer, DensityFilter, Heaviside
from .sensitivity import Sensitivity

from .mma import MMAOptimizer
from .oc import OCOptimizer

import time
import numpy as np

class SimpOptimizer(Optimizer):
    def __init__(self, descriptor, problem):
        super().__init__(descriptor, problem)
        
        self.optimizer_name = "SIMP"
        
        self.max_iter = descriptor["max_iter"]
        self.opt_tol = descriptor["opt_tol"]
        self.vol_frac = descriptor["vol_frac"]
        self.radius = descriptor["filter_radius"]
        
        self.beta_interval = descriptor["beta_interval"]
        self.beta_max = descriptor["beta_max"]
        self.move = descriptor["move"]
        
        self.verbose = 1
        
        num_consts = 1 if self.problem.objective == "compliance" else 2
        
        if descriptor["subproblem_solver"] == "mma":
            self.sub_optimizer = MMAOptimizer(problem)
            num_elems = self.problem.rho_field.vector.array.size
            rho_old1, rho_old2 = np.zeros(num_elems), np.zeros(num_elems)
            low, upp = None, None
        elif descriptor["subproblem_solver"] == "oc":
            self.sub_optimizer = OCOptimizer(problem, self.move)
        else:
            raise ValueError("Invalid subproblem_solver")
            exit(1)
            
        rho_field = self.problem.rho_field
        num_elems = rho_field.vector.array.size
        centers = rho_field.function_space.tabulate_dof_coordinates()[:num_elems].T
        solid, void = descriptor["solid_zone"](centers), descriptor["void_zone"](centers)
        rho_ini = np.full(num_elems, descriptor["vol_frac"])
        rho_ini[solid], rho_ini[void] = 0.995, 0.005
        rho_field.vector.array[:] = rho_ini
        rho_min, rho_max = np.zeros(num_elems), np.ones(num_elems)
        rho_min[solid], rho_max[void] = 0.99, 0.01
        
        self.sub_optimizer.rho_min, self.sub_optimizer.rho_max = rho_min, rho_max
            
    def solve(self):
        density_filter = DensityFilter(self.problem, self.radius)
        heaviside = Heaviside(self.problem)
        sens_problem = Sensitivity(self.problem)
        
        running_timer = time.perf_counter()
        
        self.analysis_time, self.optimization_time = 0, 0
        
        self.num_iter, beta, change = 0, 1, 2*self.opt_tol
        while self.num_iter < self.max_iter and change > self.opt_tol:
            opt_start_time = time.perf_counter()
            self.num_iter += 1
            
            density_filter.forward()
            if self.num_iter % self.beta_interval == 0 and beta < self.beta_max:
                beta *= 2
            heaviside.forward(beta)
            
            # Solve FEM
            fem_sen_time = time.perf_counter()
            self.problem.solve_prime()
            
            # Compute function values and sensitivities
            [C_value, V_value, U_value], sensitivities = sens_problem.evaluate()
            heaviside.backward(sensitivities)
            [dCdrho, dVdrho, dUdrho] = density_filter.backward(sensitivities)
            if self.problem.objective == "compliance":
                g_vec = np.array([V_value-self.vol_frac])
                dJdrho, dgdrho = dCdrho, np.vstack([dVdrho])
            else:
                g_vec = np.array([V_value-self.vol_frac, C_value-opt["compliance_bound"]])
                dJdrho, dgdrho = dUdrho, np.vstack([dVdrho, dCdrho])
            fem_sen_time = time.perf_counter() - fem_sen_time
            self.analysis_time += fem_sen_time
            
            # Update the design variables
            opt_time = time.perf_counter()
            rho_values = self.problem.rho_field.vector.array.copy()
            rho_new = self.sub_optimizer.update(rho_values, dJdrho, g_vec, dgdrho[0])
            self.problem.rho_field.vector.array = rho_new.copy()
            opt_time = time.perf_counter() - opt_time
            self.optimization_time += opt_time
            
            if self.comm.rank == 0 and self.verbose > 0:
                print(f"Iter: {self.num_iter:3d}, analysis time: {fem_sen_time:.3f} s, beta: {beta:2d}, C: {C_value:8.3f}, V: {V_value:.3f}", flush=True)
            
        self.running_time = time.perf_counter() - running_timer
        self.problem.summary()
        super().summary()
        
        self.problem.save_results()