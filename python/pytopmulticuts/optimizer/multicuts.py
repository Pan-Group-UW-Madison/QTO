from .optimizer import Optimizer, DensityFilter, Heaviside
from .sensitivity import Sensitivity

from .milp import MilpOptimizer
from .dw import DWOptimizer

import time
import numpy as np

class Cuts:
    def __init__(self):
        self.rho = None
        self.weight = None
        self.obj = None
        self.d = None

class MulticutsOptimizer(Optimizer):
    def __init__(self, descriptor, problem):
        super().__init__(descriptor, problem)
        
        self.optimizer_name = "multicuts"
        
        self.max_iter = descriptor["max_iter"]
        self.opt_tol = descriptor["opt_tol"]
        self.vol_frac = descriptor["vol_frac"]
        self.radius = descriptor["filter_radius"]
        
        self.verbose = 1
        
        num_consts = 1 if self.problem.objective == "compliance" else 2
        
        if descriptor["subproblem_solver"] == "milp":
            self.sub_optimizer = MilpOptimizer(problem)
        elif descriptor["subproblem_solver"] == "dw":
            self.sub_optimizer = DWOptimizer(problem)
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
        
    def solve_prime(self):
        self.density_filter.forward()
        
        self.problem.solve_prime()
        
        [C_value, V_value, U_value], sensitivities = self.sens_problem.evaluate()
        [dCdrho, dVdrho, dUdrho] = self.density_filter.backward(sensitivities)
        if self.problem.objective == "compliance":
            g_vec = np.array([V_value-self.vol_frac])
            dJdrho, dgdrho = dCdrho, np.vstack([dVdrho])
        else:
            g_vec = np.array([V_value-self.vol_frac, C_value-opt["compliance_bound"]])
            dJdrho, dgdrho = dUdrho, np.vstack([dVdrho, dCdrho])
            
        dJdrho = dJdrho * (self.problem.rho_field.vector.array + self.problem.eps)
            
        return C_value, V_value, dJdrho, g_vec, dgdrho[0]
            
    def solve(self):
        self.density_filter = DensityFilter(self.problem, self.radius)
        self.sens_problem = Sensitivity(self.problem)
        
        running_timer = time.perf_counter()
        
        self.analysis_time, self.optimization_time = 0, 0
        
        self.num_iter = 0
        
        stage = 1
        
        while stage < 3:
            # jump start
            fem_sen_time = time.perf_counter()
            C_value, V_value, dJdrho, g_vec, dgdrho = self.solve_prime()
            fem_sen_time = time.perf_counter() - fem_sen_time
            self.analysis_time += fem_sen_time
            
            opt_time = time.perf_counter()
            rho_values = self.problem.rho_field.vector.array.copy()
            rho_new, cost = self.sub_optimizer.update(rho_values, C_value, dJdrho, self.vol_frac, 1.0)
            self.problem.rho_field.vector.array = rho_new.copy()
            opt_time = time.perf_counter() - opt_time
            self.optimization_time += opt_time
            
            num_inner_iter = 0
            while num_inner_iter < 10:
                num_inner_iter += 1
                
                fem_sen_time = time.perf_counter()
                C_value, V_value, dJdrho, g_vec, dgdrho = self.solve_prime()
                fem_sen_time = time.perf_counter() - fem_sen_time
                self.analysis_time += fem_sen_time
                
                opt_time = time.perf_counter()
                rho_values = self.problem.rho_field.vector.array.copy()
                rho_new, cost = self.sub_optimizer.update(rho_values, C_value, dJdrho, self.vol_frac, 1.0)
                self.problem.rho_field.vector.array = rho_new.copy()
                opt_time = time.perf_counter() - opt_time
                self.optimization_time += opt_time
                
                if self.comm.rank == 0 and self.verbose > 0:
                    print(f"Iter: {self.num_iter+num_inner_iter:3d}, analysis time: {fem_sen_time:.3f} s, optimization time: {opt_time:.3f} s, C: {C_value:8.3f}, Cost: {cost:8.3f}, V: {V_value:.3f}", flush=True)
            
            stage += 1
            
            self.num_iter += num_inner_iter
            
            if stage == 2:
                self.problem.eps = 1e-9
    
        self.running_time = time.perf_counter() - running_timer
        self.problem.summary()
        super().summary()
        
        self.problem.save_results()