from mpi4py import MPI
import numpy as np

import gurobipy as gp

from .optimizer import SubOptimizer

class MilpOptimizer(SubOptimizer):
    def __init__(self, problem):
        super().__init__(problem)
        
        rho_local_size = problem.rho_field.vector.array.size
        rho_size_by_rank = self.comm.allgather(rho_local_size)
        
        self.rho_offset = np.cumsum([0] + rho_size_by_rank)
        
        self.options = {
            "outputflag": 0,
            "MIPGap": 1e-9,
            "FeasibilityTol": 1e-9,
            "OptimalityTol": 1e-9,
        }
        
    def update(self, rho, obj, weight, vol_frac, d):
        rho_global = self.comm.gather(rho, root=0)
        weight_global = self.comm.gather(weight, root=0)
        
        if self.comm.rank == 0:
            rho_global = np.concatenate(rho_global)
            weight_global = np.concatenate(weight_global)
            
            with gp.Env(params=self.options) as env, gp.Model(env=env) as model:
                # Formulate problem
                x = model.addVars(rho_global.size, vtype=gp.GRB.BINARY, name="x")
                
                # cut
                model.setObjective(gp.quicksum(weight_global[i]*x[i] for i in range(rho_global.size)), gp.GRB.MINIMIZE)
                # mass/volume constraint
                model.addConstr(gp.quicksum(x[i] for i in range(rho_global.size)) == vol_frac*rho_global.size)
                # trust region constraint
                model.addConstr(gp.quicksum((1 - 2*rho_global[i])*x[i]+rho_global[i]**2 for i in range(rho_global.size)) <= d*rho_global.size)
                
                model.optimize()
            
                rho_global_new = np.array([x[i].X for i in range(rho_global.size)])
            
        rho_global_new = self.comm.bcast(rho_global_new, root=0)
        
        rho_new = rho_global_new[self.rho_offset[self.comm.rank]:self.rho_offset[self.comm.rank+1]]
        
        cost = np.array([np.dot(weight, rho_new-rho)], dtype='d')
        self.comm.Allreduce(MPI.IN_PLACE, cost, op=MPI.SUM)
        cost = obj + cost[0]
        
        return rho_new, cost