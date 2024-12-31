from mpi4py import MPI
import numpy as np

from .optimizer import SubOptimizer

class OCOptimizer(SubOptimizer):
    def __init__(self, problem, move):
        super().__init__(problem)
        
        self.move = move
        num_elems = self.problem.rho_field.vector.size
        self.rho_min = np.zeros(num_elems)
        self.rho_max = np.ones(num_elems)
        
    def update(self, rho, dJdrho, g, dgdrho):
        lb, ub = 0.0, 1e6
        comm = MPI.COMM_WORLD
        while ub-lb > 1e-4:
            mid = (lb+ub) / 2.0
            rho_new = np.maximum.reduce([np.minimum.reduce(
                [rho*(-dJdrho/(dgdrho+1e-12)/mid)**0.5, rho+self.move, self.rho_max]), rho-self.move, self.rho_min])
            dg = comm.allreduce(dgdrho@(rho_new-rho), op=MPI.SUM)
            if g + dg > 0:
                lb = mid
            else:
                ub = mid
        return rho_new