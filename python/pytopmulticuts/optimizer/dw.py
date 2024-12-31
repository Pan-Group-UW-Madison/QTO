from mpi4py import MPI
import numpy as np

from .optimizer import SubOptimizer

class DWOptimizer(SubOptimizer):
    def __init__(self, problem):
        super().__init__(problem)
        
    def update(self):
        pass