from numpy import *

class DelaySim:
    def __init__(self,N):
        self.vec = zeros(N)
    def update(self,new):
        self.vec[0:-1] = self.vec[1:]
        self.vec[-1] = new
        return float(self.vec[0])
