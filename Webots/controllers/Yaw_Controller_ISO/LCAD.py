from numpy import *

class LCAD:
    def __init__(self):
        self.maxRoll = zeros(3)
        self.maxRollRate = zeros(3)
        self.oldRoll = 0
        self.cycle = 0
    def update(self, roll, rollRate):
        if roll > self.maxRoll[-1]:
            self.maxRoll[-1] = roll
        if rollRate < self.maxRollRate[-1]:
            self.maxRollRate[-1] = rollRate
        if self.oldRoll < 0 < roll:
            self.maxRoll[0:-1] = self.maxRoll[1:]
            self.maxRoll[-1] = 0
            self.maxRollRate[0:-1] = self.maxRollRate[1:]
            self.maxRollRate[-1] = 0
        self.oldRoll = roll
        return float(median(self.maxRoll[:])), float(median(self.maxRollRate[:]))
