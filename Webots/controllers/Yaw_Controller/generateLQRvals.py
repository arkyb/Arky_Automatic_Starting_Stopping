from RazorModel import getLQRRazor
from numpy import *

class generateLQRvals:
    def __init__(self,dT=0.005):
        self.dT = dT
    def update(self):
        vels = arange(1,7,0.5)
        for k in range(0,len(vels)):
            sys,Klqr = getLQRRazor(vels[k],0.01,self.dT)
            allGains = ravel(Klqr)
            if(k==0):
                gainmtx = hstack((vels[k],allGains))
            else:
                gainmtx = vstack((gainmtx,hstack((vels[k],allGains))))
        savetxt('razorGains_lqr.txt',gainmtx,comments='# generated by generateLQRvals.py referencing RazorModel.py')

        return gainmtx
