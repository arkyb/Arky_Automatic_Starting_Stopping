from numpy import *
from matplotlib.pyplot import *
from scipy import signal
import control

#set_printoptions(linewidth=200)

### PHYSICAL PARAMETERS
m_rider = 0.0
h_rider = 0.5
#overall
b = 0.767 #wheelbase
c = 0.02285999 #trail
alph = 1.16 # 1.16 is from Fusion 2/14/23 #1.186 #90-rake: measured from forward x
g = 9.81 #gravity
v = 3 # forward speed

#rear wheel
Rrw = 0.15875 # radius of real wheel
mrw = 2.462 #mass of rear wheel
Axx = 0.027708579 #moment of inertia of rear wheel about x
Ayy = 0.033968214 #moment of inertia of rear wheel about y
Azz = 0.027708579 #moment of inertia of rear wheel about z

#rear frame
mrf = 11.065 #kg, rear frame mass
xrf = .3386 #position of rear frame CG
yrf = 0
hrf = .25606 + (h_rider*m_rider + .25606*mrf)/(mrf+m_rider)
mrf = mrf+m_rider

#mass moments of rear frame
Bxx = 0.1+m_rider*h_rider**2
Bxz = -0.017
Byy = 0.31
Bzz = 0.249

#front frame
xff = .62218#position of front frame CG
yff = 0
hff = .46531
mff = 2.2047 #mass of front frame
Cxx = 0.0659
Cxz = 0.02565
Cyy = 0.06293
Czz = 0.03182

#front wheel
Rfw = 0.15875
mfw = 1.486
Dxx,Dyy,Dzz = 0.016724187,0.020502342,0.016724187


#intermediate terms

#total
mt = mrw+mrf+mff+mfw #total mass
xt = (xrf*mrf+xff*mff+b*mfw)/mt #CG location x
ht = (Rrw*mrw+hrf*mrf+hff*mff+Rfw*mfw)/mt #CG location z
#print("CoM: x: "+str(xt)+", z: "+str(ht))
print("Total mass: "+str(mt)+", Z CoM: "+str(ht))
Txx = Axx+Bxx+Cxx+Dxx+mrw*Rrw**2 +mrf*hrf**2 + mff*hff**2 +mfw*Rfw**2 #global inertia XX
Txz = -(Bxz+Cxz-mrf*xrf*hrf-mff*xff*hff-mfw*b*Rfw)
Tzz = Azz+Bzz+Czz+Dzz+mrf*xrf**2+mfw*b**2
#print("Txx: "+str(Txx))
#print("Txz: "+str(Txz))
#print("Tzz: "+str(Tzz))

#front frame
mf = mff+mfw
xf = (mff*xff+b*mfw)/mf
hf = (hff*mff+Rfw*mfw)/mf
Fxx = Cxx+Dxx+mff*(hff-hf)**2+mfw*(Rfw-hf)**2
Fxz =  Cxz - mff*(xff-xf)*(hff-hf)+mfw*(b-xf)*(-Rfw+hf)
Fzz = Czz+Dzz+mff*(xff-xf)**2 + mfw*(b-xf)**2
#print("Fxx: "+str(Fxx))
#print("Fxz: "+str(Fxz))
#print("Fzz: "+str(Fzz))

#steering frame
lam = pi/2 - alph #angle of steering axis with global z in the vertical plane
u = (xf-b-c)*cos(lam)+hf*sin(lam) #perp distance that CM of front is ahead of steering axis
Fll = mf*u**2+Fxx*sin(lam)**2-2*Fxz*sin(lam)*cos(lam)+Fzz*cos(lam)**2
#print("ff lambda inertia: "+str(Fll))
#print("u: "+str(u))

Flx = -(-mf*u*hf - Fxx*sin(lam) + Fxz*cos(lam))
Flz = -(mf*u*xf - Fxz*sin(lam)+Fzz*cos(lam))
#print("Flx: "+str(Flx))
#print("Flz: "+str(Flz))

f = c*cos(lam)/b #ratio of trail to wheelbase

#angular momentum divided by speed:
Sr = Ayy/Rrw
Sf = Dyy/Rfw
St = Sr+Sf

#frequently appearing static moment term:
Su = mf*u+f*mt*xt

#MDK form:
M = array([[Txx,-(Flx+f*Txz)],[-(Flx+f*Txz), Fll-2*f*Flz+f**2*Tzz]])
K0 = array([[-g*mt*ht,g*Su],[g*Su,-g*Su*sin(lam)]])
K1 = array([[0,-(St+mt*ht)*cos(lam)/b],[0,(Su+Sf*sin(lam))*cos(lam)/b]])
K = K0+K1*v**2
D1 = array([[0,-(f*St+Sf*cos(lam)+Txz*cos(lam)/b+f*mt*ht)],[f*St+Sf*cos(lam),-Flz*cos(lam)/b+f*(Su+Tzz*cos(lam)/b)]])
D = v*D1
#print("############ MDK FORM #########")
print("M: "+str(M))
print("K: "+str(K))
print("D1: "+str(D1))
# print(M)
# print(K)
# print(D1)


######### CONTROLLER DESIGN ########

def getStateSpaceForControlNoIntegrator(v):
    """ sys, A,B,C,D = getStateSpaceForControlWithIntegrator(velocity)
        Creates open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate]'
    """
    #MDK form:
    #MDK form:
    M = array([[Txx,(Flx+f*Txz)],[(Flx+f*Txz), Fll+2*f*Flz+f**2*Tzz]])
    K0 = array([[-g*mt*ht,g*Su],[g*Su,-g*Su*sin(lam)]])
    K1 = array([[0,-(St+mt*ht)*cos(lam)/b],[0,(Su+Sf*sin(lam))*cos(lam)/b]])
    K = K0+K1*v**2
    D1 = array([[0,-(f*St+Sf*cos(lam)-Txz*cos(lam)/b+f*mt*ht)],[f*St+Sf*cos(lam),Flz*cos(lam)/b+f*(Su+Tzz*cos(lam)/b)]])
    D = v*D1


    #### States are [roll steer rollrate steerrate]' by default
    ### augment the state vector to also include the integral of error between a new input goalRoll and roll
    Ass = hstack((zeros((2,2)),eye(2,2)))
    Ass = vstack((Ass,hstack((dot(-linalg.inv(M),K),dot(-linalg.inv(M),D)))))
    Bss = vstack((zeros((2,2)),dot(linalg.inv(M),eye(2))))
    #only take the second (steer) column in the B matrix
    Bss1 = vstack(Bss[:,1])

    # Css = array([[1,0,0,0,0],[0,1,0,0,0]]) #outputs of roll and steer angle.
    Css = eye(4)
    # Dss = array([[0,0],[0,0]])
    Dss = zeros((4,1))
    Dss_control = Dss
    sys = control.StateSpace(Ass,Bss1,Css,Dss_control)
    #print("####### STATE SPACE For control w/integrator #########")
    #print("A: "+str(Ass))
    #print("B: "+str(Bss1))
    #print("C: "+str(Css))
    #print("D: "+str(Dss_control))
    return sys,Ass,Bss1,Css,Dss_control



def getStateSpaceForYawControlnoactuator(velocity):
    """ sys, A,B,C,D = getStateSpaceForYawControlnoactuator(velocity)
            Creates open-loop state space model for motorcycle.
            States are x = [roll steer roll_rate steer_rate]'
        """
    # first get the baseline model
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlNoIntegrator(velocity)
    # for new A matrix, pad with zeros on right
    Aaug = hstack((sys.A, zeros((4, 1))))
    # and add the terms for yaw rate on the bottom
    Aaug = vstack((Aaug, zeros((1, 5))))
    # fill in the two entries
    Aaug[4,1] = v*sin(alph)/b
    Aaug[4,3] = c*sin(alph)/b
    # now augment B matrix with an empty row
    Baug = vstack((sys.B,0))
    # output everything, I assume?
    Caug = eye(5)
    # for D matrix, use 0
    Daug = 0
    sys = control.StateSpace(Aaug, Baug, Caug, Daug)
    return sys, Aaug, Baug, Caug, Daug

def getLQRRazor (velocity,R=0.01,dT=0.005):
    """
    sys,Klqr = getYawLQRRazor(velocity)
    Using hard-coded weights, design a linear quadratic regulator for motorcycle.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired yaw angle and roll angle as the inputs.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForYawControlnoactuator(velocity)
    import control.matlab as cnt
    Q = eye(5)/100.0
    Q[0,0] = 1.0
    Q[1,1] = 0.0
    #give yaw some beans but maybe turn back down a bit if things get too torquey.
    Q[4,4] = 1.0

    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.lqr(sys, Q, R)
    #print("######### LQR GAINS for YAW control #############")
    #print(Klqr)
    return sys,Klqr



if __name__ == '__main__':
    pass
    #main()
