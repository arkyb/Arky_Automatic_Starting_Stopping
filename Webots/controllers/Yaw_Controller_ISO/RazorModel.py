from numpy import *
from matplotlib.pyplot import *
from scipy import signal
import control

#set_printoptions(linewidth=200)

### PHYSICAL PARAMETERS
m_rider = 0.
h_rider = 0.5
#overall
w = 0.767 #wheelbase
t = 0.02285999 #trail
alph = 1.16 # 1.16 is from Fusion 2/14/23 #1.186 #90-rake: measured from forward x
g = 9.81 #gravity
v = 6.25 # forward speed

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
zrf = .25606 + (h_rider*m_rider + .25606*mrf)/(mrf+m_rider)
mrf = mrf+m_rider

#mass moments of rear frame
Bxx = 0.1+m_rider*h_rider**2
Bxz = 0.017
Byy = 0.31
Bzz = 0.249

#front frame
xff = .62218#position of front frame CG
yff = 0
zff = .46531
mff = 2.2047 #mass of front frame
Cxx = 0.0659
Cxz = -0.02565
Cyy = 0.06293
Czz = 0.03182

#front wheel
Rfw = 0.15875
mfw = 1.486
Dxx,Dyy,Dzz = 0.016724187,0.020502342,0.016724187


#intermediate terms

#total
mt = mrw+mrf+mff+mfw #total mass
xt = (xrf*mrf+xff*mff+w*mfw)/mt #CG location x
zt = (Rrw*mrw+zrf*mrf+zff*mff+Rfw*mfw)/mt #CG location z
print("CoM: x: "+str(xt)+", z: "+str(zt))
Txx = Axx+Bxx+Cxx+Dxx+mrw*Rrw**2 +mrf*zrf**2 + mff*zff**2 +mfw*Rfw**2 #global inertia XX
Txz = Bxz+Cxz+mrf*xrf*zrf+mff*xff*zff+mfw*w*Rfw
Tzz = Azz+Bzz+Czz+Dzz+mrf*xrf**2+mfw*w**2
print("Txx: "+str(Txx))
print("Txz: "+str(Txz))
print("Tzz: "+str(Tzz))

#front frame
mf = mff+mfw
xf = (mff*xff+w*mfw)/mf
zf = (zff*mff+Rfw*mfw)/mf
Fxx = Cxx+Dxx+mff*(zff-zf)**2+mfw*(Rfw-zf)**2
Fxz = Cxz - mff*(xff-xf)*(-zff+zf)+mfw*(w-xf)*(Rfw-zf)
Fzz = Czz+Dzz+mff*(xff-xf)**2 + mfw*(w-xf)**2
print("Fxx: "+str(Fxx))
print("Fxz: "+str(Fxz))
print("Fzz: "+str(Fzz))

#steering frame
lam = pi/2 - alph #angle of steering axis with global z in the vertical plane
u = (xf-w-t)*cos(lam)+zf*sin(lam) #perp distance that CM of front is ahead of steering axis
Fll = mf*u**2+Fxx*sin(lam)**2-2*Fxz*sin(lam)*cos(lam)+Fzz*cos(lam)**2
print("ff lambda inertia: "+str(Fll))
print("u: "+str(u))

Flx = (mf*u*zf + Fxx*sin(lam) + Fxz*cos(lam))
Flz = mf*u*xf + Fxz*sin(lam)+Fzz*cos(lam)
print("Flx: "+str(Flx))
print("Flz: "+str(Flz))

f = t*cos(lam)/w #ratio of trail to wheelbase

#angular momentum divided by speed:
Sr = Ayy/Rrw
Sf = Dyy/Rfw
St = Sr+Sf

#frequently appearing static moment term:
Su = mf*u+f*mt*xt

#MDK form:
M = array([[Txx,-(Flx+f*Txz)],[-(Flx+f*Txz), Fll+2*f*Flz+f**2*Tzz]])
K0 = array([[-g*mt*zt,g*Su],[g*Su,-g*Su*sin(lam)]])
K1 = array([[0,-(St+mt*zt)*cos(lam)/w],[0,(Su+Sf*sin(lam))*cos(lam)/w]])
K = K0+K1*v**2
D1 = array([[0,-(f*St+Sf*cos(lam)+Txz*cos(lam)/w+f*mt*zt)],[f*St+Sf*cos(lam),-Flz*cos(lam)/w+f*(Su+Tzz*cos(lam)/w)]])
D = v*D1
print("############ MDK FORM #########")
print("M: "+str(M))
print("K: "+str(K))
print("D1: "+str(D1))
# print(M)
# print(K)
# print(D1)

#now translate the model into state space, given this forward velocity
def getStateSpace(v):
    """ sys, A,B,C,D = getStateSpace(velocity)
        Creates open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate]'
    """
    #MDK form:
    M = array([[Txx,Flx+f*Txz],[Flx+f*Txz, Fll+2*f*Flz+f**2*Tzz]])
    K0 = array([[g*mt*zt,-g*Su],[-g*Su,-g*Su*sin(lam)]])
    K1 = array([[0,(St-mt*zt)*cos(lam)/w],[0,(Su+Sf*sin(lam))*cos(lam)/w]])
    K = K0+K1*v**2
    D1 = array([[0,f*St+Sf*cos(lam)+Txz*cos(lam)/w-f*mt*zt],[-(f*St+Sf*cos(lam)),Flz*cos(lam)/w+f*(Su+Tzz*cos(lam)/w)]])
    D = v*D1
    Ass = hstack((zeros((2,2)),eye(2,2)))
    Ass = vstack((Ass,hstack((dot(-linalg.inv(M),K),dot(-linalg.inv(M),D)))))
    Bss = vstack((zeros((2,2)),dot(linalg.inv(M),eye(2))))
    Css = array([[1,0,0,0],[0,1,0,0]]) #output is roll.
    Dss = array([[0,0],[0,0]])
    sys = control.StateSpace(Ass,Bss,Css,Dss)
    print("####### STATE SPACE #########")
    print("A: "+str(Ass))
    print("B: "+str(Bss))
    print("C: "+str(Css))
    print("D: "+str(Dss))
    return sys,Ass,Bss,Css,Dss



######### CONTROLLER DESIGN ########

#now translate the model into state space, given this forward velocity
def getStateSpaceForControl(v):
    """ sys, A,B,C,D = getStateSpaceForControl(velocity)
        Creates open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate]'
    """
    #MDK form:
    Mmdk = array([[Txx,Flx+f*Txz],[Flx+f*Txz, Fll+2*f*Flz+f**2*Tzz]])
    K0mdk = array([[g*mt*zt,-g*Su],[-g*Su,-g*Su*sin(lam)]])
    K1mdk = array([[0,(St-mt*zt)*cos(lam)/w],[0,(Su+Sf*sin(lam))*cos(lam)/w]])
    Kmdk = K0mdk+K1mdk*v**2
    D1mdk = array([[0,f*St+Sf*cos(lam)+Txz*cos(lam)/w-f*mt*zt],[-(f*St+Sf*cos(lam)),Flz*cos(lam)/w+f*(Su+Tzz*cos(lam)/w)]])
    Dmdk = v*D1mdk


    #### States are [roll steer rollrate steerrate]'
    Ass = hstack((zeros((2,2)),eye(2,2)))
    Ass = vstack((Ass,hstack((dot(-linalg.inv(Mmdk),Kmdk),dot(-linalg.inv(Mmdk),Dmdk)))))
    Bss = vstack((zeros((2,2)),dot(linalg.inv(Mmdk),eye(2))))
    #get rid of the first column, which is roll torque (don't have without rider)
    Bss_control = vstack(Bss[:,1])
    Css = array([[1,0,0,0],[0,1,0,0]]) #outputs of roll and steer angle.
    Dss = array([[0,0],[0,0]])
    Dss_control = 0
    print("####### STATE SPACE for control #########")
    print("A: "+str(Ass))
    print("B: "+str(Bss_control))
    print("C: "+str(Css))
    print("D: "+str(Dss_control))
    sys = control.StateSpace(Ass,Bss_control,Css,Dss_control)

    return sys,Ass,Bss_control,Css,Dss_control

def getStateSpaceForControlWithIntegrator(v):
    """ sys, A,B,C,D = getStateSpaceForControlWithIntegrator(velocity)
        Creates open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate integral_roll_error]'
    """
    #MDK form:
    M = array([[Txx,Flx+f*Txz],[Flx+f*Txz, Fll+2*f*Flz+f**2*Tzz]])
    K0 = array([[g*mt*zt,-g*Su],[-g*Su,-g*Su*sin(lam)]])
    K1 = array([[0,(St-mt*zt)*cos(lam)/w],[0,(Su+Sf*sin(lam))*cos(lam)/w]])
    K = K0+K1*v**2
    D1 = array([[0,f*St+Sf*cos(lam)+Txz*cos(lam)/w-f*mt*zt],[-(f*St+Sf*cos(lam)),Flz*cos(lam)/w+f*(Su+Tzz*cos(lam)/w)]])
    D = v*D1


    #### States are [roll steer rollrate steerrate]' by default
    ### augment the state vector to also include the integral of error between a new input goalRoll and roll
    Ass = hstack((zeros((2,2)),eye(2,2)))
    Ass = vstack((Ass,hstack((dot(-linalg.inv(M),K),dot(-linalg.inv(M),D)))))
    #now add the integrator on the right end:
    Ass = hstack((Ass,zeros((4,1))))
    #now add the last row: the derivative of this integrator state is the goalRoll - roll
    Ass = vstack((Ass,array([-1, 0, 0, 0, 0])))

    Bss = vstack((zeros((2,2)),dot(linalg.inv(M),eye(2))))

    #get rid of the first column, which is roll torque (don't have without rider)
    # print(Bss[:,1])
    Bss_control = vstack((vstack(Bss[:,1]),array([0])))
    #now add back on an input to be the goal roll angle:
    #Bss_control = hstack((Bss_control,array([[0],[0],[0],[0],[1]])))

    # Css = array([[1,0,0,0,0],[0,1,0,0,0]]) #outputs of roll and steer angle.
    Css = eye(5)
    # Dss = array([[0,0],[0,0]])
    Dss = zeros((5,1))
    Dss_control = Dss
    sys = control.StateSpace(Ass,Bss_control,Css,Dss_control)
    print("####### STATE SPACE For control w/integrator #########")
    print("A: "+str(Ass))
    print("B: "+str(Bss_control))
    print("C: "+str(Css))
    print("D: "+str(Dss_control))
    return sys,Ass,Bss_control,Css,Dss_control

def getStateSpaceForControlWithActuator(v,tau_motor=0.025):
    """ sys, A,B,C,D = getStateSpaceForControlWithIntegrator(velocity)
        Creates open-loop state space model for motorcycle.
        simulates motor as a first order filter between commanded and actual torque.
        States are x = [roll steer roll_rate steer_rate command_torque]'
    """
    #MDK form:
    M = array([[Txx,Flx+f*Txz],[Flx+f*Txz, Fll+2*f*Flz+f**2*Tzz]])
    K0 = array([[g*mt*zt,-g*Su],[-g*Su,-g*Su*sin(lam)]])
    K1 = array([[0,(St-mt*zt)*cos(lam)/w],[0,(Su+Sf*sin(lam))*cos(lam)/w]])
    K = K0+K1*v**2
    D1 = array([[0,f*St+Sf*cos(lam)+Txz*cos(lam)/w-f*mt*zt],[-(f*St+Sf*cos(lam)),Flz*cos(lam)/w+f*(Su+Tzz*cos(lam)/w)]])
    D = v*D1


    #### States are [roll steer rollrate steerrate]' by default
    ### augment the state vector to also include the integral of error between a new input goalRoll and roll
    Ass = hstack((zeros((2,2)),eye(2,2)))
    Ass = vstack((Ass,hstack((dot(-linalg.inv(M),K),dot(-linalg.inv(M),D)))))
    Bss = vstack((zeros((2,2)),dot(linalg.inv(M),eye(2))))
    #only take the second (steer) column in the B matrix
    Bss1 = vstack(Bss[:,1])
    #now stack this B matrix alongside A to represent adding true torque as a state
    Ass = hstack((Ass,Bss1))
    #now add the last row: torquedot = 1/tau*(command - torque)
    Ass = vstack((Ass,array([0, 0, 0, 0, -1/tau_motor])))
    # print(Bss[:,1])
    Bss_control = vstack(array([0,0,0,0,1/tau_motor]))
    #now add back on an input to be the goal roll angle:
    #Bss_control = hstack((Bss_control,array([[0],[0],[0],[0],[1]])))

    # Css = array([[1,0,0,0,0],[0,1,0,0,0]]) #outputs of roll and steer angle.
    Css = eye(5)
    # Dss = array([[0,0],[0,0]])
    Dss = zeros((5,1))
    Dss_control = Dss
    sys = control.StateSpace(Ass,Bss_control,Css,Dss_control)
    #print("####### STATE SPACE For control w/integrator #########")
    #print("A: "+str(Ass))
    #print("B: "+str(Bss_control))
    #print("C: "+str(Css))
    #print("D: "+str(Dss_control))
    return sys,Ass,Bss_control,Css,Dss_control

def getStateSpaceForControlWithActuatorMotorInertia(v,tau_motor=0.025,motor_inertia = 2e-4, motor_damping = 0.0017):
    """ sys, A,B,C,D = getStateSpaceForControlWithIntegrator(velocity)
        Creates open-loop state space model for motorcycle.
        simulates motor as a first order filter between commanded and actual torque.
        States are x = [roll steer roll_rate steer_rate command_torque]'
    """

    #MDK form:
    M = array([[Txx,Flx+f*Txz],[Flx+f*Txz, motor_inertia + Fll+2*f*Flz+f**2*Tzz]])
    K0 = array([[g*mt*zt,-g*Su],[-g*Su,-g*Su*sin(lam)]])
    K1 = array([[0,(St-mt*zt)*cos(lam)/w],[0,(Su+Sf*sin(lam))*cos(lam)/w]])
    K = K0+K1*v**2
    D1 = array([[0,f*St+Sf*cos(lam)+Txz*cos(lam)/w-f*mt*zt],[-(f*St+Sf*cos(lam)),Flz*cos(lam)/w+f*(Su+Tzz*cos(lam)/w)]])
    D = v*D1
    D[1,1]+=motor_damping


    #### States are [roll steer rollrate steerrate]' by default
    ### augment the state vector to also include the integral of error between a new input goalRoll and roll
    Ass = hstack((zeros((2,2)),eye(2,2)))
    Ass = vstack((Ass,hstack((dot(-linalg.inv(M),K),dot(-linalg.inv(M),D)))))
    Bss = vstack((zeros((2,2)),dot(linalg.inv(M),eye(2))))
    #only take the second (steer) column in the B matrix
    Bss1 = vstack(Bss[:,1])
    #now stack this B matrix alongside A to represent adding true torque as a state
    Ass = hstack((Ass,Bss1))
    #now add the last row: torquedot = 1/tau*(command - torque)
    Ass = vstack((Ass,array([0, 0, 0, 0, -1/tau_motor])))
    # print(Bss[:,1])
    Bss_control = vstack(array([0,0,0,0,1/tau_motor]))
    #now add back on an input to be the goal roll angle:
    #Bss_control = hstack((Bss_control,array([[0],[0],[0],[0],[1]])))

    # Css = array([[1,0,0,0,0],[0,1,0,0,0]]) #outputs of roll and steer angle.
    Css = eye(5)
    # Dss = array([[0,0],[0,0]])
    Dss = zeros((5,1))
    Dss_control = Dss
    sys = control.StateSpace(Ass,Bss_control,Css,Dss_control)
    print("####### STATE SPACE For control w/integrator #########")
    print("A: "+str(Ass))
    print("B: "+str(Bss_control))
    print("C: "+str(Css))
    print("D: "+str(Dss_control))
    return sys,Ass,Bss_control,Css,Dss_control

def getStateSpaceForControlWithActuatorLane(velocity,tau_motor=0.05):
    """ sys, A,B,C,D = getStateSpaceForControlWithIntegratorLane(velocity,tau_motor=0.05)
        Creates open-loop state space model for motorcycle.
        simulates motor as a first order filter between commanded and actual torque.
        States are x = [roll steer roll_rate steer_rate command_torque yaw y]'
    """
    #first get the baseline model
    sys1,A,B,C,D = getStateSpaceForControlWithActuator(velocity,tau_motor)
    #for new A matrix, pad with zeros on right
    Aaug = hstack((sys1.A,zeros((5,2))))
    #pad with zeros on bottom
    Aaug = vstack((Aaug,zeros((2,7))))
    #fill in the two entries for y and yaw.
    Aaug[5,1] = velocity/w*sin(alph) ### AAB ADDED SIN(alph) JULY 2023!!
    Aaug[6,5] = velocity

    #B matrix needs three more zeros on the end
    Baug = vstack((sys1.B,zeros((2,1))))


    #for C matrix, output everything:
    Caug = eye(7)
    #for D matrix, use 0
    Daug = 0

    sys = control.StateSpace(Aaug,Baug,Caug,Daug)

    return sys,Aaug,Baug,Caug,Daug

def getStateSpaceForControlNoIntegrator(v):
    """ sys, A,B,C,D = getStateSpaceForControlWithIntegrator(velocity)
        Creates open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate]'
    """
    #MDK form:
    M = array([[Txx,-(Flx+f*Txz)],[-(Flx+f*Txz), Fll+2*f*Flz+f**2*Tzz]])
    K0 = array([[-g*mt*zt,g*Su],[g*Su,g*Su*cos(lam)]])
    K1 = array([[0,-(St+mt*zt)*sin(lam)/w],[0,(Su+Sf*cos(lam))*sin(lam)/w]])
    K = K0+K1*v**2
    D1 = array([[0,-(f*St+Sf*sin(lam)+Txz*sin(lam)/w+f*mt*zt)],[f*St+Sf*sin(lam),Flz*sin(lam)/w+f*(Su+Tzz*sin(lam)/w)]])
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
    print("####### STATE SPACE For control w/integrator #########")
    print("A: "+str(Ass))
    print("B: "+str(Bss1))
    print("C: "+str(Css))
    print("D: "+str(Dss_control))
    return sys,Ass,Bss1,Css,Dss_control

def getStateSpaceForControlWithRoadCoords(velocity):
    """ sys, A,B,C,D = getStateSpaceForControlWithRoadError(velocity)
        Creates open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate integral_roll_error yaw y integral_y_error]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
    """
    #first get the baseline model
    sys1,A,B,C,D = getStateSpaceForControlWithIntegrator(velocity)
    #now augment that model to include states for yaw and y.
    #yaw rate is VERY APPROXIMATELY, at steady state, equal to g/U*roll
    #this is because at SS, moment balance gives m*g*sin(roll)*h = m*cos(roll)*yawrate*U*h
    #that leads to tan(roll) = (U*yawrate)/g, and under small roll, the above falls out.

    #HOWEVER it's also true that even NOT at SS, yaw rate is more directly given by
    # delta = wheelbase/R, and U/R = yawrate so R = U/yawrate, and
    #therefore delta = wheelbase/U*yawrate so yawrate = delta*U/wheelbase

    #then, dy/dt = U*sin(yaw) and again small angle gives us U*yaw = dy/dt.

    #for new A matrix, pad with zeros on right
    Aaug = hstack((sys1.A,zeros((5,3))))
    #pad with zeros on bottom
    Aaug = vstack((Aaug,zeros((3,8))))
    #fill in the two entries for y and yaw.
    Aaug[5,1] = velocity/w*sin(alph) ### AAB ADDED SIN(alph) JULY 2023!!
    Aaug[6,5] = velocity
    #this last one we'll use as integral of lane error.
    Aaug[7,6] = -1

    #B matrix needs three more zeros on the end
    Baug = vstack((sys1.B,zeros((3,1))))


    #for C matrix, output everything:
    Caug = eye(8)
    #for D matrix, use 0
    Daug = 0

    sys = control.StateSpace(Aaug,Baug,Caug,Daug)

    return sys,Aaug,Baug,Caug,Daug

def getStateSpaceForPreviewControl(velocity,dT,Np):
    """ sys, A,B,C,D = getStateSpaceForPreviewControl(velocity)
        Creates a DISCRETE open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate yaw y yr0 yr1 ... yrnp-1]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
    """
    #get nominal lane error system
    sysc,Ac,Bc,Cc,Dc = getStateSpaceForControl(velocity)
    #we want to add yaw rate and lateral position as states:
    #for new A matrix, pad with zeros on right
    Aaug = hstack((sysc.A,zeros((4,2))))
    #pad with zeros on bottom
    Aaug = vstack((Aaug,zeros((2,6))))
    #fill in the two entries for y and yaw.
    Aaug[4,1] = velocity/w
    Aaug[5,4] = velocity
    #B matrix needs tWO more zeros on the end
    Baug = vstack((sysc.B,zeros((2,1))))
    #for C matrix, output everything:
    Caug = eye(6)
    #for D matrix, use 0
    Daug = 0
    sysa = control.StateSpace(Aaug,Baug,Caug,Daug)
    #now discretize the system
    import control.matlab as cnt
    sysd = cnt.c2d(sysa,dT)
    #now augment the system with the road states.
    Aprev = hstack((sysd.A,zeros((6,Np))))
    Aprev = vstack((Aprev,zeros((Np,Np+6))))
    #fill in the identity portion of the preview controller:
    #identity matrix should be of dimention Np-1, and begin at 8th column, 7th row.
    Aprev[6:-1,7:] = eye(Np-1)

    B1prev = vstack((sysd.B,zeros((Np,1))))
    B2prev = zeros((Np+6,1))
    B2prev[-1,0] = 1 #this is the lone 1 for road input.
    Bprev = hstack((B1prev,B2prev)) #captures both road input and torque input
    Cprev = zeros((6,Np+6))
    Cprev[0:6,0:6] = eye(6)
    Dprev = 0
    sysprev = control.StateSpace(Aprev,Bprev,Cprev,Dprev,dT)

    return sysprev,Aprev,Bprev,Cprev,Dprev

def getPreviewControlLQRRazor(velocity, dT, Np):
    #start with the open-loop system
    sys,A0,B0,C0,D0 = getStateSpaceForPreviewControl(velocity,dT,Np)
    #set up the cost function
    Ccost = zeros((2,Np+6))
    #cost should weight lane error NOW, or y-yr0
    #cost should also weight yaw error NOW, or y -(yr1-yr0)/(U*dT) since the latter is approx psi of road
    #reminder: state vector is [roll steer droll dsteer psi y yr0 yr1 ... yrnp-1]
    # Ccost[:,0:8] = array([[0,0,0,0,0,0,1,0],[0,0,0,0,1,0,-1/(velocity*dT),1/(velocity*dT)]])
    Ccost[:,0:8] = array([[0,0,0,0,0,-1,1,0],[1,0,0,0,0,0,0,0]])

    weights = array([[1,0],[0,1]])
    Q = np.dot(Ccost.T,weights)
    Q = np.dot(Q,Ccost)
    #cost on steering effort
    R = 100
    #now, find the LQR for non-preview portion of the system.
    #need to pull portion of original B that is 8x1 and deals with steer torque
    Bprev = vstack(sys.B[:,0])
    Klqr,Slqr,Elqr = control.dlqr(sys.A,Bprev,Q,R)
    return sys,Klqr

def getStateSpaceForPreviewControlLocal(velocity,dT,Np):
    """ sys, A,B,C,D = getStateSpaceForPreviewControl(velocity)
        Creates a DISCRETE open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate yaw yr0 yr1 ... yrnp-1]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
    """
    #get nominal lane error system
    sysc,Ac,Bc,Cc,Dc = getStateSpaceForControl(velocity)
    #now discretize the system
    import control.matlab as cnt
    sysd = cnt.c2d(sysc,dT)
    #now augment the system with the road states.
    Aprev = hstack((sysd.A,zeros((4,Np))))
    Aprev = vstack((Aprev,zeros((Np,Np+4))))
    #fill in the identity portion of the preview controller:
    #identity matrix should be of dimention Np-1, and begin at 8th column, 7th row.
    Aprev[4:-1,5:] = eye(Np-1)
    #now fill in the -U^2*dT*psidot portion to account for movement of preview points due to yaw misalign.
    #note that psidot in our model is directly velocity*steer/wheelbase
    Aprev[4:-1,1] = -velocity**2*dT/w*vstack(arange(1,Np)) #i think velo should have - on it?

    B1prev = vstack((sysd.B,zeros((Np,1))))
    B2prev = zeros((Np+4,1))
    B2prev[-1,0] = 1 #this is the lone 1 for road input.
    Bprev = hstack((B1prev,B2prev)) #captures both road input and torque input
    Cprev = zeros((4,Np+4))
    Cprev[0:4,0:4] = eye(4)
    Dprev = 0
    sysprev = control.StateSpace(Aprev,Bprev,Cprev,Dprev,dT)

    return sysprev,Aprev,Bprev,Cprev,Dprev

def getPreviewControlLQRRazorLocal(velocity, dT, Np):
    #start with the open-loop system
    sys,A0,B0,C0,D0 = getStateSpaceForPreviewControlLocal(velocity,dT,Np)
    #set up the cost function
    Ccost = zeros((2,Np+4))
    #cost should weight lane error NOW, or y-yr0
    #cost should also weight yaw error NOW, or y -(yr1-yr0)/(U*dT) since the latter is approx psi of road
    #reminder: state vector is [roll steer droll dsteer psi y yr0 yr1 ... yrnp-1]
    # Ccost[:,0:8] = array([[0,0,0,0,0,0,1,0],[0,0,0,0,1,0,-1/(velocity*dT),1/(velocity*dT)]])
    Ccost[:,0:5] = array([[0,0,0,0,1],[1,0,0,0,0]])

    weights = array([[1,0],[0,1]])
    Q = np.dot(Ccost.T,weights)
    Q = np.dot(Q,Ccost)
    #cost on steering effort
    R = 10
    #now, find the LQR for non-preview portion of the system.
    #need to pull portion of original B that is 8x1 and deals with steer torque
    #matlab: [Kf,Sf,Ef] = dlqr(Aaug,[Baug],C'*Q*C,[1])%steering only
    Bprev = vstack(sys.B[:,0])
    Klqr,Slqr,Elqr = control.dlqr(sys.A,Bprev,Q,R)
    return sys,Klqr

def getStateSpaceForPreviewControlRoll(velocity,dT,Np):
    """ sys, A,B,C,D = getStateSpaceForPreviewControl(velocity)
        Creates a DISCRETE open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate roll0 rollp1 ... rollnp-1]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
        The previewed road feature is REQUIRED ROLL
    """
    #get nominal lane error system
    sysc,Ac,Bc,Cc,Dc = getStateSpaceForControl(velocity)
    import control.matlab as cnt
    sysd = cnt.c2d(sysc,dT)
    nstates = sysd.A.shape[0]
    #now augment the system with the road states.
    Aprev = hstack((sysd.A,zeros((nstates,Np))))
    Aprev = vstack((Aprev,zeros((Np,Np+nstates))))
    #fill in the identity portion of the preview controller:
    #identity matrix should be of dimention Np-1, and begin at 8th column, 7th row.
    Aprev[nstates:-1,(nstates+1):] = eye(Np-1)

    B1prev = vstack((sysd.B,zeros((Np,1))))
    B2prev = zeros((Np+nstates,1))
    B2prev[-1,0] = 1 #this is the lone 1 for road input.
    Bprev = hstack((B1prev,B2prev)) #captures both road input and torque input
    Cprev = zeros((nstates,Np+nstates))
    Cprev[0:nstates,0:nstates] = eye(nstates)
    Dprev = 0
    sysprev = control.StateSpace(Aprev,Bprev,Cprev,Dprev,dT)

    return sysprev,Aprev,Bprev,Cprev,Dprev

def getPreviewControlLQRRazorRoll(velocity, dT, Np):
    #start with the open-loop system
    sys,A0,B0,C0,D0 = getStateSpaceForPreviewControlRoll(velocity,dT,Np)
    #set up the cost function
    Ccost = zeros((2,Np+4))
    #cost should weight roll error NOW, or roll-roll0
    #reminder: state vector is [roll steer droll dsteer psi y roll0 roll1 ... rollnp-1]
    # Ccost[:,0:8] = array([[0,0,0,0,0,0,1,0],[0,0,0,0,1,0,-1/(velocity*dT),1/(velocity*dT)]])
    Ccost[:,0:6] = array([[-1,0,0,0,1,0],[0,velocity**2/(w*g),0,0,1,0]])
    #first weight is on "roll angle of road"
    #second weight is on "steering angle of road"
    weights = array([[1,0],[0,1]])
    Q = np.dot(Ccost.T,weights)
    Q = np.dot(Q,Ccost)
    #cost on steering effort
    R = 0.1
    #now, find the LQR for non-preview portion of the system.
    #need to pull portion of original B that is 8x1 and deals with steer torque
    Bprev = vstack(sys.B[:,0])
    Klqr,Slqr,Elqr = control.dlqr(sys.A,Bprev,Q,R)
    return sys,Klqr

def getStateSpaceForPreviewControlRollLane(velocity,dT,Np):
    """ sys, A,B,C,D = getStateSpaceForPreviewControl(velocity)
        Creates a DISCRETE open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate yaw y roll0 rollp1 ... rollnp-1]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
        The previewed road feature is REQUIRED ROLL
    """
    #get nominal lane error system
    sysc,Ac,Bc,Cc,Dc = getStateSpaceForControl(velocity)
    #we want to add yaw rate and lateral position as states:
    #for new A matrix, pad with zeros on right
    Aaug = hstack((sysc.A,zeros((4,2))))
    #pad with zeros on bottom
    Aaug = vstack((Aaug,zeros((2,6))))
    #fill in the two entries for y and yaw.
    Aaug[4,1] = velocity/w
    Aaug[5,4] = velocity
    #B matrix needs tWO more zeros on the end
    Baug = vstack((sysc.B,zeros((2,1))))
    #for C matrix, output everything:
    Caug = eye(6)
    #for D matrix, use 0
    Daug = 0
    sysa = control.StateSpace(Aaug,Baug,Caug,Daug)
    #now discretize the system
    import control.matlab as cnt
    sysd = cnt.c2d(sysa,dT)
    nstates = sysd.A.shape[0]
    #now augment the system with the road states.
    Aprev = hstack((sysd.A,zeros((nstates,Np))))
    Aprev = vstack((Aprev,zeros((Np,Np+nstates))))
    #fill in the identity portion of the preview controller:
    #identity matrix should be of dimention Np-1, and begin at 8th column, 7th row.
    Aprev[nstates:-1,(nstates+1):] = eye(Np-1)

    B1prev = vstack((sysd.B,zeros((Np,1))))
    B2prev = zeros((Np+nstates,1))
    B2prev[-1,0] = 1 #this is the lone 1 for road input.
    Bprev = hstack((B1prev,B2prev)) #captures both road input and torque input
    Cprev = zeros((nstates,Np+nstates))
    Cprev[0:nstates,0:nstates] = eye(nstates)
    Dprev = 0
    sysprev = control.StateSpace(Aprev,Bprev,Cprev,Dprev,dT)

    return sysprev,Aprev,Bprev,Cprev,Dprev

def getPreviewControlLQRRazorRollLane(velocity, dT, Np):
    #start with the open-loop system
    sys,A0,B0,C0,D0 = getStateSpaceForPreviewControlRollLane(velocity,dT,Np)
    #set up the cost function
    Ccost = zeros((3,Np+6))
    #cost should weight roll error NOW, or roll-roll0
    #reminder: state vector is [roll steer droll dsteer psi y roll0 roll1 ... rollnp-1]
    # Ccost[:,0:8] = array([[0,0,0,0,0,0,1,0],[0,0,0,0,1,0,-1/(velocity*dT),1/(velocity*dT)]])
    # Ccost[:,0:8] = array([[-1,0,0,0,0,0,1,0],[0,velocity**2/(w*g),0,0,0,0,1,0]])
    Ccost[:,0:8] = array([[-1,0,0,0,0,0,1,0],[0,0,0,0,0,1,0,0],[0,0,0,1,0,0,0,0]])

    #best design was 3, 1, .5
    weights = array([[5,0,0],[0,5,0],[0,0,1]])
    Q = np.dot(Ccost.T,weights)
    Q = np.dot(Q,Ccost)
    #cost on steering effort
    R = 1
    #now, find the LQR for non-preview portion of the system.
    #need to pull portion of original B that is 8x1 and deals with steer torque
    Bprev = vstack(sys.B[:,0])
    Klqr,Slqr,Elqr = control.dlqr(sys.A,Bprev,Q,R)
    return sys,Klqr

def getStateSpaceForPreviewControlRollFutureLane(velocity,dT,Np):
    """ sys, A,B,C,D = getStateSpaceForPreviewControl(velocity)
        Creates a DISCRETE open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate yaw y roll0 rollp1 ... rollnp-1]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
        The previewed road feature is REQUIRED ROLL
    """
    #get nominal lane error system
    sysc,Ac,Bc,Cc,Dc = getStateSpaceForControl(velocity)
    #we want to add yaw rate and lateral position as states:
    #for new A matrix, pad with zeros on right
    Aaug = hstack((sysc.A,zeros((4,2))))
    #pad with zeros on bottom
    Aaug = vstack((Aaug,zeros((2,6))))
    #fill in the two entries for y and yaw.
    Aaug[4,1] = velocity/w
    #here we use yaw angle to predict y vel. BUT
    #now, we also say that yaw rate * preview arm contributes!
    typrev = 0.5
    yprev = velocity*typrev
    #so now ydot is the rate of change of a previewed y coord.
    Aaug[5,4] = velocity + yprev*velocity/w
    #B matrix needs tWO more zeros on the end
    Baug = vstack((sysc.B,zeros((2,1))))
    #for C matrix, output everything:
    Caug = eye(6)
    #for D matrix, use 0
    Daug = 0
    sysa = control.StateSpace(Aaug,Baug,Caug,Daug)
    #now discretize the system
    import control.matlab as cnt
    sysd = cnt.c2d(sysa,dT)
    nstates = sysd.A.shape[0]
    #now augment the system with the road states.
    Aprev = hstack((sysd.A,zeros((nstates,Np))))
    Aprev = vstack((Aprev,zeros((Np,Np+nstates))))
    #fill in the identity portion of the preview controller:
    #identity matrix should be of dimention Np-1, and begin at 8th column, 7th row.
    Aprev[nstates:-1,(nstates+1):] = eye(Np-1)

    B1prev = vstack((sysd.B,zeros((Np,1))))
    B2prev = zeros((Np+nstates,1))
    B2prev[-1,0] = 1 #this is the lone 1 for road input.
    Bprev = hstack((B1prev,B2prev)) #captures both road input and torque input
    Cprev = zeros((nstates,Np+nstates))
    Cprev[0:nstates,0:nstates] = eye(nstates)
    Dprev = 0
    sysprev = control.StateSpace(Aprev,Bprev,Cprev,Dprev,dT)

    return sysprev,Aprev,Bprev,Cprev,Dprev

def getPreviewControlLQRRazorRollFutureLane(velocity, dT, Np):
    #start with the open-loop system
    sys,A0,B0,C0,D0 = getStateSpaceForPreviewControlRollLane(velocity,dT,Np)
    #set up the cost function
    Ccost = zeros((4,Np+6))
    #cost should weight roll error NOW, or roll-roll0
    #reminder: state vector is [roll steer droll dsteer psi y roll0 roll1 ... rollnp-1]
    # Ccost[:,0:8] = array([[0,0,0,0,0,0,1,0],[0,0,0,0,1,0,-1/(velocity*dT),1/(velocity*dT)]])
    # Ccost[:,0:8] = array([[-1,0,0,0,0,0,1,0],[0,velocity**2/(w*g),0,0,0,0,1,0]])
    Ccost[:,0:8] = array([[-1,0,0,0,0,0,1,0],[0,0,0,0,0,1,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,0,0,0,0]])

    #best design was 3, 1, .5
    qrollE,qlaneE,qsteerrate,qrollrate = 5,.5,1,1
    weights = array([[qrollE,0,0,0],[0,qlaneE,0,0],[0,0,qsteerrate,0],[0,0,0,qrollrate]])
    Q = np.dot(Ccost.T,weights)
    Q = np.dot(Q,Ccost)
    #cost on steering effort
    R = 1
    #now, find the LQR for non-preview portion of the system.
    #need to pull portion of original B that is 8x1 and deals with steer torque
    Bprev = vstack(sys.B[:,0])
    Klqr,Slqr,Elqr = control.dlqr(sys.A,Bprev,Q,R)
    return sys,Klqr

def getStateSpaceForPreviewControlRollLaneI(velocity,dT,Np):
    """ sys, A,B,C,D = getStateSpaceForPreviewControl(velocity)
        Creates a DISCRETE open-loop state space model for motorcycle.
        States are x = [roll steer roll_rate steer_rate yaw y roll_error_int roll0 rollp1 ... rollnp-1]'
        yaw is the Northing-Easting angle of the vehicle measured CCW from positive East (X)
        y is the linearized northing coordinate of the vehicle.
        The previewed road feature is REQUIRED ROLL
    """
    #get nominal lane error system
    sysc,Ac,Bc,Cc,Dc = getStateSpaceForControl(velocity)
    #we want to add yaw rate and lateral position as states:
    #for new A matrix, pad with zeros on right
    Aaug = hstack((sysc.A,zeros((4,3))))
    #pad with zeros on bottom
    Aaug = vstack((Aaug,zeros((3,7))))
    #fill in the two entries for y and yaw.
    Aaug[4,1] = velocity/w
    Aaug[5,4] = velocity
    #fill in last entry to be integral of roll
    Aaug[6,0] = 1
    #B matrix needs tWO more zeros on the end
    Baug = vstack((sysc.B,zeros((3,1))))
    #for C matrix, output everything:
    Caug = eye(7)
    #for D matrix, use 0
    Daug = 0
    sysa = control.StateSpace(Aaug,Baug,Caug,Daug)
    #now discretize the system
    import control.matlab as cnt
    sysd = cnt.c2d(sysa,dT)
    nstates = sysd.A.shape[0]
    #now augment the system with the road states.
    Aprev = hstack((sysd.A,zeros((nstates,Np))))
    Aprev = vstack((Aprev,zeros((Np,Np+nstates))))
    #fill in the identity portion of the preview controller:
    #identity matrix should be of dimention Np-1, and begin at 8th column, 7th row.
    Aprev[nstates:-1,(nstates+1):] = eye(Np-1)

    B1prev = vstack((sysd.B,zeros((Np,1))))
    B2prev = zeros((Np+nstates,1))
    B2prev[-1,0] = 1 #this is the lone 1 for road input.
    Bprev = hstack((B1prev,B2prev)) #captures both road input and torque input
    Cprev = zeros((nstates,Np+nstates))
    Cprev[0:nstates,0:nstates] = eye(nstates)
    Dprev = 0
    sysprev = control.StateSpace(Aprev,Bprev,Cprev,Dprev,dT)

    return sysprev,Aprev,Bprev,Cprev,Dprev

def getPreviewControlLQRRazorRollLaneI(velocity, dT, Np):
    #start with the open-loop system
    sys,A0,B0,C0,D0 = getStateSpaceForPreviewControlRollLaneI(velocity,dT,Np)
    #set up the cost function
    Ccost = zeros((4,Np+7))
    #cost should weight roll error NOW, or roll-roll0
    #reminder: state vector is [roll steer droll dsteer psi y roll0 roll1 ... rollnp-1]
    # Ccost[:,0:8] = array([[0,0,0,0,0,0,1,0],[0,0,0,0,1,0,-1/(velocity*dT),1/(velocity*dT)]])
    # Ccost[:,0:8] = array([[-1,0,0,0,0,0,1,0],[0,velocity**2/(w*g),0,0,0,0,1,0]])
    Ccost[:,0:9] = array([[-1,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,1,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0]])


    weights = array([[3,0,0,0],[0,1,0,0],[0,0,.5,0],[0,0,0,1]])
    Q = np.dot(Ccost.T,weights)
    Q = np.dot(Q,Ccost)
    #cost on steering effort
    R = 1
    #now, find the LQR for non-preview portion of the system.
    #need to pull portion of original B that is 8x1 and deals with steer torque
    Bprev = vstack(sys.B[:,0])
    Klqr,Slqr,Elqr = control.dlqr(sys.A,Bprev,Q,R)
    return sys,Klqr

def getLaneErrorLQRRazor(velocity):
    """
    sys,Klqr = getRollLQRRazor(velocity)
    Using hard-coded weights, design an LQR for motorcycle.
    The final CL system will take two inputs: desired roll angle and desired lateral position
    """
    #get open-loop system
    sys,Ass,Bss,Css,Dss = getStateSpaceForControlWithRoadCoords(velocity)
    #set up cost function
    # Q = zeros((8,8))
    Q = eye(8)/100.
    Q[7,7] = 1
    R = .1

    Klqr,Slqr,Elqr = control.lqr(sys,Q,R)
    print("######### LQR GAINS for Lateral Position control #############")
    print(Klqr)
    return sys,Klqr


def getRollLQRRazor(velocity,R=0.01):
    """
    sys,Klqr = getRollLQRRazor(velocity)
    Using hard-coded weights, design a linear quadratic regulator for motorcycle.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired roll angle as the only input.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlWithIntegrator(velocity)
    Q = eye(5)/100.0
    Q[1,1] =0
    Q[4,4] =1
    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.lqr(sys, Q, R)
    print("######### LQR GAINS for ROLL control #############")
    print(Klqr)
    return sys,Klqr

def getRollDLQRRazor(velocity,R=0.01,dT=0.02):
    """
    sys,Klqr = getRollLQRRazor(velocity)
    Using hard-coded weights, design a linear quadratic regulator for motorcycle.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired roll angle as the only input.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlWithIntegrator(velocity)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(5)/100.0
    Q[1,1] =0
    Q[4,4] =1
    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    print("######### LQR GAINS for ROLL control #############")
    print(Klqr)
    return sys,Klqr

def getRollDLQRRazorNoIntegrator(velocity,R=0.01,dT=0.02):
    """
    sys,Klqr = getRollLQRRazor(velocity)
    Using hard-coded weights, design a linear quadratic regulator for motorcycle.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired roll angle as the only input.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlNoIntegrator(velocity)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(4)/100.0
    Q[0,0] = 1.0
    Q[1,1] =1.0
    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    print("######### LQR GAINS for ROLL control #############")
    print(Klqr)
    return sys,Klqr

def getRollDLQRRazorActuator(velocity,R=0.01,dT=0.02,tau_motor = .025):
    """
    sys,Klqr = getRollLQRRazorActuator(velocity,R,dT,tau_motor)
    R is LQR weight on control effort
    dT is control loop timestep
    tau_motor is estimated time constant of motor dynamics from commanded to actual torque.
    Using hard-coded weights, design a linear quadratic regulator for motorcycle with actuator dynamics.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired roll angle as the only input.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlWithActuator(velocity,tau_motor)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(5)/100.0
    Q[0,0] = 1
    Q[1,1] =0
    Q[4,4] =0
    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    print("######### LQR GAINS for ROLL control #############")
    print(Klqr)
    return sys,Klqr

def getRollDLQRRazorActuatorMotorInertia(velocity,R=0.01,dT=0.02,tau_motor = .025):
    """
    sys,Klqr = getRollLQRRazorActuator(velocity,R,dT,tau_motor)
    R is LQR weight on control effort
    dT is control loop timestep
    tau_motor is estimated time constant of motor dynamics from commanded to actual torque.
    Using hard-coded weights, design a linear quadratic regulator for motorcycle with actuator dynamics.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired roll angle as the only input.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlWithActuatorMotorInertia(velocity,tau_motor)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(5)/100.0
    Q[0,0] = 1
    Q[1,1] =0
    Q[4,4] =0
    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    print("######### LQR GAINS for ROLL control #############")
    print(Klqr)
    return sys,Klqr

def getRollDLQRRazorActuatorLane(velocity,R=0.01,dT=0.02,tau_motor = .025):
    """
    sys,Klqr = getRollLQRRazorActuatorLane(velocity,R,dT,tau_motor)
    R is LQR weight on control effort
    dT is control loop timestep
    tau_motor is estimated time constant of motor dynamics from commanded to actual torque.
    Using hard-coded weights, design a linear quadratic regulator for motorcycle with actuator dynamics.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired roll angle as the only input.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlWithActuatorLane(velocity,tau_motor)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(7)/100.0
    Q[0,0] = 1
    Q[1,1] =0
    Q[4,4] =0
    Q[6,6] = 1
    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    print("######### LQR GAINS for ROLL control #############")
    print(Klqr)
    return sys,Klqr

def getClosedLoopRazorForLaneError(sys,Klqr):
    """
        syscl, eigs = getClosedLoopRazor(sys,Klqr)

        Given a set of state feedback gains for the lane-error system, construct the final closed loop
        state space model. Consider as an input ONLY the desired lateral motorcycle position.

    """

    Acl = sys.A-dot(sys.B,Klqr)
    Bcl = zeros((8,1))
    Bcl[7,0] = 1
    Ccl = vstack((eye(8),Klqr))
    Dcl = 0
    print("####### CLOSED LOOP EIGS for LANE POSITION control: ######")
    eigs,vecs = linalg.eig(Acl)
    print(eigs)
    syscl = control.ss(Acl,Bcl,Ccl,Dcl)
    return syscl, eigs


def getClosedLoopRazor(sys,Klqr):
    """
    syscl, eigs = getClosedLoopRazor(sys,Klqr)

    Given a set of state feedback gains and a matrix of gains (presumably from LQR,
    but could be any gains), return the closed loop system. This function can be called
    to do an eigenvalue study by designing ONE set of gains and testing them at different speeds

    """
    #create a closed loop feedback system with this state feedback CONTROLLER
    Acl = sys.A - dot(sys.B,Klqr)
    # Bcl = dot(sys.B,Klqr)
    #only keep the first column of this, since we only care about roll as an input
    # Bcl = Bcl[:,0]

    Bcl = array([[0],[0],[0],[0],[1]])


    Ccl = vstack((eye(5),Klqr))
    Dcl = 0

    print("####### CLOSED LOOP DYNAMICS [roll steer rollrate steerrate int_error] #########")
    print("A: "+str(Acl))
    print("B: "+str(Bcl))
    print("C: "+str(Ccl))
    print("D: "+str(Dcl))

    print("####### CLOSED LOOP EIGS: for ROLL control ######")
    eigs,vecs = linalg.eig(Acl)
    print(eigs)
    syscl = control.ss(Acl,Bcl,Ccl,Dcl)
    return syscl, eigs

def designClosedLoopRazor(velocity):
    """
    syscl, Klqr = designClosedLoopRazor(velocity)
    This is a wrapper that builds and designs a closed-loop LQR for the motorcycle.
    The input to the closed loop system is a desired roll angle

    """
    sys,Klqr = getRollLQRRazor(velocity)
    syscl,eigs = getClosedLoopRazor(sys,Klqr)
    return syscl,Klqr

def designClosedLoopRazorLaneError(velocity):
    """
    syscl, Klqr = designClosedLoopRazorLaneError(velocity)
    This is a wrapper that builds and designs a closed-loop LQR for the motorcycle.
    The input to the closed loop system is a desired lateral (left-right) position in the lane
    """

    sys,Klqr = getLaneErrorLQRRazor(velocity)
    syscl,eigs = getClosedLoopRazorForLaneError(sys,Klqr)
    return syscl,Klqr

def designPreviewController(velocity, Np,dT):
    """
    syscl,Klqr = designPreviewController(velocity, Np,dT)
    Np is the number of preview points.
    dT is the discretization timestep of the preview controller.
    """
    pass

def getStateSpaceForYawControl(velocity):
    """ sys, A,B,C,D = getStateSpaceForYawControl(velocity)
            Creates open-loop state space model for motorcycle.
            States are x = [roll steer roll_rate steer_rate yaw]'
        """
    # first get the baseline model
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForControlWithActuator(velocity)
    # for new A matrix, pad with zeros on right
    Aaug = hstack((sys.A, zeros((5, 1))))
    # and add the terms for yaw rate on the bottom
    Aaug = vstack((Aaug, zeros((1, 6))))
    # fill in the two entries
    Aaug[5,1] = v*sin(lam)/w
    Aaug[5,3] = t*sin(lam)/w
    # now augment B matrix with an empty row
    Baug = vstack((sys.B,0))
    # output everything, I assume?
    Caug = eye(6)
    # for D matrix, use 0
    Daug = 0
    sys = control.StateSpace(Aaug, Baug, Caug, Daug)
    return sys, Aaug, Baug, Caug, Daug

def getYawDLQRRazor (velocity,R=0.01,dT=0.005):
    """
    sys,Klqr = getYawDLQRRazor(velocity)
    Using hard-coded weights, design a linear quadratic regulator for motorcycle.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired yaw angle and roll angle as the inputs.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForYawControl(velocity)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(6)/100.0
    Q[0,0] = 1.0
    Q[1,1] = 0.0
    Q[4,4] = 0.0
    #give yaw some beans but maybe turn back down a bit if things get too torquey.
    Q[5,5] = 2.0

    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    #print("######### LQR GAINS for YAW control #############")
    #print(Klqr)
    return sys,Klqr

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
    Aaug[4,1] = v*sin(lam)/w
    Aaug[4,3] = t*sin(lam)/w
    # now augment B matrix with an empty row
    Baug = vstack((sys.B,0))
    # output everything, I assume?
    Caug = eye(5)
    # for D matrix, use 0
    Daug = 0
    sys = control.StateSpace(Aaug, Baug, Caug, Daug)
    return sys, Aaug, Baug, Caug, Daug

def getYawDLQRRazornoactuator (velocity,R=0.01,dT=0.005):
    """
    sys,Klqr = getYawDLQRRazor(velocity)
    Using hard-coded weights, design a linear quadratic regulator for motorcycle.
    returns open-loop system (as a control object) and the matrix of LQR gains.
    This LQR takes desired yaw angle and roll angle as the inputs. No actuator compensation.
    """
    sys,Ass,Bss_control,Css,Dss_control = getStateSpaceForYawControlnoactuator(velocity)
    import control.matlab as cnt
    sysd = cnt.c2d(sys,dT)
    Q = eye(5)/100.0
    Q[0,0] = 1.0
    Q[1,1] = 0.0
    #give yaw some beans but maybe turn back down a bit if things get too torquey.
    Q[4,4] = 1.0

    #set up LQR weight on input matrix, which is
    #R = array([[1,0],[0,.01]])
    # R = .01
    Klqr, Slqr, Elqr = control.dlqr(sysd, Q, R)
    #print("######### LQR GAINS for YAW control #############")
    #print(Klqr)
    return sys,Klqr
#this is now a function that will run only if we use this file as a script.
def main():
    sys,Ass,Bss,Css,Dss = getStateSpace(5)

    step_torque_mag = 0.1
    tsim,yout = control.step_response(sys)
    yout = yout*step_torque_mag
    figure()
    subplot(2,1,1)
    plot(tsim,yout[0][1][:],'k')
    xlabel('time (s)')
    ylabel('roll angle (rad)')
    title('Open-Loop Step response from steer torque to roll angle')
    subplot(2,1,2)
    plot(tsim,yout[1][1][:],'k')
    xlabel('time (s)')
    ylabel('steer angle (rad)')
    title('Open-Loop Step response from steer torque to steer angle')


    #now set up loop to calculate eigenvalues at various speeds
    vvec = arange(.01,10,.01)
    vvec2 = zeros((4,len(vvec)))
    eigs_re = zeros((4,len(vvec)))
    eigs_im = zeros((4,len(vvec)))





    ###### eigenvalue study for the open loop bike:

    for k in range(0,len(vvec)):
        v = vvec[k]
        K = K0+K1*v**2
        D = v*D1

        #now translate the model into state space
        A = hstack((zeros((2,2)),eye(2,2)))
        A = vstack((A,hstack((dot(-linalg.inv(M),K),dot(-linalg.inv(M),D)))))
        eigs,vecs = linalg.eig(A)
        eigs = sort(eigs)
        eigs_re[:,k] = real(eigs)
        eigs_im[:,k] = imag(eigs)
        vvec2[:,k] = [v,v,v,v]

    figure()
    plot(vvec,eigs_re[0,:],'k',vvec,abs(eigs_im[0,:]),'k-.')
    xlabel('Speed (m/s)')
    ylabel('Eig (1/s)')
    legend(['real','imag'])
    plot(vvec,eigs_re[1,:],'k',vvec,abs(eigs_im[1,:]),'k-.')
    plot(vvec,eigs_re[2,:],'k',vvec,abs(eigs_im[2,:]),'k-.')
    plot(vvec,eigs_re[3,:],'k',vvec,abs(eigs_im[3,:]),'k-.')


    axis([0,10,-10,10])

    # https://python-control.readthedocs.io/en/latest/generated/control.lqr.html#control.lqr
    velocity = 4#6.25
    sys,Ass,Bss,Css,Dss = getStateSpaceForControlWithIntegrator(velocity)
    #set up LQR weight matrix on each state
    # Q = array([[1,0,0,0,0],[0,1,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,100]])

    #use the function above to get a closed loop model of the razor
    #I put this in a function so we can use this from another file.
    syscl,Klqr = designClosedLoopRazor(velocity)

    #do a simulation of the closed loop behavior of the system under a step in desired roll

    goalAccel = -0.1*g #goal lateral acceleration is U^2/R
    goalRadius = velocity**2/goalAccel
    goalYawRate = velocity/goalRadius

    #now what does that mean for roll?
    goalRoll = -arctan(velocity*goalYawRate/9.81)
    # #what does it mean for goal steering? delta = L/R, ignoring trig stuff.
    # goalSteer = w/goalRadius
    tsim = linspace(0,5,1000)

    xdesired = zeros((len(tsim),1))

    xdesired[:,0] = goalRoll


    import control.matlab as cnt
    ycl,tsim_out,xcl = cnt.lsim(syscl,xdesired,tsim)
    #compute steer torque

    figure()

    subplot(3,1,1)
    title("Closed Loop Step Response: Desired Roll = "+"{:.2f}".format(goalRoll*180/pi)+" degrees")
    plot(tsim,goalRoll*ones((len(tsim),1)),'k--',tsim,ycl[:,0],'k')
    xlabel('Time (s)')
    ylabel('Roll Angle (rad)')
    legend(['desired','actual'])
    subplot(3,1,2)
    plot(tsim,ycl[:,1],'k')
    ylabel('Steer Angle (rad)')
    subplot(3,1,3)
    plot(tsim,ycl[:,5],'k')
    xlabel('Time (s)')
    ylabel('Steer Torque (Nm)')



    ########## DO A DESIGN FOR LANE ERROR TOO!
    syscl_lane, Klqr_lane = designClosedLoopRazorLaneError(velocity)
    goalPos = 1.0 #do a 1m lane change
    tsim = linspace(0,5,1000)
    xdesired_lane = zeros((len(tsim),1))
    xdesired_lane[:,0] = goalPos
    ycl,tsim_out,xcl = cnt.lsim(syscl_lane,xdesired_lane,tsim)
    #compute steer torque
    figure()

    subplot(2,2,1)
    title("Closed Loop Step Response (lane error)")
    plot(tsim,goalPos*ones((len(tsim),1)),'k--',tsim,ycl[:,6],'k')
    xlabel('Time (s)')
    ylabel('lateral position (m)')
    legend(['desired','actual'])
    subplot(2,2,2)
    plot(tsim,ycl[:,1],'k')
    xlabel('Time (s)')
    ylabel('Steer Angle (rad)')
    subplot(2,2,3)
    plot(tsim,ycl[:,0],'k')
    ylabel('Roll Angle (rad)')
    xlabel('Time (s)')
    subplot(2,2,4)
    plot(tsim,ycl[:,8],'k')
    xlabel('Time (s)')
    ylabel('Steer Torque (Nm)')


    show()




if __name__ == '__main__':
    pass
    #main()
