from controller import Robot, Motor, InertialUnit, Display
from numpy import *
from Rollover import Rollover
from DelaySim import DelaySim
from LCAD import LCAD
from generateLQRvals import generateLQRvals
import time
import math
import datetime

now = datetime.datetime.now()
fname = 'Data/'+str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'_'+str(now.hour)+'_'+str(now.minute)+'_'+str(now.second)+'.txt'
dtxt = open('delay.txt')
delayTime = float(dtxt.read())
ktxt = open('ksum.txt')
ksum = float(ktxt.read())
recordData = True
if recordData:
    #start a file we can use to collect data
    f = open(fname,'w')
    f.write("#time, goalRoll, Torque, speed, roll, steer, rollrate, steerrate, intE, fsmstate, kickangle, Tdelayed\r\n")

class Timer:
    def __init__(self,preset):
        #current "state" of the timer
        self.state = False
        #current elapsed time
        self.elapsed = 0
        #timer will go true if it has counted for more than 1 second
        self.preset = preset
    def update(self,ENABLE,dt):
        #dt is the timestep by which we should count up. make sure its units match your preset!
        #don't set the preset in seconds and increment the elapsed time in milliseconds, for example!
        #ENABLE is a boolean. When it is true, we run up the timer. When it is not, the time resets and we stop counting.
        if(ENABLE):
            #increment time by dt
            self.elapsed+=dt
            self.state=self.elapsed>=self.preset
        else:
            self.elapsed=0
            self.state=False

#create the Robot instance.
robot = Robot()
yawCorr = Rollover()

#get the time step of the current world.
timestep = int(robot.getBasicTimeStep())

#CONTROL PARAMETERS
lastControlTime = 0
dTcontrol = 0.005
delayTime = 0.005

#get LQR gains
generator = generateLQRvals(dTcontrol)
gainmtx=generator.update()
lqrspeeds = gainmtx[:,0]
lqrgains = gainmtx[:,1:]

Ndelay = int(delayTime/(dTcontrol)) #length of delay queue
tDelay = DelaySim(Ndelay) #create DelaySim object

lcad = LCAD() #create an LCAD object

#SIMULATION SETUP
driveVelocity = 3.0
stopRollVal = -0.05
T = 0
fricComp = 0
simtime = 0.0
yawRate = 0
oldYaw = 0
rollInt = 0
inteYawRate = 0
oldRoll = 0
steerangle = 0
oldsteer = 0
steerRate = 0
rollMax = 0
rollRateMax = 0

#STATE VARIABLES
driving = True
stopping = False
stopped = False
sup = False

#TIMERS
driveTimer = Timer(5000)
stopTimer = Timer(1000)
stoppedTimer = Timer(2000) 
supTimer = Timer(500)

motor = robot.getDevice('drive_motor')
steer = robot.getDevice('steering_motor')
kick = robot.getDevice('kick_motor')
motor.setPosition(float('inf'))
imu = robot.getDevice('imu')
imu.enable(timestep)
gps = robot.getDevice('gps')
gps.enable(timestep)
gyro = robot.getDevice('gyro')
gyro.enable(timestep)
steergyro = robot.getDevice('steergyro')
steergyro.enable(timestep)
steersensor = robot.getDevice('steer_angle')
steersensor.enable(timestep)
kicksensor = robot.getDevice('kick_angle')
kicksensor.enable(timestep)

oldRoll,oldPitch,oldYaw = imu.getRollPitchYaw()

firstLoop = True

#set the simulation forward speed and calculate rear wheel omega
Rrw = 0.15875
driveOmega = driveVelocity/Rrw
swingDownAcc = 1.1*driveOmega
swingUpAcc = 2*driveOmega
def clamp(val,min,max):
    assert min<max
    if(val<min):
        val=min
    elif val>max:
        val=max
    return val

def setDriveMotorTorque(self,motor,command,omega):
    kt = 0.86/13.5 #published for MY1016-C1 at 24V max torque .86 Nm, max curr 13.5A
    Vbatt = 24.0 #total voltage
    """! sets motor torque in simulation based on a physic-based model (not for user use)"""
    Vcommand = clamp(command,-Vbatt,Vbatt)
    #this uses a physics-based model to calculate torque based on motor params.
    torque = self.kt/(self.R)*(Vcommand-self.kt*omega)
    #set motor force
    motor.setTorque(torque)

def find_nearest_index(array, value):
    array = asarray(array)
    idx = (abs(array - value)).argmin()
    return idx

def getCurrentGains(speed):
    #find the speed in the gain array closest to ours
    idx = find_nearest_index(lqrspeeds,speed)
    #find the gain set at this index
    Klqr = lqrgains[idx,:]
    return Klqr

#Main loop:
#- perform simulation steps until Webots is stopping the controller
while robot.step(timestep) != -1:
    simtime+=timestep/1000.0
    if(firstLoop):
        oldRoll,oldPitch,oldYaw = imu.getRollPitchYaw()
        oldYaw = yawCorr.update(oldYaw)
        oldsteer = steersensor.getValue()
        firstLoop=False
    #BLOCK 1
    #FSM SETUP
    #Set timer states
    driveTimer.update(driving,timestep)
    stopTimer.update(stopping,timestep)
    stoppedTimer.update(stopped,timestep)
    supTimer.update(sup,timestep)
    #Read the sensors:
    #get current fwd speed
    U = gps.getSpeed();
    #get IMU values and process to get yaw, roll rates
    #read IMU value
    rpy = imu.getRollPitchYaw()
    gyros = gyro.getValues()
    yaw = rpy[2]
    yaw = yawCorr.update(yaw)
    yawRate = gyros[2]
    oldYaw = yaw
    roll = rpy[0]
    rollRate = gyros[0]
    oldRoll = roll
    #get kickstand angle
    kickangle = kicksensor.getValue()
    #now get steer values and calculate steer rate.
    #WEBOTS is in ISO (z up) but Papa is in SAE (z down) so need to flip dir.
    steerangle = steersensor.getValue()
    steergyros = steergyro.getValues()
    steerRate = steergyros[2]
    oldsteer = steerangle

    #BLOCK 2
    #SET TRANSITION BEHAVIOR
    if(rollMax>0.01):
        T1 = driving and not (driveTimer.state and rollRate<=0.5*rollRateMax and roll>=0.5*rollMax)
        T2 = driving and driveTimer.state and rollRate<=0.5*rollRateMax and roll>=0.5*rollMax
    else:
        T1 = driving and not (driveTimer.state and rollRate<=0.5*rollRateMax)
        T2 = driving and driveTimer.state and rollRate<=0.5*rollRateMax
    T3 = stopping and not U<0.1
    T4 = stopping and U<0.1
    T5 = stopped and not stoppedTimer.state
    T6 = stopped and stoppedTimer.state
    T7 = sup and not roll>= 0
    T8 = sup and roll>=0
    #BLOCK 3
    #SET THE STATES
    driving = T1 or T8
    stopping = T2 or T3
    stopped = T4 or T5
    sup = T6 or T7

    #BLOCK 4
    #SET THE STATE BEHAVIOR
    if(driving):
        state = "driving"
        goalYaw = 0
        goalRoll=0
        motor.setVelocity(driveOmega)
        kick.setPosition(1.2)
        if(driveTimer.elapsed>2000):
            rollMax,rollRateMax = lcad.update(roll,rollRate)

    elif(stopping):
        rollMax = 0
        rollRateMax = 0
        state = "stopping"
        motor.setVelocity(clamp(driveOmega - (1.1*swingDownAcc*stopTimer.elapsed/stopTimer.preset),0,driveOmega))
        kick.setPosition(1.2-clamp((stopTimer.elapsed*1.2/(500)),0,1.2))
        if(stopTimer.elapsed<0.4*stopTimer.preset):
            goalRoll = 0
            goalRollRate=-0.05
        else:
            goalRoll = stopRollVal
    elif(stopped):
        state = "stopped"
        motor.setVelocity(0)
        if(stoppedTimer.elapsed>stoppedTimer.preset*0.5):
            goalRoll = -0.12
        else:
            goalRoll = stopRollVal
    elif(sup):
        goalRoll = 0.1
        kick.setPosition(clamp(3*1.2*(supTimer.elapsed)/supTimer.preset,0,1.2))
        motor.setVelocity(clamp(swingUpAcc*supTimer.elapsed/supTimer.preset,0,driveOmega))
        state = "swingup"

    eYaw = goalYaw - yaw
    eRoll = goalRoll - roll
    rollInt = rollInt + eRoll*(dTcontrol)
    Klqr = getCurrentGains(driveVelocity)

    T = Klqr[4]*eYaw + Klqr[0]*(eRoll) - Klqr[1]*steerangle - Klqr[2]*rollRate - Klqr[3]*steerRate

    if(T<0):
        T-=fricComp
    elif(T>0):
        T+=fricComp

    Tdelayed = tDelay.update(T)

    Tfinal = Tdelayed

    Tlim = 1.5916
    if(Tfinal>Tlim):
        Tfinal = Tlim
    elif(T<-Tlim):
        Tfinal = -Tlim

    Tfinal = Tfinal*ksum

    steer.setControlPID(0.0001,0,0)
    steer.setTorque(Tfinal)
    lastControlTime = simtime
    if(recordData):
        f.write(str(simtime)+","+str(goalRoll)+","+str(T)+","+str(U)+","+str(roll)+","+str(steerangle)+","+str(rollRate)+","+str(steerRate)+","+str(rollInt)+","+str(state)+","+str(kickangle)+","+str(Tdelayed)+"\r\n")
