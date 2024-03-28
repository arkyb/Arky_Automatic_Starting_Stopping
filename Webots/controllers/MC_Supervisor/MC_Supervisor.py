from controller import Supervisor
from controller import Node
import os
import datetime
import csv
import numpy as np

robot = Supervisor()  # create Supervisor instance

#open Ksum and Delay files, set initial values
f = open('../Yaw_Controller_ISO/ksum.txt','w')
g = open('../Yaw_Controller_ISO/delay.txt','w')
f.write('1.0')
g.write('0.005')
g.close()
f.close()
success = np.ones([18,7])

# get base nodes
mc_node = robot.getFromDef('MOTORCYCLE')
st_node = robot.getFromDef('steer_joint')
world_node = robot.getFromDef('circle_arena')

#get fields
trans_field = mc_node.getField('translation')
rot_field = mc_node.getField('rotation')
fric_field = st_node.getField('staticFriction')
world_rot_field = world_node.getField('rotation')

#get values from fields
trans_init = trans_field.getSFVec3f()
rot_init = rot_field.getSFRotation()

#get current time
now = datetime.datetime.now()

#create first video name. uncomment if recording videos
#fname = 'Videos/'+str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'_'+str(now.hour)+'_'+str(now.minute)+'_'+str(now.second)+'.mp4'

#initialize indices, and intial values
i = 0
j = 0
k=1
m=1
delay = 0.005
initVel = 3
TIME_STEP = 5

print('Iteration number '+str((k-1)*7+m)+': Ksum = '+str(1-0.05*(k-1))+', delay = '+str(delay))
while robot.step(TIME_STEP) != -1:
    rot = rot_field.getSFRotation()
    roll = rot[3]*np.dot([1,0,0],[rot[0],rot[1],rot[2]])
    #fail condition
    if 1 < roll or roll < -1:
        success[k-1,m-1]=0
        print('FAIL')
    if (i<1):
        #uncomment if recording videos
        #robot.movieStartRecording(fname,1280,720,1,40,1,1>2)
        rot_field.setSFRotation([0.05,0.02,0,-0.25])
    if (i < 10):
        mc_node.setVelocity([initVel,0,0,0,0,0])
    i += 1
    #uncomment if recording videos
    #if (robot.getTime()>95 and j==0):
        #j+=1
        #robot.movieStopRecording()
    if (robot.getTime()>12):
        #uncomment if recording videos 
        #now = datetime.datetime.now()
        #fname = 'Videos/'+str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'_'+str(now.hour)+'_'+str(now.minute)+'_'+str(now.second)+'.mp4'
        if m==7:
           k += 1
           m=0
        m+=1
        if m==1:
            delay = 0.005
        elif m==2:
            delay = 0.05
        elif m==3:
            delay = 0.06
        elif m==4:
            delay = 0.07
        elif m==5:
            delay = 0.08
        elif m==6:
            delay= 0.09
        elif m==7:
            delay = 0.1
        g = open('../Yaw_Controller_ISO/delay.txt','w')
        f = open('../Yaw_Controller_ISO/ksum.txt','w')
        print('Iteration number '+str((k-1)*7+m)+': Ksum = '+str(1-0.05*(k-1))+', delay = '+str(delay))
        f.write(str(1-0.05*(k-1)))
        g.write(str(delay))
        with open('../Yaw_Controller_ISO/success_tab.csv',mode='w') as success_file:
            success_writer = csv.writer(success_file, delimiter=',')
            success_writer.writerow(success)
        f.close()
        g.close()
        robot.simulationReset()
        mc_node.restartController()
        mc_node.resetPhysics()
        trans_field.setSFVec3f(trans_init)
        i=0
        j=0
