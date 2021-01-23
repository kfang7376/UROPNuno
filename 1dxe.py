#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 23:06:29 2021

@author: kellyfang
"""

#1/25 allow arbitrary # of electrons

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
from random import randrange
#import matplotlib.animation as animation

#constants
eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C

#create electrons
class electron():
    
    def __init__(self,x,v,f):   #added new attribute that tells what the net force on the electron is
        self.x=x
        self.v=v
        self.f=f
        self.poslist=np.zeros(n)
        self.vellist=np.zeros(n)

       
    def ChangePos(self,xnew,i):     #sets the new current position of the electron
        self.x=xnew
        self.AddPos(xnew,i)
    
    def ChangeVel(self,vnew,i):     #sets the new current vel of the electron
        self.v=vnew
        self.AddVel(vnew,i)
    
    def ChangeFor(self,fnew):
        self.f=fnew
 
    def getFor(self):
        return self.f
    
    def getPos(self):       #return the electron's position
        return self.x
    
    def getVel(self):       #returns the electron's velocity
        return self.v
    
    
    def AddPos(self,value,i):       #adds positions to the position array
       self.poslist[i]=value
       
       
    def AddVel(self,value,i):       #adds velocities to the velocity array
       self.vellist[i]=value
       
       
    def getPosListi(self,i):        #returns a sepecific index from the position array
        return self.poslist[i]
    
    def getVelListi(self,i):         #returns a sepecific index from the velocity array
        return self.vellist[i]
    
    def getPosList(self):       #returns the position array
        return self.poslist
    
    def getVelList(self):       #returns the velocity array
        return self.vellist
    
    def getPosVel(self):        #returns the current position and velocity as an array
        return [self.getPos(),self.getVel()]
    
    def CheckPos(self,pos2):        #checks if the electrons is to the left or right of another position
        if(self.x-pos2>0):
            return 1
        else:
            return 0
        

def getDistance(d1,d2):       #finds the distance from some other position
    return (d1-d2)

#sets the timepoints
n=10000
t=np.linspace(0,5,n)

initpos = [5,10,-12,-3,14]
    
num_e = len(initpos) #number of electrons in system

electrons = [] #list of electrons

for i in range (num_e):  #creates num_e number of electrons w/ rand pos and 0 vel
    x0 = initpos[i]
    v0 = 0
    f0 = 0
    electrons.append(electron(x0,v0,f0)) #adds new instance of electron to an array
    print(x0)
    

def model1(z,t):
        v=z[1] #sets v equal to 2nd imputed array position
        dvdt= electron.getFor()/mass
        dxdt=v
        return[dxdt,dvdt]

for x in range (0,n):
    pos_now = []
    for electron in electrons:
        pos_now.append(electron.getPos())
    
    for electron in electrons:
        fnet = 0
        distances = np.zeros(num_e)
        for i in range (num_e):
            distances[i] = getDistance(electron.getPos(),pos_now[i])
        #print(distances)
        
        for i in range (num_e):
            if (distances[i]==0):
                force = 0
            if (distances[i]>0):
                force = k*(elecCharge**2)/((distances[i]**2)) #force will push electron right
            if (distances[i]<0):
                force = -k*(elecCharge**2)/((distances[i]**2)) #force will push electron left
            
            fnet = fnet + force
                
        electron.ChangeFor(fnet)
    
    
    for electron in electrons:
        z = odeint(model1,electron.getPosVel(),t)
        z0=z[1]
        electron.ChangePos(z0[0],x)
        electron.ChangeVel(z0[1],x)
 
for electron in electrons:
    plt.plot(t,electron.getPosList())
    plt.xlabel('time')
    plt.ylabel('position')

plt.show()
    
for electron in electrons:
    plt.plot(t,electron.getVelList())
    plt.xlabel('time')
    plt.ylabel('velocity')

plt.show()
        

    

    

    
