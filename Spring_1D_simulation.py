#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Spring_1D:
    def __init__(self,mass=1,stiffness=1,damping=0.03,gravity=2,time=200,dt=0.1,natural_length=4,yinitial=-2,velocity=0,frame_interval=15,figsize=(10,8),repeat=True):
        self.m=mass
        self.k=stiffness
        self.damp=damping
        self.gravity=gravity
        self.yi=yinitial
        self.vyi=velocity
        self.l0=natural_length
        self.tf=time
        self.dt=dt
        self.y=[]
        self.vy=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.frame_interval=frame_interval
        self.figsize=figsize
        self.count=0
        self.__ani=None
        self.repeat=repeat
    def f(self,y,vy):
        return (-self.k*(self.l0+y)-self.gravity-self.damp*vy)/self.m
    def __solve_system(self):
        self.t.append(0)
        self.y.append(self.yi)
        self.vy.append(self.vyi)
        for i in range(self.nt):
            t1=self.f(self.y[i],self.vy[i])
            self.t.append((i+1)*self.dt)
            self.y.append(self.y[i]+self.dt*self.vy[i]+self.dt**2/2*t1)
            self.vy.append(self.vy[i]+self.dt*self.f(self.y[i]+self.dt/2*self.vy[i],self.vy[i]+self.dt/2*t1))
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-4,4), ylim=(min(self.y)-3,1)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        cline1,=ax.plot([], [], lw=2,c='r')
        plt.hlines(0,-2,2)
        def animate(i):
            cx1,cy1=circle(0,self.y[i],0.2)
            tn1,tn2=springcoil(0,0,0,self.y[i],10,0.1,70)
            if self.y[i]**2<self.l0**2:
                pline1.set_color('lime')
            else:
                pline1.set_color('b')
            pline1.set_data(tn1,tn2)
            cline1.set_data(cx1,cy1)
            return cline1,pline1
        self.__ani = animation.FuncAnimation(fig,animate,frames=self.nt+1,interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()
    def show_phase_space_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        plt.title('Phase Space Plot')
        plt.xlabel('Position')
        plt.ylabel('Momentum')
        plt.plot(self.y,self.vy,c='r')
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig = plt.figure(figsize=self.figsize) 
        plt.title('Position time plot')
        plt.xlabel('time')
        plt.ylabel('position')
        plt.plot(self.t,self.y,c='r')
        plt.show()