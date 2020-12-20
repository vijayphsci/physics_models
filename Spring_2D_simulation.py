#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Spring_2D:
    def __init__(self,mass=1,stiffness=1,damping=0,gravity=2,time=200,dt=0.1,natural_length=4,xinitial=2,yinitial=-2,velocity_x=0,velocity_y=0,frame_interval=15,figsize=(10,8),repeat=True):
        self.m=mass
        self.k=stiffness
        self.damp=damping
        self.gravity=gravity
        self.xi=xinitial
        self.yi=yinitial
        self.vxi=velocity_x
        self.vyi=velocity_y
        self.l0=natural_length
        self.tf=time
        self.dt=dt
        self.x=[]
        self.y=[]
        self.vx=[]
        self.vy=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.frame_interval=frame_interval
        self.figsize=figsize
        self.count=0
        self.__ani=None
        self.repeat=repeat
    def f(self,x,y,vx):
        return (self.k*x*(self.l0/np.sqrt(x**2+y**2)-1)-self.damp*vx)/self.m
    def g(self,x,y,vy):
        return (self.k*y*(self.l0/np.sqrt(x**2+y**2)-1)-self.damp*vy-self.gravity)/self.m
    def __solve_system(self):
        self.t.append(0)
        self.y.append(self.yi)
        self.x.append(self.xi)
        self.vy.append(self.vyi)
        self.vx.append(self.vxi)
        for i in range(self.nt):
            t1=self.f(self.x[i],self.y[i],self.vx[i])
            t2=self.g(self.x[i],self.y[i],self.vy[i])
            self.t.append((i+1)*self.dt)
            self.x.append(self.x[i]+self.dt*self.vx[i]+self.dt**2/2*t1)
            self.y.append(self.y[i]+self.dt*self.vy[i]+self.dt**2/2*t2)
            self.vx.append(self.vx[i]+self.dt*self.f(self.x[i]+self.dt/2*self.vx[i],self.y[i]+self.dt/2*self.vy[i],self.vx[i]+self.dt/2*t1))
            self.vy.append(self.vy[i]+self.dt*self.g(self.x[i]+self.dt/2*self.vx[i],self.y[i]+self.dt/2*self.vy[i],self.vy[i]+self.dt/2*t2))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min(self.x)-2,max(self.x)+2), ylim=(min(self.y)-2,max(self.y)+4)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        cline1,=ax.plot([], [], lw=2,c='r')
        tline1,=ax.plot([], [], lw=1,c='m')
        plt.hlines(0,-2,2)
        tj1x=[]
        tj1y=[]
        def animate(i):
            tj1x.append(self.x[i])
            tj1y.append(self.y[i])
            cx1,cy1=circle(self.x[i],self.y[i],0.2)
            #tn1,tn2=sline(0,0,self.x[i],self.y[i])
            tn1,tn2=springcoil(0,0,self.x[i],self.y[i],10,0.1,60)
            if self.x[i]**2+self.y[i]**2<self.l0**2:
                pline1.set_color('lime')
            else:
                pline1.set_color('b')
            pline1.set_data(tn1,tn2)
            cline1.set_data(cx1,cy1)
            tline1.set_data(tj1x,tj1y)
            return cline1,pline1,tline1
        self.__ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show() 
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min(self.x)-2,max(self.x)+2), ylim=(min(self.y)-2,1)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        cline1,=ax.plot([], [], lw=2,c='r')
        tline1,=ax.plot([], [], lw=1,c='m')
        plt.hlines(0,-2,2)
        def animate(i):
            cx1,cy1=circle(self.x[i],self.y[i],0.2)
            #tn1,tn2=sline(0,0,self.x[i],self.y[i])
            tn1,tn2=springcoil(0,0,self.x[i],self.y[i],10,0.1,60)
            if self.x[i]**2+self.y[i]**2<self.l0**2:
                pline1.set_color('lime')
            else:
                pline1.set_color('b')
            pline1.set_data(tn1,tn2)
            cline1.set_data(cx1,cy1)
            return cline1,pline1
        self.__ani= animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show() 
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x,self.y,c='m',label='spring')
        plt.legend()
        plt.show()
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig, axs = plt.subplots(2,figsize=(12,12))
        axs[0].plot(self.t,self.x,c='r')
        axs[1].plot(self.t,self.y,c='b')
        axs[0].set_title('x coordinate')
        axs[1].set_title('y coodinate')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='position')