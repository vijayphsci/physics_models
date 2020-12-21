#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Inverse_Squared_Law:
    def __init__(self,mass=1,attraction_const=-13,position=(5,0),velocity=(0,1),time=200,dt=0.1,frame_interval=15,figsize=(10,8),repeat=True):
        self.m=mass
        self.k=attraction_const
        self.ini_pos=position
        self.ini_vel=velocity
        self.frame_interval=frame_interval
        self.figsize=figsize
        self.tf=time
        self.dt=dt
        self.x=[]
        self.y=[]
        self.vx=[]
        self.vy=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.count=0
        self.ani=None
        self.repeat=repeat
    def fx(self,x,y):
        return self.k*x/(self.m*np.power(x**2+y**2,3/2))
    def fy(self,x,y):
        return self.k*y/(self.m*np.power(x**2+y**2,3/2))
    def __solve_system(self):
        self.t.append(0)
        self.x.append(self.ini_pos[0])
        self.y.append(self.ini_pos[1])
        self.vx.append(self.ini_vel[0])
        self.vy.append(self.ini_vel[1])
        for i in range(self.nt):
            t1=self.fx(self.x[i],self.y[i])
            t2=self.fy(self.x[i],self.y[i])
            self.t.append((i+1)*self.dt)
            self.x.append(self.x[i]+self.dt*self.vx[i]+self.dt**2/2*t1)
            self.y.append(self.y[i]+self.dt*self.vy[i]+self.dt**2/2*t2)
            self.vx.append(self.vx[i]+self.dt*self.fx(self.x[i]+self.dt/2*self.vx[i],self.y[i]+self.dt/2*self.vy[i]))
            self.vy.append(self.vy[i]+self.dt*self.fy(self.x[i]+self.dt/2*self.vx[i],self.y[i]+self.dt/2*self.vy[i]))
            
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        xlim_min=min(self.x)-2
        xlim_max=max(self.x)+2
        ylim_min=min(self.y)-2
        ylim_max=max(self.y)+2
        ax = plt.axes(xlim=(xlim_min,xlim_max), ylim=(ylim_min,ylim_max)) 
        pline1,=ax.plot([], [], lw=3,c='k',marker='o') 
        tline1,=ax.plot([], [], lw=1,c='m')
        plt.hlines(0,xlim_min,xlim_max,linestyle='--',color='r')
        plt.vlines(0,ylim_min,ylim_max,linestyle='--',color='r')
        def animate(i):
            pline1.set_data(self.x[i],self.y[i])
            tline1.set_data(self.x[:i],self.y[:i])
            return pline1,tline1
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show() 
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        xlim_min=min(self.x)-2
        xlim_max=max(self.x)+2
        ylim_min=min(self.y)-2
        ylim_max=max(self.y)+2
        ax = plt.axes(xlim=(xlim_min,xlim_max), ylim=(ylim_min,ylim_max)) 
        pline1,=ax.plot([], [], lw=3,c='k',marker='o') 
        plt.hlines(0,xlim_min,xlim_max,linestyle='--',color='r')
        plt.vlines(0,ylim_min,ylim_max,linestyle='--',color='r')
        def animate(i):
            pline1.set_data(self.x[i],self.y[i])
            return pline1,
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show() 
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        xlim_min=min(self.x)-2
        xlim_max=max(self.x)+2
        ylim_min=min(self.y)-2
        ylim_max=max(self.y)+2
        ax = plt.axes(xlim=(xlim_min,xlim_max), ylim=(ylim_min,ylim_max)) 
        plt.hlines(0,xlim_min,xlim_max,linestyle='--',color='b')
        plt.vlines(0,ylim_min,ylim_max,linestyle='--',color='b')
        plt.plot(self.x,self.y,c='r')
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