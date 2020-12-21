#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Double_Pendulum:
    def __init__(self,angle1=40,angle2=80,length1=3,length2=3,mass1=1,mass2=0.1,gravity=3,time=200,dt=0.1,frame_interval=15,figsize=(10,8),repeat=True):
        self.m1=mass1
        self.m2=mass2
        self.l1=length1
        self.l2=length2
        self.gravity=gravity
        self.theta1i_deg=angle1
        self.theta2i_deg=angle2
        self.theta1i_rad=self.theta1i_deg*np.pi/180
        self.theta2i_rad=self.theta2i_deg*np.pi/180
        self.tf=time
        self.dt=dt
        self.frame_interval=frame_interval
        self.theta1=[]
        self.theta2=[]
        self.dtheta1=[]
        self.dtheta2=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.x1=[]
        self.x2=[]
        self.y1=[]
        self.y2=[]
        self.figsize=figsize
        self.count=0
        self.ani=None
        self.repeat=repeat
    def f(self,th1,th2,dth1,dth2):
        a1= self.m2*self.l1*self.l2*dth2**2*mt.sin(th1-th2)+(self.m1+self.m2)*self.gravity*self.l1*mt.sin(th1)
        a2=2*self.m2*self.l1*self.l2*dth1*dth2*mt.sin(th1-th2)-self.m2*self.l1*self.l2*(dth1**2)*mt.sin(th1-th2)+self.m2*self.gravity*self.l2*mt.sin(th2)
        k1=-a1*self.m2*self.l2**2+a2*self.m2*self.l1*self.l2*mt.cos(th1-th2)
        k2=self.m1*self.m2*self.l1**2*self.l2**2+self.m2**2*self.l1**2*self.l2**2*(mt.sin(th1-th2))**2
        return k1/k2 
    def g(self,th1,th2,dth1,dth2):
        a1= self.m2*self.l1*self.l2*dth2**2*mt.sin(th1-th2)+(self.m1+self.m2)*self.gravity*self.l1*mt.sin(th1)
        a2=2*self.m2*self.l1*self.l2*dth1*dth2*mt.sin(th1-th2)-self.m2*self.l1*self.l2*(dth1**2)*mt.sin(th1-th2)+self.m2*self.gravity*self.l2*mt.sin(th2)
        k3=-a2*(self.m1+self.m2)*self.l1**2+a1*self.m2*self.l1*self.l2*mt.cos(th1-th2)
        k4=self.m1*self.m2*self.l1**2*self.l2**2+self.m2**2*self.l1**2*self.l2**2*(mt.sin(th1-th2))**2
        return k3/k4
    def __solve_system(self):
        self.theta1.append(self.theta1i_rad)
        self.theta2.append(self.theta2i_rad)
        self.dtheta1.append(0)
        self.dtheta2.append(0)
        self.t.append(0)
        self.x1.append(self.l1*mt.sin(self.theta1[0]))
        self.y1.append(-self.l1*mt.cos(self.theta1[0]))
        self.x2.append(self.l1*mt.sin(self.theta1[0])+self.l2*mt.sin(self.theta2[0]))
        self.y2.append(-self.l1*mt.cos(self.theta1[0])-self.l2*mt.cos(self.theta2[0]))
        for i in range(self.nt):
            t1=self.f(self.theta1[i],self.theta2[i],self.dtheta1[i],self.dtheta2[i])
            t2=self.g(self.theta1[i],self.theta2[i],self.dtheta1[i],self.dtheta2[i])
            self.t.append((i+1)*self.dt)
            self.theta1.append(self.theta1[i]+self.dt*self.dtheta1[i]+self.dt**2/2*t1)
            self.dtheta1.append(self.dtheta1[i]+self.dt*self.f(self.theta1[i]+self.dt/2*self.dtheta1[i],self.theta2[i]+self.dt/2*self.dtheta2[i],self.dtheta1[i]+self.dt/2*t1,self.dtheta2[i]+self.dt/2*t2))
            self.theta2.append(self.theta2[i]+self.dt*self.dtheta2[i]+self.dt**2/2*t2)
            self.dtheta2.append(self.dtheta2[i]+self.dt*self.g(self.theta1[i]+self.dt/2*self.dtheta1[i],self.theta2[i]+self.dt/2*self.dtheta2[i],self.dtheta1[i]+self.dt/2*t1,self.dtheta2[i]+self.dt/2*t2))
            self.x1.append(self.l1*mt.sin(self.theta1[i]))
            self.y1.append(-self.l1*mt.cos(self.theta1[i]))
            self.x2.append(self.l1*mt.sin(self.theta1[i])+self.l2*mt.sin(self.theta2[i]))
            self.y2.append(-self.l1*mt.cos(self.theta1[i])-self.l2*mt.cos(self.theta2[i]))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l1-self.l2-2,self.l1+self.l2+2), ylim=(-self.l1-self.l2-2,2)) 
        line, = ax.plot([], [], lw=2,c='b') 
        line2,=ax.plot([], [], lw=2,c='r')
        line3,=ax.plot([], [], lw=1,c='g')
        line4,=ax.plot([], [], lw=1,c='m')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.15)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.15)
            t1,t2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=sline(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            line.set_data(t1,t2)
            line2.set_data(s1,s2)
            line3.set_data(self.x1[:i],self.y1[:i])
            line4.set_data(self.x2[:i],self.y2[:i])
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            return line,line2,line3,line4,cline1,cline2,
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()   
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l1-self.l2-2,self.l1+self.l2+2), ylim=(-self.l1-self.l2-2,2)) 
        line, = ax.plot([], [], lw=2,c='b')
        line2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.15)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.15)
            t1,t2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=sline(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            line.set_data(t1,t2)
            line2.set_data(s1,s2)
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            return line,line2,cline1,cline2
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()   
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x1,self.y1,c='b',label='pendulum1')
        plt.plot(self.x2,self.y2,c='m',label='pendulum2')
        plt.legend()
        plt.show()
    def show_phase_space_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig, axs = plt.subplots(2,figsize=(12,12))
        axs[0].plot(self.theta1,self.dtheta1,c='r')
        axs[1].plot(self.theta2,self.dtheta2,c='b')
        axs[0].set_title('phase space plot pendulm1')
        axs[1].set_title('phase space plot pendulum2')
        for ax in axs.flat:
            ax.set(xlabel='position', ylabel='momentum')
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig, axs = plt.subplots(2,figsize=(12,12))
        axs[0].plot(self.t,self.theta1,c='r')
        axs[1].plot(self.t,self.theta2,c='b')
        axs[0].set_title('position time plot pendulm1')
        axs[1].set_title('position time plot pendulum2')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='angle')