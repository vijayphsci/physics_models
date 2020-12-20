import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Simple_Pendulum:
    def __init__(self,angle=40,length=5,gravity=2,damping=0.0,time=200,dt=0.1,velocity=0,frame_interval=10,figsize=(10,8),repeat=True):
        self.thetai_deg=angle
        self.thetai_rad=self.thetai_deg*np.pi/180
        self.l=length
        self.gravity=gravity
        self.damp=damping
        self.tf=time
        self.dt=dt
        self.vthetai=velocity
        self.frame_interval=frame_interval
        self.theta=[]
        self.vtheta=[]
        self.x=[]
        self.y=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.figsize=figsize
        self.count=0
        self.__ani=None
        self.repeat=repeat
    def f(self,th,vth):
        return -self.gravity*np.sin(th)/self.l-self.damp*vth
    def __solve_system(self):
        self.t.append(0)
        self.theta.append(self.thetai_rad)
        self.vtheta.append(self.vthetai)
        self.x.append(self.l*np.sin(self.theta[0]))
        self.y.append(-self.l*np.cos(self.theta[0]))
        for i in range(self.nt):
            t1=self.f(self.theta[i],self.vtheta[i])
            self.t.append((i+1)*self.dt)
            self.theta.append(self.theta[i]+self.dt*self.vtheta[i]+self.dt**2/2*t1)
            self.vtheta.append(self.vtheta[i]+self.dt*self.f(self.theta[i]+self.dt/2*self.vtheta[i],self.vtheta[i]+self.dt/2*t1))
            self.x.append(self.l*np.sin(self.theta[i]))
            self.y.append(-self.l*np.cos(self.theta[i]))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l-2,self.l+2), ylim=(-self.l-2,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        cline1,=ax.plot([], [], lw=2,c='r')
        tline1,=ax.plot([], [], lw=1,c='m')
        tj1x=[]
        tj1y=[]
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x[i],self.y[i],0.25)
            tj1x.append(self.x[i])
            tj1y.append(self.y[i])
            tn1,tn2=sline(0,0,self.x[i],self.y[i])
            pline1.set_data(tn1,tn2)
            tline1.set_data(tj1x,tj1y)
            cline1.set_data(cx1,cy1)
            return pline1,tline1,cline1
        self.__ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()           
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l-2,self.l+2), ylim=(-self.l-2,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        cline1,=ax.plot([], [], lw=2,c='r')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x[i],self.y[i],0.25)
            tn1,tn2=sline(0,0,self.x[i],self.y[i])
            pline1.set_data(tn1,tn2)
            cline1.set_data(cx1,cy1)
            return pline1,cline1
        self.__ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()   
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x,self.y,c='m',label='simple pendulum')
        plt.legend()
        plt.show()
    def show_phase_space_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        plt.title('Phase Space Plot')
        plt.xlabel('Position')
        plt.ylabel('Momentum')
        plt.plot(self.theta,self.vtheta,c='r')
        plt.show()
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig = plt.figure(figsize=self.figsize) 
        plt.title('Position time plot')
        plt.xlabel('time')
        plt.ylabel('angle')
        plt.plot(self.t,self.theta,c='r')
        plt.show()