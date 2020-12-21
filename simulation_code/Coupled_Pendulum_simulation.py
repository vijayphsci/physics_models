#This Simulation is Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Coupled_Pendulum:
    """
    mass particle 1 = mass1
    mass particle 2 = mass2
    initial angle particle 1 = angle1
    initial angle particle 2 = angle2
    gravity , g  = gravity
    stiffness constant of spring , k = stiffness  
    damping constant , damp = damping
    length of pendulum , l = pendulum_length
    ratio of length at spring attached to length of pendulum  = spring_position
    sepration of pendulum , l0 = sepration
    position of spring , l1 = spring_position*pendulum_length
    final time ,tf = time
    angle particle 1 = theta1
    angle particle 2 = theta2
    velocity particle 1 , vtheta1=  d(theta1)/dt
    velocity particle 1 , vtheta2=  d(theta2)/dt
    repeat =True  (animation is repeated again)
    Lagrangian of system:
        s1 = position at which spring attached to pendulum1 , (l1*sin(theta1),-l1*cos(theta1))
        s2 = position at which spring attached to pendulum1 , (l0+l1*sin(theta1),-l1*cos(theta2))
        s=distance(s1,s2)
        s=sqrt(2*l1**2*(1-cos(theta1-theta2))+l0**2+2*l1*l0*(sin(theta1)-sin(theta2)))
        L= 1/2*m1*(l*vtheta1)**2+1/2*m2*(l*vtheta2)**2+m1*g*l*cos(theta1)+m2*g*l*cos(theta2)-1/2*k*(s-l0)**2
    equation of motion:
        D - > exact diffrential 
        d - > partial derivative
        dissipative function ,F = damping/2*(vtheta1**2+vtheta2**2)
        D/dt(dL/d(vtheta1))-dL/d(theta1)+dF/(dvtheta1) = 0
        D/dt(dL/d(vtheta2))-dL/d(theta2)+dF/(dvtheta2) = 0
        
    """
    def __init__(self,mass1=1,mass2=1,angle1=0,angle2=40,gravity=3,stiffness=0.5,damping=0,pendulum_length=5,spring_position=0.3,sepration=4,time=200,dt=0.1,frame_interval=15,figsize=(10,8),repeat=True):
        self.m1=mass1
        self.m2=mass2
        self.theta1_ini_deg=angle1
        self.theta2_ini_deg=angle2
        self.theta1_ini_rad=self.theta1_ini_deg*np.pi/180
        self.theta2_ini_rad=self.theta2_ini_deg*np.pi/180
        self.g=gravity
        self.k=stiffness
        self.tf=time
        self.dt=dt
        self.nt=int(self.tf/self.dt)
        self.l=pendulum_length
        self.l1=self.l*spring_position
        self.l0=sepration
        self.figsize=figsize
        self.frame_interval=frame_interval
        self.damp=damping
        self.t=[]
        self.theta1=[]
        self.theta2=[]
        self.vtheta1=[]
        self.vtheta2=[]
        self.x1=[]
        self.y1=[]
        self.x2=[]
        self.y2=[]
        self.s1x=[]
        self.s1y=[]
        self.s2x=[]
        self.s2y=[]
        self.count=0
        self.ani=None
        self.repeat=repeat
    def f1(self,theta1,theta2,vtheta1):
        s=np.sqrt(2*self.l1**2*(1-np.cos(theta1-theta2))+self.l0**2+2*self.l1*self.l0*(np.sin(theta1)-np.sin(theta2)))
        return -self.g/self.l*np.sin(theta1)+self.k/(self.m1*self.l**2)*(self.l0/s-1)*(self.l0*self.l1*np.cos(theta1)+self.l1**2*np.sin(theta1-theta2))-self.damp/self.m1*vtheta1        
    def f2(self,theta1,theta2,vtheta2):
        s=np.sqrt(2*self.l1**2*(1-np.cos(theta1-theta2))+self.l0**2+2*self.l1*self.l0*(np.sin(theta1)-np.sin(theta2)))
        return -self.g/self.l*np.sin(theta2)+self.k/(self.m2*self.l**2)*(self.l0/s-1)*(-self.l0*self.l1*np.cos(theta2)+self.l1**2*np.sin(theta2-theta1))-self.damp/self.m2*vtheta2
    def __solve_system(self):
        self.t.append(0)
        self.theta1.append(self.theta1_ini_rad)
        self.theta2.append(self.theta2_ini_rad)
        self.vtheta1.append(0)
        self.vtheta2.append(0)
        self.x1.append(self.l*np.sin(self.theta1[0]))
        self.y1.append(-self.l*np.cos(self.theta1[0]))
        self.x2.append(self.l0+self.l*np.sin(self.theta2[0]))
        self.y2.append(-self.l*np.cos(self.theta2[0]))
        self.s1x.append(self.l1*np.sin(self.theta1[0]))
        self.s1y.append(-self.l1*np.cos(self.theta1[0]))
        self.s2x.append(self.l0+self.l1*np.sin(self.theta2[0]))
        self.s2y.append(-self.l1*np.cos(self.theta2[0]))
        for i in range(self.nt):
            self.t.append((i+1)*self.dt)
            t1=self.f1(self.theta1[i],self.theta2[i],self.vtheta1[i])
            t2=self.f2(self.theta1[i],self.theta2[i],self.vtheta2[i])
            self.theta1.append(self.theta1[i]+self.dt*self.vtheta1[i]+self.dt**2/2*t1)
            self.vtheta1.append(self.vtheta1[i]+self.dt*self.f1(self.theta1[i]+self.dt/2*self.vtheta1[i],self.theta2[i]+self.dt/2*self.vtheta2[i],self.vtheta1[i]+self.dt/2*t1))
            self.theta2.append(self.theta2[i]+self.dt*self.vtheta2[i]+self.dt**2/2*t2)
            self.vtheta2.append(self.vtheta2[i]+self.dt*self.f2(self.theta1[i]+self.dt/2*self.vtheta1[i],self.theta2[i]+self.dt/2*self.vtheta2[i],self.vtheta2[i]+self.dt/2*t2))
            self.x1.append(self.l*np.sin(self.theta1[i]))
            self.y1.append(-self.l*np.cos(self.theta1[i]))
            self.x2.append(self.l0+self.l*np.sin(self.theta2[i]))
            self.y2.append(-self.l*np.cos(self.theta2[i]))
            self.s1x.append(self.l1*np.sin(self.theta1[i]))
            self.s1y.append(-self.l1*np.cos(self.theta1[i]))
            self.s2x.append(self.l0+self.l1*np.sin(self.theta2[i]))
            self.s2y.append(-self.l1*np.cos(self.theta2[i]))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l-2,self.l0+self.l+2), ylim=(-self.l-2,2)) 
        line, = ax.plot([], [], lw=2,c='b') 
        line2,=ax.plot([], [], lw=2,c='r')
        line3,=ax.plot([], [], lw=1,c='g')
        line4,=ax.plot([], [], lw=1,c='m')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        spline,=ax.plot([], [], lw=2,c='c')
        plt.hlines(0,-2,self.l0+2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.25)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.25)
            t1,t2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=sline(self.l0,0,self.x2[i],self.y2[i])
            sp1,sp2=springcoil(self.s1x[i],self.s1y[i],self.s2x[i],self.s2y[i],10,0.1,60)
            if np.sqrt((self.s1x[i]-self.s2x[i])**2+(self.s1y[i]-self.s2y[i])**2)>self.l0:
                spline.set_color('lime')
            else:
                spline.set_color('c')            
            line.set_data(t1,t2)
            line2.set_data(s1,s2)
            line3.set_data(self.x1[:i],self.y1[:i])
            line4.set_data(self.x2[:i],self.y2[:i])
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            spline.set_data(sp1,sp2)
            return line,line2,line3,line4,cline1,cline2,spline
        self.ani = animation.FuncAnimation(fig, animate,
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()   
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l-2,self.l0+self.l+2), ylim=(-self.l-2,2)) 
        line, = ax.plot([], [], lw=2,c='b') 
        line2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        spline,=ax.plot([], [], lw=2,c='c')
        plt.hlines(0,-2,self.l0+2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.2)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.2)
            t1,t2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=sline(self.l0,0,self.x2[i],self.y2[i])
            sp1,sp2=springcoil(self.s1x[i],self.s1y[i],self.s2x[i],self.s2y[i],10,0.1,60)
            if np.sqrt((self.s1x[i]-self.s2x[i])**2+(self.s1y[i]-self.s2y[i])**2)>self.l0:
                spline.set_color('lime')
            else:
                spline.set_color('c')
            line.set_data(t1,t2)
            line2.set_data(s1,s2)
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            spline.set_data(sp1,sp2)
            return line,line2,cline1,cline2,spline
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
        axs[0].plot(self.theta1,self.vtheta1,c='r')
        axs[1].plot(self.theta2,self.vtheta2,c='b')
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
