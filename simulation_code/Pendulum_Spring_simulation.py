#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Pendulum_and_Spring:
    """
    particle 1 pendulum
    particle 2 spring
    coordinate particle 1 = (x1,y1)
    coordinate particle 2 = (x2,y2)
    velocity particle 1 = (vx1,vy1)
    velocity particle 2 = (vx2,vy2)
    mass particle 1 = m1
    mass particle 2 = m2
    length of pendulum = l
    stiffness of spring = k
    spring natural length = l0
    gravity = g
    angle of pendulum along negative y axis = theta
    angular velocity vtheta=d(theta)/dt
    Lagrangian L of system:
        L= 1/2*m1*(l*vtheta)**2+1/2*m2*(x2**2+y2**2)+m1*g*l*cos(theta)-m2*g*y2-1/2*k*(sqrt((x2-l*sin(theta))**2+(y2+l*cos(theta))**2)-l0)**2
    dissipative term F:
        F = damping*((l*vtheta)**2+vx2**2+vy2**2)    
    equation of motion:
        d --> partial derivative
        dL/d(vtheta)-dL/d(theta)+dF/d(vtheta) = 0
        dL/d(vx2)-dL/d(x2)+dF/d(vx2) = 0
        dL/d(vy2)-dL/d(y2)+dF/d(vy2) = 0
    """
    def __init__(self,mass1=1,mass2=1,pendulum_length=5,angle=40,spring_position=(4,-3),stiffness=1,damping=0.01,gravity=2,time=200,dt=0.1,natural_length=4,spring_velocity=(0,0),pendulum_velocity=0,frame_interval=15,figsize=(10,8),repeat=True):
        self.m1=mass1
        self.m2=mass2
        self.l=pendulum_length
        self.thetai_deg=angle
        self.thetai_rad=self.thetai_deg*np.pi/180
        self.spring_pos=spring_position
        self.spring_vel=spring_velocity
        self.vthetai=pendulum_velocity
        self.k=stiffness
        self.g=gravity
        self.figsize=figsize
        self.frame_interval=frame_interval
        self.tf=time
        self.dt=dt
        self.l0=natural_length
        self.damp=damping
        self.theta=[]
        self.x2=[]
        self.y2=[]
        self.vx2=[]
        self.vy2=[]
        self.vtheta=[]
        self.x1=[]
        self.y1=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.count=0
        self.ani=None
        self.repeat=repeat
    def f1(self,x2,y2,theta,vtheta):
        t=np.sqrt(x2**2+y2**2+self.l**2-2*x2*self.l*np.sin(theta)+2*y2*self.l*np.cos(theta))
        return -self.g/self.l*np.sin(theta)+self.k*(x2*np.cos(theta)+y2*np.sin(theta))*(1-self.l0/t)/(self.m1*self.l)-2*self.damp*vtheta/self.m1
    def f2(self,x2,y2,theta,vx2):
        t=np.sqrt(x2**2+y2**2+self.l**2-2*x2*self.l*np.sin(theta)+2*y2*self.l*np.cos(theta))
        return self.k/self.m2*(self.l*np.sin(theta)-x2)*(1-self.l0/t)-2*self.damp*vx2/self.m2
    def f3(self,x2,y2,theta,vy2):
        t=np.sqrt(x2**2+y2**2+self.l**2-2*x2*self.l*np.sin(theta)+2*y2*self.l*np.cos(theta))
        return -self.g+self.k/self.m2*(self.l0/t-1)*(y2+self.l*np.cos(theta))-2*self.damp*vy2/self.m2
    def __solve_system(self):
        self.t.append(0)
        self.theta.append(self.thetai_rad)
        self.vtheta.append(self.vthetai)
        self.x2.append(self.spring_pos[0])
        self.y2.append(self.spring_pos[1])
        self.vx2.append(self.spring_vel[0])
        self.vy2.append(self.spring_vel[1])
        self.x1.append(self.l*np.sin(self.theta[0]))
        self.y1.append(-self.l*np.cos(self.theta[0]))
        for i in range(self.nt):
            t1=self.f1(self.x2[i],self.y2[i],self.theta[i],self.vtheta[i])
            t2=self.f2(self.x2[i],self.y2[i],self.theta[i],self.vx2[i])
            t3=self.f3(self.x2[i],self.y2[i],self.theta[i],self.vy2[i])
            self.t.append((i+1)*self.dt)
            self.theta.append(self.theta[i]+self.dt*self.vtheta[i]+self.dt**2/2*t1)
            self.x2.append(self.x2[i]+self.dt*self.vx2[i]+self.dt**2/2*t2)
            self.y2.append(self.y2[i]+self.dt*self.vy2[i]+self.dt**2/2*t3)
            self.vtheta.append(self.vtheta[i]+self.dt*self.f1(self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i],self.theta[i]+self.dt/2*self.vtheta[i],self.vtheta[i]+self.dt/2*t1))
            self.vx2.append(self.vx2[i]+self.dt*self.f2(self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i],self.theta[i]+self.dt/2*self.vtheta[i],self.vx2[i]+self.dt/2*t2))
            self.vy2.append(self.vy2[i]+self.dt*self.f3(self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i],self.theta[i]+self.dt/2*self.vtheta[i],self.vy2[i]+self.dt/2*t3))
            self.x1.append(self.l*np.sin(self.theta[i]))
            self.y1.append(-self.l*np.cos(self.theta[i]))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l+min(self.x2)-2,self.l+max(self.x2)+2), ylim=(-self.l+min(self.y2)-2,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        pline2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        tline1,=ax.plot([], [], lw=1,c='g')
        tline2,=ax.plot([], [], lw=1,c='m')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.3)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.3)
            tn1,tn2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=springcoil(self.x1[i],self.y1[i],self.x2[i],self.y2[i],10,0.1,60)
            if (self.x2[i]-self.x1[i])**2+(self.y2[i]-self.y1[i])**2<self.l0**2:
                pline2.set_color('lime')
            else:
                pline2.set_color('r')
            pline1.set_data(tn1,tn2)
            pline2.set_data(s1,s2)
            tline1.set_data(self.x1[:i],self.y1[:i])
            tline2.set_data(self.x2[:i],self.y2[:i])
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            return pline1,pline2,tline1,tline2,cline1,cline2
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()  
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l+min(self.x2)-2,self.l+max(self.x2)+2), ylim=(-self.l+min(self.y2)-2,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        pline2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.3)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.3)
            tn1,tn2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=springcoil(self.x1[i],self.y1[i],self.x2[i],self.y2[i],10,0.1,60)
            if (self.x2[i]-self.x1[i])**2+(self.y2[i]-self.y1[i])**2<self.l0**2:
                pline2.set_color('lime')
            else:
                pline2.set_color('r')
            pline1.set_data(tn1,tn2)
            pline2.set_data(s1,s2)
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            return pline1,pline2,cline1,cline2
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()  
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x1,self.y1,c='b',label='pendulum')
        plt.plot(self.x2,self.y2,c='m',label='spring')
        plt.legend()
        plt.show()        
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1 
        fig, axs = plt.subplots(2, 2,figsize=(15,15))
        axs[0, 0].plot(self.t, self.x1)
        axs[0, 0].set_title('pendulum particle x coordinarte')
        axs[0, 1].plot(self.t, self.y1, 'tab:orange')
        axs[0, 1].set_title('pendulum particle y coordinarte')
        axs[1, 0].plot(self.t,self.x2, 'tab:green')
        axs[1, 0].set_title('spring particle x coordinarte')
        axs[1, 1].plot(self.t, self.y2, 'tab:red')
        axs[1, 1].set_title('spring particle y coordinarte')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='position')
