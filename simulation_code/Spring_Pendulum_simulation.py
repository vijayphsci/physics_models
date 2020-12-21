#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Spring_and_Pendulum:
    def __init__(self,mass1=1,mass2=1,pendulum_length=3,angle=40,spring_position=(2,-2),stiffness=1,damping=0.0,gravity=2,time=200,dt=0.1,natural_length=4,spring_velocity=(0,0),pendulum_velocity=0,frame_interval=15,figsize=(10,8),repeat=True):
        self.m1=mass1
        self.m2=mass2
        self.k=stiffness
        self.tf=time
        self.gravity=gravity
        self.damp=damping
        self.dt=dt
        self.l0=natural_length
        self.l=pendulum_length
        self.spring_pos=spring_position
        self.x1i=self.spring_pos[0]
        self.y1i=self.spring_pos[1]
        self.spring_vel=spring_velocity
        self.vx1i=self.spring_vel[0]
        self.vy1i=self.spring_vel[1]
        self.vthetai=pendulum_velocity
        self.thetai_deg=angle
        self.thetai_rad=self.thetai_deg*np.pi/180
        self.frame_interval=frame_interval
        self.x1=[]
        self.y1=[]
        self.x2=[]
        self.y2=[]
        self.vx1=[]
        self.vy1=[]
        self.theta=[]
        self.vtheta=[]
        self.t=[]
        self.nt=int(self.tf/self.dt)
        self.figsize=figsize
        self.count=0
        self.ani=None
        self.repeat=repeat
    def eqn(self,n,x1,y1,th,vx1,vy1,vth):
        a1=self.m1+self.m2
        c1=self.m2*self.l*np.cos(th)
        c2=self.m2*self.l*np.sin(th)
        d1=self.k*x1*(self.l0/np.sqrt(x1**2+y1**2)-1)-self.damp*(vx1+2*vth*self.l*np.cos(th))+c2*vth**2
        b2=a1
        d2=self.k*y1*(self.l0/np.sqrt(x1**2+y1**2)-1)-self.damp*(vy1+2*vth*self.l*np.sin(th))-c1*vth**2-a1*self.gravity
        a3=c1
        b3=c2
        c3=self.m2*self.l**2
        d3=-self.gravity*c2-2*self.damp*(self.l**2*vth+vx1*self.l*np.cos(th)+vy1*self.l*np.sin(th))
        temp=(d3-a3*d1/a1-b3*d2/b2)/(c3-c2*b3/b2-c1*a3/a1)                             
        if n==3:
            return temp
        elif n==2:
            return (d2-c2*temp)/b2
        else:
            return (d1-c1*temp)/a1
    def __solve_system(self):
        self.x1.append(self.x1i)
        self.y1.append(self.y1i)
        self.vx1.append(self.vx1i)
        self.vy1.append(self.vy1i)
        self.t.append(0)
        self.theta.append(self.thetai_rad)
        self.vtheta.append(self.vthetai)
        self.x2.append(self.x1[0]+self.l*np.sin(self.theta[0]))
        self.y2.append(self.y1[0]-self.l*np.cos(self.theta[0]))
        for i in range(self.nt):
            self.t.append((i+1)*self.dt)
            e1=self.eqn(1,self.x1[i],self.y1[i],self.theta[i],self.vx1[i],self.vy1[i],self.vtheta[i])
            e2=self.eqn(2,self.x1[i],self.y1[i],self.theta[i],self.vx1[i],self.vy1[i],self.vtheta[i])
            e3=self.eqn(3,self.x1[i],self.y1[i],self.theta[i],self.vx1[i],self.vy1[i],self.vtheta[i])
            self.x1.append(self.x1[i]+self.dt*self.vx1[i]+self.dt**2/2*e1)
            self.vx1.append(self.vx1[i]+self.dt*self.eqn(1,self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.theta[i]+self.dt/2*self.vtheta[i],self.vx1[i]+self.dt/2*e1,self.vy1[i]+self.dt/2*e2,self.vtheta[i]+self.dt/2*e3))                                                
            self.y1.append(self.y1[i]+self.dt*self.vy1[i]+self.dt**2/2*e2)
            self.vy1.append(self.vy1[i]+self.dt*self.eqn(2,self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.theta[i]+self.dt/2*self.vtheta[i],self.vx1[i]+self.dt/2*e1,self.vy1[i]+self.dt/2*e2,self.vtheta[i]+self.dt/2*e3))
            self.theta.append(self.theta[i]+self.dt*self.vtheta[i]+self.dt**2/2*e3)
            self.vtheta.append(self.vtheta[i]+self.dt*self.eqn(3,self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.theta[i]+self.dt/2*self.vtheta[i],self.vx1[i]+self.dt/2*e1,self.vy1[i]+self.dt/2*e2,self.vtheta[i]+self.dt/2*e3))
            self.x2.append(self.x1[i]+self.l*np.sin(self.theta[i]))
            self.y2.append(self.y1[i]-self.l*np.cos(self.theta[i]))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(-self.l+min(self.x1)-2,self.l+max(self.x1)+2), ylim=(-self.l+min(self.y1)-2,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        pline2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        tline1,=ax.plot([], [], lw=1,c='g')
        tline2,=ax.plot([], [], lw=1,c='m')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.25)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.25)
            tn1,tn2=springcoil(0,0,self.x1[i],self.y1[i],10,0.1,60)
            s1,s2=sline(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            if self.x1[i]**2+self.y1[i]**2<self.l0**2:
                pline1.set_color('lime')
            else:
                pline1.set_color('b')
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
        ax = plt.axes(xlim=(-self.l+min(self.x1)-2,self.l+max(self.x1)+2), ylim=(-self.l+min(self.y1)-2,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        pline2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.25)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.25)
            tn1,tn2=springcoil(0,0,self.x1[i],self.y1[i],10,0.1,60)
            s1,s2=sline(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            pline1.set_data(tn1,tn2)
            if self.x1[i]**2+self.y1[i]**2<self.l0**2:
                pline1.set_color('lime')
            else:
                pline1.set_color('b')
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
        plt.plot(self.x1,self.y1,c='b',label='spring')
        plt.plot(self.x2,self.y2,c='m',label='pendulum')
        plt.legend()
        plt.show()
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1 
        fig, axs = plt.subplots(2, 2,figsize=(15,15))
        axs[0, 0].plot(self.t, self.x1)
        axs[0, 0].set_title('spring particle x coordinarte')
        axs[0, 1].plot(self.t, self.y1, 'orange')
        axs[0, 1].set_title('spring particle y coordinarte')
        axs[1, 0].plot(self.t,self.x2, 'green')
        axs[1, 0].set_title('pendulum particle x coordinarte')
        axs[1, 1].plot(self.t, self.y2, 'red')
        axs[1, 1].set_title('pendulum particle y coordinarte')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='position')
