#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Triple_Pendulum:
    def __init__(self,angle1=60,angle2=70,angle3=80,length1=4,length2=2.5,length3=2,mass1=1,mass2=0.5,mass3=0.2,gravity=2,time=200,dt=0.1,frame_interval=15,figsize=(10,8)):         
         self.m1=mass1
         self.m2=mass2
         self.m3=mass3
         self.l1=length1
         self.l2=length2
         self.l3=length3
         self.g=gravity
         self.theta1i_deg=angle1
         self.theta2i_deg=angle2
         self.theta3i_deg=angle3
         self.dtheta1i=0
         self.dtheta2i=0
         self.dtheta3i=0
         self.theta1i_rad=self.theta1i_deg*np.pi/180
         self.theta2i_rad=self.theta2i_deg*np.pi/180
         self.theta3i_rad=self.theta3i_deg*np.pi/180
         self.tf=time
         self.dt=dt
         self.frame_interval=frame_interval
         self.n=int(self.tf/dt)
         self.theta1=[]
         self.theta2=[]
         self.theta3=[]
         self.dtheta1=[]
         self.dtheta2=[]
         self.dtheta3=[]
         self.x1=[]
         self.y1=[]
         self.x2=[]
         self.y2=[]
         self.x3=[]
         self.y3=[]
         self.t=[]
         self.figsize=figsize
         self.count=0
         self.ani=None
    def eqn(self,f,th1,th2,th3,dth1,dth2,dth3):
         a1=(self.m1+self.m2+self.m3)*self.l1**2
         a2=(self.m2+self.m3)*self.l1*self.l2*np.cos(th1-th2)
         a3=self.m3*self.l1*self.l3*np.cos(th1-th3)
         b1=(self.m2+self.m3)*self.l1*self.l2*np.cos(th1-th2)
         b2=(self.m2+self.m3)*self.l2**2
         b3=self.m3*self.l2*self.l3*np.cos(th2-th3)
         c1=self.m3*self.l1*self.l3*np.cos(th1-th3)
         c2=self.m3*self.l2*self.l3*np.cos(th2-th3)
         c3=self.m3*self.l3**2
         d1=-(self.m1+self.m2+self.m3)*self.l1*self.g*np.sin(th1)+(self.m2+self.m3)*self.l1*self.l2*dth2**2*np.sin(th2-th1)+self.m3*self.l1*self.l3*dth3**2*np.sin(th3-th1)
         d2=-(self.m2+self.m3)*self.l2*self.g*np.sin(th2)+(self.m2+self.m3)*self.l1*self.l2*dth1**2*np.sin(th1-th2)+self.m3*self.l2*self.l3*dth3**2*np.sin(th3-th2)
         d3=-self.m3*self.l3*self.g*np.sin(th3)+self.m3*self.l2*self.l3*dth2**2*np.sin(th2-th3)+self.m3*self.l1*self.l3*dth1**2*np.sin(th1-th3)
         D=np.array([[a1,b1,c1],[a2,b2,c2],[a3,b3,c3]])
         Det=np.linalg.det(D)
         if f==1:
             D1=np.array([[d1,b1,c1],[d2,b2,c2],[d3,b3,c3]])
             return np.linalg.det(D1)/Det
         elif f==2:
             D2=np.array([[a1,d1,c1],[a2,d2,c2],[a3,d3,c3]])
             return np.linalg.det(D2)/Det
         else:
             D3=np.array([[a1,b1,d1],[a2,b2,d2],[a3,b3,d3]])
             return np.linalg.det(D3)/Det
    def __solve_system(self):
        self.t.append(0)
        self.theta1.append(self.theta1i_rad)
        self.theta2.append(self.theta2i_rad)
        self.theta3.append(self.theta3i_rad)
        self.dtheta1.append(self.dtheta1i)
        self.dtheta2.append(self.dtheta2i)
        self.dtheta3.append(self.dtheta3i)
        self.x1.append(self.l1*np.sin(self.theta1[0]))
        self.y1.append(-self.l1*np.cos(self.theta1[0]))
        self.x2.append(self.x1[0]+self.l2*np.sin(self.theta2[0]))
        self.y2.append(self.y1[0]-self.l2*np.cos(self.theta2[0]))
        self.x3.append(self.x2[0]+self.l3*np.sin(self.theta3[0]))
        self.y3.append(self.y2[0]-self.l3*np.cos(self.theta3[0]))
        for i in range(self.n):
            self.t.append((i+1)*self.dt)
            t1=self.eqn(1,self.theta1[i],self.theta2[i],self.theta3[i],self.dtheta1[i],self.dtheta2[i],self.dtheta3[i])
            t2=self.eqn(2,self.theta1[i],self.theta2[i],self.theta3[i],self.dtheta1[i],self.dtheta2[i],self.dtheta3[i])
            t3=self.eqn(3,self.theta1[i],self.theta2[i],self.theta3[i],self.dtheta1[i],self.dtheta2[i],self.dtheta3[i])
            self.theta1.append(self.theta1[i]+self.dt*self.dtheta1[i]+self.dt**2/2*t1)
            self.theta2.append(self.theta2[i]+self.dt*self.dtheta2[i]+self.dt**2/2*t2)
            self.theta3.append(self.theta3[i]+self.dt*self.dtheta3[i]+self.dt**2/2*t3)
            self.dtheta1.append(self.dtheta1[i]+self.dt*self.eqn(1,self.theta1[i]+self.dt/2*self.dtheta1[i],self.theta2[i]+self.dt/2*self.dtheta2[i],self.theta3[i]+self.dt/2*self.dtheta3[i],self.dtheta1[i]+self.dt/2*t1,self.dtheta2[i]+self.dt/2*t2,self.dtheta3[i]+self.dt/2*t3))
            self.dtheta2.append(self.dtheta2[i]+self.dt*self.eqn(2,self.theta1[i]+self.dt/2*self.dtheta1[i],self.theta2[i]+self.dt/2*self.dtheta2[i],self.theta3[i]+self.dt/2*self.dtheta3[i],self.dtheta1[i]+self.dt/2*t1,self.dtheta2[i]+self.dt/2*t2,self.dtheta3[i]+self.dt/2*t3))
            self.dtheta3.append(self.dtheta3[i]+self.dt*self.eqn(3,self.theta1[i]+self.dt/2*self.dtheta1[i],self.theta2[i]+self.dt/2*self.dtheta2[i],self.theta3[i]+self.dt/2*self.dtheta3[i],self.dtheta1[i]+self.dt/2*t1,self.dtheta2[i]+self.dt/2*t2,self.dtheta3[i]+self.dt/2*t3))
            self.x1.append(self.l1*np.sin(self.theta1[i]))
            self.y1.append(-self.l1*np.cos(self.theta1[i]))
            self.x2.append(self.x1[i]+self.l2*np.sin(self.theta2[i]))
            self.y2.append(self.y1[i]-self.l2*np.cos(self.theta2[i]))
            self.x3.append(self.x2[i]+self.l3*np.sin(self.theta3[i]))
            self.y3.append(self.y2[i]-self.l3*np.cos(self.theta3[i]))
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize)
        ax = plt.axes(xlim=(-self.l1-self.l2-self.l3-2,self.l1+self.l2+self.l3+2), ylim=(-self.l1-self.l2-self.l3-4,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        pline2,=ax.plot([], [], lw=2,c='r')
        pline3,=ax.plot([], [], lw=2,c='g')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        cline3,=ax.plot([], [], lw=2,c='k')
        tline1,=ax.plot([], [], lw=1,c='b')
        tline2,=ax.plot([], [], lw=1,c='r')
        tline3,=ax.plot([], [], lw=1,c='g')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.2)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.2)
            cx3,cy3=circle(self.x3[i],self.y3[i],0.2)
            tn1,tn2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=sline(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            u1,u2=sline(self.x2[i],self.y2[i],self.x3[i],self.y3[i])
            pline1.set_data(tn1,tn2)
            pline2.set_data(s1,s2)
            pline3.set_data(u1,u2)
            tline1.set_data(self.x1[:i],self.y1[:i])
            tline2.set_data(self.x2[:i],self.y2[:i])
            tline3.set_data(self.x3[:i],self.y3[:i])
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            cline3.set_data(cx3,cy3)
            return pline1,pline2,pline3,tline1,tline2,tline3,cline1,cline2,cline3
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.n+1, interval=self.frame_interval, blit=True) 
        plt.show()
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize)
        ax = plt.axes(xlim=(-self.l1-self.l2-self.l3-2,self.l1+self.l2+self.l3+2), ylim=(-self.l1-self.l2-self.l3-4,2)) 
        pline1,=ax.plot([], [], lw=2,c='b') 
        pline2,=ax.plot([], [], lw=2,c='r')
        pline3,=ax.plot([], [], lw=2,c='g')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        cline3,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.2)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.2)
            cx3,cy3=circle(self.x3[i],self.y3[i],0.2)
            tn1,tn2=sline(0,0,self.x1[i],self.y1[i])
            s1,s2=sline(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            u1,u2=sline(self.x2[i],self.y2[i],self.x3[i],self.y3[i])
            pline1.set_data(tn1,tn2)
            pline2.set_data(s1,s2)
            pline3.set_data(u1,u2)
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            cline3.set_data(cx3,cy3)
            return pline1,pline2,pline3,cline1,cline2,cline3
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.n+1, interval=self.frame_interval, blit=True) 
        plt.show()
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x1,self.y1,c='b',label='pendulum1')
        plt.plot(self.x2,self.y2,c='r',label='pendulum2')
        plt.plot(self.x3,self.y3,c='g',label='pendulum3')
        plt.legend()
        plt.show() 
    def show_phase_space_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig, axs = plt.subplots(3,figsize=(12,12))
        axs[0].plot(self.theta1,self.dtheta1,c='b')
        axs[1].plot(self.theta2,self.dtheta2,c='r')
        axs[2].plot(self.theta3,self.dtheta3,c='g')
        for ax in axs.flat:
            ax.set(xlabel='position', ylabel='momentum')
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1  
        fig, axs = plt.subplots(3,figsize=(12,12))
        axs[0].plot(self.t,self.theta1,c='b')
        axs[1].plot(self.t,self.theta2,c='r')
        axs[2].plot(self.t,self.theta3,c='g')
        axs[0].set_title('position time plot pendulm1')
        axs[1].set_title('position time plot pendulum2')
        axs[2].set_title('position time plot pendulum3')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='angle')