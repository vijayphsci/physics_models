#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Double_Spring:
    """
    coordinate particle 1 = (x1,y1)
    coordinate particle 2 = (x2,y2)
    velocity particle 1 = (vx1,vy1)
    velocity particle 2 = (vx2,vy2)
    stifness of spring 1 = k1
    stifness of speing 2 = k2
    natural length spring 1 = l1
    natural length spring 2 = l2
    mass particle 1 = m1
    mass particle 2 = m2
    gravity = g
    Lagrangian L of system:
        L = 1/2*m1*(vx1**2+vy1**2)+1/2*m2*(vx2**2+vy2**2)-1/2*k1*(sqrt(x1**2+y1**2)-l1)**2-k2*(sqrt((x2-x1)**2+(y2-y1)**2)-l2)**2-m1*g*y1-m2*g*y2
    dissipative term F:
        F = damping/2*(vx1**2+vy1**2+vx2**2+vy2**2)
    equation of motion:
        d --> partial derivative
        dL/d(vx1)-dL/d(x1)+dF/d(vx1) = 0
        dL/d(vy1)-dL/d(y1)+dF/d(vy1) = 0
        dL/d(vx2)-dL/d(x2)+dF/d(vx2) = 0
        dL/d(vy2)-dL/d(y2)+dF/d(vy2) = 0
    
    """
    def __init__(self,mass1=1,mass2=2,stiffness1=1,stiffness2=2,gravity=3,damping=0.0,position1=(2,-2),position2=(4,-5),velocity1=(0,0),velocity2=(0,0),natural_length1=2,natural_length2=4,time=200,dt=0.1,frame_interval=15,figsize=(10,8),repeat=True):
        self.m1=mass1
        self.m2=mass2
        self.k1=stiffness1
        self.k2=stiffness2
        self.g=gravity
        self.damp=damping
        self.inipos1=position1
        self.inipos2=position2
        self.inivel1=velocity1
        self.inivel2=velocity2
        self.l1=natural_length1
        self.l2=natural_length2
        self.tf=time
        self.dt=dt
        self.frame_interval=frame_interval
        self.figsize=figsize
        self.x1=[]
        self.y1=[]
        self.x2=[]
        self.y2=[]
        self.t=[]
        self.vx1=[]
        self.vy1=[]
        self.vx2=[]
        self.vy2=[]
        self.nt=int(self.tf/self.dt)
        self.count=0
        self.ani=None
        self.repeat=repeat
    def f1(self,x1,x2,y1,y2,vx1):
        return 1/self.m1*(self.k1*x1*(self.l1/np.sqrt(x1**2+y1**2)-1)+2*self.k2*(x2-x1)*(1-self.l2/np.sqrt((x2-x1)**2+(y2-y1)**2))-self.damp*vx1)
    def f2(self,x1,x2,y1,y2,vy1):
        return 1/self.m1*(-self.m1*self.g+self.k1*y1*(self.l1/np.sqrt(x1**2+y1**2)-1)+2*self.k2*(y2-y1)*(1-self.l2/np.sqrt((x2-x1)**2+(y2-y1)**2))-self.damp*vy1)
    def f3(self,x1,x2,y1,y2,vx2):
        return 1/self.m2*(2*self.k2*(x1-x2)*(1-self.l2/np.sqrt((x2-x1)**2+(y2-y1)**2))-self.damp*vx2)
    def f4(self,x1,x2,y1,y2,vy2):
        return 1/self.m2*(-self.m2*self.g+2*self.k2*(y1-y2)*(1-self.l2/np.sqrt((x2-x1)**2+(y2-y1)**2))-self.damp*vy2)
    def __solve_system(self):
        self.t.append(0)
        self.x1.append(self.inipos1[0])
        self.y1.append(self.inipos1[1])
        self.x2.append(self.inipos2[0])
        self.y2.append(self.inipos2[1])
        self.vx1.append(self.inivel1[0])
        self.vy1.append(self.inivel1[1])
        self.vx2.append(self.inivel2[0])
        self.vy2.append(self.inivel2[1])
        for i in range(self.nt):
            self.t.append((i+1)*self.dt)
            t1=self.f1(self.x1[i],self.x2[i],self.y1[i],self.y2[i],self.vx1[i])
            t2=self.f2(self.x1[i],self.x2[i],self.y1[i],self.y2[i],self.vy1[i])
            t3=self.f3(self.x1[i],self.x2[i],self.y1[i],self.y2[i],self.vx2[i])
            t4=self.f4(self.x1[i],self.x2[i],self.y1[i],self.y2[i],self.vy2[i])
            self.x1.append(self.x1[i]+self.dt*self.vx1[i]+self.dt**2/2*t1)
            self.vx1.append(self.vx1[i]+self.dt*self.f1(self.x1[i]+self.dt/2*self.vx1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y1[i]+self.dt/2*self.vy1[i],self.y2[i]+self.dt/2*self.vy2[i],self.vx1[i]+self.dt/2*t1))
            self.y1.append(self.y1[i]+self.dt*self.vy1[i]+self.dt**2/2*t2)
            self.vy1.append(self.vy1[i]+self.dt*self.f2(self.x1[i]+self.dt/2*self.vx1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y1[i]+self.dt/2*self.vy1[i],self.y2[i]+self.dt/2*self.vy2[i],self.vy1[i]+self.dt/2*t2))
            self.x2.append(self.x2[i]+self.dt*self.vx2[i]+self.dt**2/2*t3)
            self.vx2.append(self.vx2[i]+self.dt*self.f3(self.x1[i]+self.dt/2*self.vx1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y1[i]+self.dt/2*self.vy1[i],self.y2[i]+self.dt/2*self.vy2[i],self.vx2[i]+self.dt/2*t3))
            self.y2.append(self.y2[i]+self.dt*self.vy2[i]+self.dt**2/2*t4)
            self.vy2.append(self.vy2[i]+self.dt*self.f4(self.x1[i]+self.dt/2*self.vx1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y1[i]+self.dt/2*self.vy1[i],self.y2[i]+self.dt/2*self.vy2[i],self.vy2[i]+self.dt/2*t4))
    
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min(min(self.x1),min(self.x2))-2,max(max(self.x1),max(self.x2))+2), ylim=(min(min(self.y1),min(self.y2))-2,4)) 
        line, = ax.plot([], [], lw=2,c='b') 
        line2,=ax.plot([], [], lw=2,c='r')
        line3,=ax.plot([], [], lw=1,c='g')
        line4,=ax.plot([], [], lw=1,c='m')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.3)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.3)
            t1,t2=springcoil(0,0,self.x1[i],self.y1[i],10,0.1,60)
            s1,s2=springcoil(self.x1[i],self.y1[i],self.x2[i],self.y2[i],10,0.1,60)
            line.set_data(t1,t2)
            line2.set_data(s1,s2)
            line3.set_data(self.x1[:i],self.y1[:i])
            line4.set_data(self.x2[:i],self.y2[:i])
            cline1.set_data(cx1,cy1)
            cline2.set_data(cx2,cy2)
            return line,line2,line3,line4,cline1,cline2
        self.ani = animation.FuncAnimation(fig, animate, 
						frames=self.nt+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.show()   
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min(min(self.x1),min(self.x2))-2,max(max(self.x1),max(self.x2))+2), ylim=(min(min(self.y1),min(self.y2))-2,4))
        line, = ax.plot([], [], lw=2,c='b')
        line2,=ax.plot([], [], lw=2,c='r')
        cline1,=ax.plot([], [], lw=2,c='k')
        cline2,=ax.plot([], [], lw=2,c='k')
        plt.hlines(0,-2,2)
        def animate(i): 
            cx1,cy1=circle(self.x1[i],self.y1[i],0.3)
            cx2,cy2=circle(self.x2[i],self.y2[i],0.3)
            t1,t2=springcoil(0,0,self.x1[i],self.y1[i],10,0.1,60)
            s1,s2=springcoil(self.x1[i],self.y1[i],self.x2[i],self.y2[i],10,0.1,60)
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
        plt.plot(self.x1,self.y1,c='b',label='spring1')
        plt.plot(self.x2,self.y2,c='m',label='spring2')
        plt.legend()
        plt.show()
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1 
        fig, axs = plt.subplots(2, 2,figsize=(15,15))
        axs[0, 0].plot(self.t, self.x1)
        axs[0, 0].set_title('particle1 x coordinarte')
        axs[0, 1].plot(self.t, self.y1, 'tab:orange')
        axs[0, 1].set_title('particle1 y coordinarte')
        axs[1, 0].plot(self.t,self.x2, 'tab:green')
        axs[1, 0].set_title('particle2 x coordinarte')
        axs[1, 1].plot(self.t, self.y2, 'tab:red')
        axs[1, 1].set_title('particle2 y coordinarte')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='position')
