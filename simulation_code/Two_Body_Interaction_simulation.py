#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Two_Body_Interaction:
    def __init__(self,interaction_const=20,mass1=1,mass2=1,position1=(0,0),position2=(5,0),velocity1=(0,1),velocity2=(0,-1),time=150,dt=0.1,frame_interval=15,figsize=(10,8)):
        self.m1=mass1
        self.m2=mass2
        self.k12=interaction_const
        self.r1i=position1
        self.r2i=position2
        self.v1i=velocity1
        self.v2i=velocity2
        self.tf=time
        self.dt=dt
        self.n=int(self.tf/self.dt)
        self.t=[]
        self.x1=[]
        self.y1=[]
        self.x2=[]
        self.y2=[]
        self.vx1=[]
        self.vy1=[]
        self.vx2=[]
        self.vy2=[]
        self.frame_interval=frame_interval
        self.figsize=figsize
        self.count=0
        self.ani=None
    def f1(self,x1,y1,x2,y2):
        t1=mt.pow((x2-x1)**2+(y2-y1)**2,1.5)
        f1x=self.k12*(x2-x1)/t1
        f1y=self.k12*(y2-y1)/t1
        return f1x/self.m1,f1y/self.m1
    def f2(self,x1,y1,x2,y2):
        t1=mt.pow((x2-x1)**2+(y2-y1)**2,1.5)
        f2x= -self.k12*(x2-x1)/t1
        f2y= -self.k12*(y2-y1)/t1
        return f2x/self.m2,f2y/self.m2    
    def __solve_system(self):
        self.t.append(0)
        self.x1.append(self.r1i[0])
        self.y1.append(self.r1i[1])
        self.x2.append(self.r2i[0])
        self.y2.append(self.r2i[1])
        self.vx1.append(self.v1i[0])
        self.vy1.append(self.v1i[1])
        self.vx2.append(self.v2i[0])
        self.vy2.append(self.v2i[1])
        for i in range(self.n):
            ax1,ay1=self.f1(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            ax2,ay2=self.f2(self.x1[i],self.y1[i],self.x2[i],self.y2[i])
            nax1,nay1=self.f1(self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i])
            nax2,nay2=self.f2(self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i])
            self.t.append((i+1)*self.dt)    
            self.x1.append(self.x1[i]+self.dt*self.vx1[i]+self.dt**2/2*ax1)
            self.x2.append(self.x2[i]+self.dt*self.vx2[i]+self.dt**2/2*ax2)
            self.y1.append(self.y1[i]+self.dt*self.vy1[i]+self.dt**2/2*ay1)
            self.y2.append(self.y2[i]+self.dt*self.vy2[i]+self.dt**2/2*ay2)
            self.vx1.append(self.vx1[i]+self.dt*nax1)
            self.vx2.append(self.vx2[i]+self.dt*nax2)
            self.vy1.append(self.vy1[i]+self.dt*nay1)
            self.vy2.append(self.vy2[i]+self.dt*nay2)
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        min_x1,max_x1=min(self.x1),max(self.x1)
        min_y1,max_y1=min(self.y1),max(self.y1)
        min_x2,max_x2=min(self.x2),max(self.x2)
        min_y2,max_y2=min(self.y2),max(self.y2)
        min_x=min(min_x1,min_x2)
        min_y=min(min_y1,min_y2)
        max_x=max(max_x1,max_x2)
        max_y=max(max_y1,max_y2)
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min_x-2,max_x+2), ylim=(min_y-2,max_y+2)) 
        tline1,=ax.plot([], [], lw=1,c='r') 
        tline2,=ax.plot([], [], lw=1,c='m')
        p1,=ax.plot([], [], lw=2,c='r',marker='o',label='particle1') 
        p2,=ax.plot([], [], lw=2,c='m',marker='o',label='particle2') 
        def animate(i): 
            tline1.set_data(self.x1[:i],self.y1[:i])
            tline2.set_data(self.x2[:i],self.y2[:i])
            p1.set_data(self.x1[i],self.y1[i])
            p2.set_data(self.x2[i],self.y2[i])
            return tline1,tline2,p1,p2
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.n+1, interval=self.frame_interval, blit=True) 
        plt.legend()
        plt.show()
    def show_particle(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        min_x1,max_x1=min(self.x1),max(self.x1)
        min_y1,max_y1=min(self.y1),max(self.y1)
        min_x2,max_x2=min(self.x2),max(self.x2)
        min_y2,max_y2=min(self.y2),max(self.y2)
        min_x=min(min_x1,min_x2)
        min_y=min(min_y1,min_y2)
        max_x=max(max_x1,max_x2)
        max_y=max(max_y1,max_y2)
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min_x-2,max_x+2), ylim=(min_y-2,max_y+2)) 
        p1,=ax.plot([], [], lw=3,c='r',marker='o',label='particle1') 
        p2,=ax.plot([], [], lw=3,c='g',marker='o',label='particle2') 
        def animate(i): 
            p1.set_data(self.x1[i],self.y1[i])
            p2.set_data(self.x2[i],self.y2[i])
            return p1,p2
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.n+1, interval=self.frame_interval, blit=True) 
        plt.legend()
        plt.show()
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x1,self.y1,c='r',label='particle1')
        plt.plot(self.x2,self.y2,c='g',label='particle2')
        plt.legend()
        plt.show()         
    def show_position_time_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1 
        fig, axs = plt.subplots(2, 2,figsize=(15,15))
        axs[0, 0].plot(self.t, self.x1)
        axs[0, 0].set_title('particle1 x coordinarte')
        axs[0, 1].plot(self.t, self.y1, 'orange')
        axs[0, 1].set_title('particle1 y coordinarte')
        axs[1, 0].plot(self.t,self.x2, 'green')
        axs[1, 0].set_title('particle2 x coordinarte')
        axs[1, 1].plot(self.t, self.y2, 'red')
        axs[1, 1].set_title('particle2 y coordinarte')
        for ax in axs.flat:
            ax.set(xlabel='time', ylabel='position')