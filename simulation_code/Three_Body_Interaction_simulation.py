#This Simulation Created by Vijay Kag 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as mt
from physical_objects import circle,sline,springcoil
class Three_Body_Interaction:
    def __init__(self,k12=-0.1,k13=-0.2,k23=-0.1,mass1=1.5,mass2=7,mass3=9,position1=(0,0),position2=(5,0),position3=(2.5,4),velocity1=(0,0),velocity2=(0,0),velocity3=(0,0),time=200,dt=0.1,frame_interval=15,figsize=(10,8),repeat=True):
        self.m1=mass1
        self.m2=mass2
        self.m3=mass3
        self.k12=k12
        self.k13=k13
        self.k23=k23
        self.r1i=position1
        self.r2i=position2
        self.r3i=position3
        self.v1i=velocity1
        self.v2i=velocity2
        self.v3i=velocity3
        self.tf=time
        self.dt=dt
        self.n=int(self.tf/self.dt)
        self.t=[]
        self.x1=[]
        self.y1=[]
        self.x2=[]
        self.y2=[]
        self.x3=[]
        self.y3=[]
        self.vx1=[]
        self.vy1=[]
        self.vx2=[]
        self.vy2=[]
        self.vx3=[]
        self.vy3=[]
        self.frame_interval=frame_interval
        self.figsize=figsize
        self.count=0
        self.ani=None
        self.repeat=repeat
    def f1(self,x1,y1,x2,y2,x3,y3):
        t1=mt.pow((x2-x1)**2+(y2-y1)**2,1.5)
        t2=mt.pow((x3-x1)**2+(y3-y1)**2,1.5)
        f1x= self.k12*(x1-x2)/t1+self.k13*(x1-x3)/t2
        f1y= self.k12*(y1-y2)/t1+self.k13*(y1-y3)/t2
        return f1x/self.m1,f1y/self.m1
    def f2(self,x1,y1,x2,y2,x3,y3):
        t1=mt.pow((x2-x1)**2+(y2-y1)**2,1.5)
        t2=mt.pow((x3-x2)**2+(y3-y2)**2,1.5)
        f2x= self.k12*(x2-x1)/t1+self.k23*(x2-x3)/t2
        f2y= self.k12*(y2-y1)/t1+self.k23*(y2-y3)/t2
        return f2x/self.m2,f2y/self.m2
    def f3(self,x1,y1,x2,y2,x3,y3):
        t1=mt.pow((x3-x1)**2+(y3-y1)**2,1.5)
        t2=mt.pow((x3-x2)**2+(y3-y2)**2,1.5)
        f3x= self.k13*(x3-x1)/t1+self.k23*(x3-x2)/t2
        f3y= self.k13*(y3-y1)/t1+self.k23*(y3-y2)/t2
        return f3x/self.m3,f3y/self.m3
    def __solve_system(self):
        self.t.append(0)
        self.x1.append(self.r1i[0])
        self.y1.append(self.r1i[1])
        self.x2.append(self.r2i[0])
        self.y2.append(self.r2i[1])
        self.x3.append(self.r3i[0])
        self.y3.append(self.r3i[1])
        self.vx1.append(self.v1i[0])
        self.vy1.append(self.v1i[1])
        self.vx2.append(self.v2i[0])
        self.vy2.append(self.v2i[1])
        self.vx3.append(self.v3i[0])
        self.vy3.append(self.v3i[1])
        for i in range(self.n):
            ax1,ay1=self.f1(self.x1[i],self.y1[i],self.x2[i],self.y2[i],self.x3[i],self.y3[i])
            ax2,ay2=self.f2(self.x1[i],self.y1[i],self.x2[i],self.y2[i],self.x3[i],self.y3[i])
            ax3,ay3=self.f3(self.x1[i],self.y1[i],self.x2[i],self.y2[i],self.x3[i],self.y3[i])
            nax1,nay1=self.f1(self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i],self.x3[i]+self.dt/2*self.vx3[i],self.y3[i]+self.dt/2*self.vy3[i])
            nax2,nay2=self.f2(self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i],self.x3[i]+self.dt/2*self.vx3[i],self.y3[i]+self.dt/2*self.vy3[i])
            nax3,nay3=self.f3(self.x1[i]+self.dt/2*self.vx1[i],self.y1[i]+self.dt/2*self.vy1[i],self.x2[i]+self.dt/2*self.vx2[i],self.y2[i]+self.dt/2*self.vy2[i],self.x3[i]+self.dt/2*self.vx3[i],self.y3[i]+self.dt/2*self.vy3[i])
            self.x1.append(self.x1[i]+self.dt*self.vx1[i]+self.dt**2/2*ax1)
            self.x2.append(self.x2[i]+self.dt*self.vx2[i]+self.dt**2/2*ax2)
            self.x3.append(self.x3[i]+self.dt*self.vx3[i]+self.dt**2/2*ax3)
            self.y1.append(self.y1[i]+self.dt*self.vy1[i]+self.dt**2/2*ay1)
            self.y2.append(self.y2[i]+self.dt*self.vy2[i]+self.dt**2/2*ay2)
            self.y3.append(self.y3[i]+self.dt*self.vy3[i]+self.dt**2/2*ay3)
            self.vx1.append(self.vx1[i]+self.dt*nax1)
            self.vx2.append(self.vx2[i]+self.dt*nax2)
            self.vx3.append(self.vx3[i]+self.dt*nax3)
            self.vy1.append(self.vy1[i]+self.dt*nay1)
            self.vy2.append(self.vy2[i]+self.dt*nay2)
            self.vy3.append(self.vy3[i]+self.dt*nay3)
    def show_particle_and_trajectory(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        min_x1,max_x1=min(self.x1),max(self.x1)
        min_y1,max_y1=min(self.y1),max(self.y1)
        min_x2,max_x2=min(self.x2),max(self.x2)
        min_y2,max_y2=min(self.y2),max(self.y2)
        min_x3,max_x3=min(self.x3),max(self.x3)
        min_y3,max_y3=min(self.y3),max(self.y3)
        min_x=min(min_x1,min_x2,min_x3)
        min_y=min(min_y1,min_y2,min_y3)
        max_x=max(max_x1,max_x2,max_x3)
        max_y=max(max_y1,max_y2,max_y3)
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min_x-2,max_x+2), ylim=(min_y-2,max_y+2)) 
        tline1,=ax.plot([], [], lw=1,c='r') 
        tline2,=ax.plot([], [], lw=1,c='g')
        tline3,=ax.plot([], [], lw=1,c='b')
        p1,=ax.plot([], [], lw=2,c='r',marker='o',label='particle1') 
        p2,=ax.plot([], [], lw=2,c='g',marker='o',label='particle2') 
        p3,=ax.plot([], [], lw=2,c='b',marker='o',label='particle3') 
        def animate(i): 
            p1.set_data(self.x1[i],self.y1[i])
            p2.set_data(self.x2[i],self.y2[i])
            p3.set_data(self.x3[i],self.y3[i])
            tline1.set_data(self.x1[:i],self.y1[:i])
            tline2.set_data(self.x2[:i],self.y2[:i])
            tline3.set_data(self.x3[:i],self.y3[:i])
            return tline1,tline2,tline3,p1,p2,p3
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
        min_x3,max_x3=min(self.x3),max(self.x3)
        min_y3,max_y3=min(self.y3),max(self.y3)
        min_x=min(min_x1,min_x2,min_x3)
        min_y=min(min_y1,min_y2,min_y3)
        max_x=max(max_x1,max_x2,max_x3)
        max_y=max(max_y1,max_y2,max_y3)
        fig = plt.figure(figsize=self.figsize) 
        ax = plt.axes(xlim=(min_x-2,max_x+2), ylim=(min_y-2,max_y+2)) 
        p1,=ax.plot([], [], lw=2,c='r',marker='o',label='particle1') 
        p2,=ax.plot([], [], lw=2,c='g',marker='o',label='particle2')
        p3,=ax.plot([], [], lw=2,c='b',marker='o',label='particle3')
        def animate(i): 
            p1.set_data(self.x1[i],self.y1[i])
            p2.set_data(self.x2[i],self.y2[i])
            p3.set_data(self.x3[i],self.y3[i])
            return p1,p2,p3
        self.ani = animation.FuncAnimation(fig, animate, 
							frames=self.n+1, interval=self.frame_interval, blit=True,repeat=self.repeat) 
        plt.legend()
        plt.show()
    def show_trajectory_plot(self):
        if self.count==0:
            self.__solve_system()
        self.count=self.count+1
        plt.figure(figsize=self.figsize)
        plt.plot(self.x1,self.y1,c='r',label='particle1')
        plt.plot(self.x2,self.y2,c='g',label='particle2')
        plt.plot(self.x3,self.y3,c='b',label='particle3')
        plt.legend()
        plt.show()         
