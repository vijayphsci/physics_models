import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import animation
c1=0.5
c2=2
vel=3
def x1(t):
    return c1*t
def u1(t):
    return c1
def y1(t):
    return c2*t
def v1(t):
    return c2

def f(p1,q1,p2,q2):
    return vel*(p1-p2)/np.sqrt((p1-p2)**2+(q1-q2)**2)

def g(p1,q1,p2,q2):
    return vel*(q1-q2)/np.sqrt((p1-p2)**2+(q1-q2)**2)
dt=0.1
ti=0
tf=5
t=[]
p1x=[]
p1y=[]
p2x=[]
p2y=[]
p2_xini=1
p2_yini=0
t.append(ti)
p2x.append(p2_xini)
p2y.append(p2_yini)
p1x.append(x1(ti))
p1y.append(y1(ti))
n=int((tf-ti)/dt)

def dist(x1,y1,x2,y2):
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

'''
for i in range(n):
    t.append(ti+(i+1)*dt)
    p1x.append(x1(t[i+1]))
    p1y.append(y1(t[i+1]))
    t1=p2x[i]+dt*f(p1x[i]+dt/2*u1(t[i]),p1y[i]+dt/2*v1(t[i]),p2x[i]+dt/2*f(p1x[i],p1y[i],p2x[i],p2y[i]),p2y[i]+dt/2*g(p1x[i],p1y[i],p2x[i],p2y[i]))
    t2=p2y[i]+dt*g(p1x[i]+dt/2*u1(t[i]),p1y[i]+dt/2*v1(t[i]),p2x[i]+dt/2*f(p1x[i],p1y[i],p2x[i],p2y[i]),p2y[i]+dt/2*g(p1x[i],p1y[i],p2x[i],p2y[i]))
    p2x.append(t1)
    p2y.append(t2)
'''
epsilon=5e-2
temp=2*epsilon
i=0
while temp>epsilon:
    t.append(ti+(i+1)*dt)
    p1x.append(x1(t[i+1]))
    p1y.append(y1(t[i+1]))
    t1=p2x[i]+dt*f(p1x[i]+dt/2*u1(t[i]),p1y[i]+dt/2*v1(t[i]),p2x[i]+dt/2*f(p1x[i],p1y[i],p2x[i],p2y[i]),p2y[i]+dt/2*g(p1x[i],p1y[i],p2x[i],p2y[i]))
    t2=p2y[i]+dt*g(p1x[i]+dt/2*u1(t[i]),p1y[i]+dt/2*v1(t[i]),p2x[i]+dt/2*f(p1x[i],p1y[i],p2x[i],p2y[i]),p2y[i]+dt/2*g(p1x[i],p1y[i],p2x[i],p2y[i]))
    p2x.append(t1)
    p2y.append(t2)
    i=i+1
    temp=dist(p1x[i],p1y[i],p2x[i],p2y[i])


fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0,10))
#ax.set_facecolor((0, 1, 1))
line, = ax.plot([], [], lw=2,c='r',marker='o')
line2, = ax.plot([], [], lw=2,c='b',marker='o')
def init():
    line2.set_data([],[])
    line.set_data([], [])
    return line,line2

x=[]
y=[]
xn=[]
yn=[]
def animate2(i):
    #x.append(p1x[i])
    #y.append(p1y[i])
    line2.set_data(p1x[i], p1y[i])
    #xn.append(p2x[i])
    #yn.append(p2y[i])
    line.set_data(p2x[i],p2y[i])
    return line , line2
anim = animation.FuncAnimation(fig, func=animate2, init_func=init,
                               frames=i, interval=20, blit=True,repeat=False)

plt.show()