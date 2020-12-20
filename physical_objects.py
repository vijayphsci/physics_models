import numpy as np
def circle(xi,yi,radius):
    ang=np.linspace(0,2*np.pi,20)
    x=xi+radius*np.cos(ang)
    y=yi+radius*np.sin(ang)
    return x,y
def sline(xi,yi,xf,yf):
    x=[xi,xf]
    y=[yi,yf]
    return x,y
def springcoil(x1,y1,x2,y2,cycle,amplitude,linspaces):
    l=np.sqrt((x2-x1)**2+(y2-y1)**2)
    xp=np.linspace(0,l,linspaces)
    yp=amplitude*np.sin((2*cycle)*np.pi*xp/l)
    newx=(xp*(x2-x1)-yp*(y2-y1))/l+x1
    newy=(xp*(y2-y1)+yp*(x2-x1))/l+y1
    return newx,newy