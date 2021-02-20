import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import models as ps
plt.rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg'
#%%
def magsq(tup):
    return tup[0]**2+tup[1]**2
def mag(tup):
    return np.sqrt(magsq(tup))
def diff(tup1,tup2):
    return (tup1[0]-tup2[0],tup1[1]-tup2[1])
#%%
mass1=1000
mass2=10
mass3=0.1
velocity1=(0,0)
velocity2=(0,1)
velocity3=(0,2)
position1=(0,0)
position2=(5,0)
time=50
dt=0.01
#%%
print('distance of 2,3 must excced',mag(diff(position1,position2))*mass3*magsq(velocity3)/(mass2*magsq(velocity2)))
print('position1',position1)
#%%
position3=(6,0)
#%%
print('k12',-mass2*magsq(velocity2)*mag(diff(position1,position2)))
print('k23',-mass3*magsq(velocity3)*mag(diff(position2,position3)))
#%%
k12=-50
k23=-0.25
k13=0
#%%
m=ps.Three_Body_Interaction(time=150,dt=0.01,k12=k12,k13=k13,k23=k23,mass1=mass1,mass2=mass2,mass3=mass3,position1=position1,position2=position2,position3=position3,velocity1=velocity1,velocity2=velocity2,velocity3=velocity3)
m.show_trajectory_plot()
#%%
m.show_particle_and_trajectory()
#%%
f1 = r"D:\pcode\git_upload\particle.mp4" 
writermp4 = animation.FFMpegWriter(fps=100) 
m.ani.save(f1, writer=writermp4)
#%%
m.show_particle()
#%%
f1 = r"D:\pcode\git_upload\trajectory_particle.mp4" 
writermp4 = animation.FFMpegWriter(fps=60) 
m.ani.save(f1, writer=writermp4)


