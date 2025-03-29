import numpy as np
import matplotlib.pyplot as pl
from astropy.time import Time

data = np.loadtxt('inic.txt')

mass = data[:,0]
pos = data[:,1:4]
vel = data[:,4:]

N = len(mass)

def dot(x,y):
    return np.sum(x*y)/(np.sum(x*x)*np.sum(y*y))**(1/2)

def pred(t,c):
    t = Time(t+2460000,format='jd')
    if c[1] < -.9998 and c[1] == min(c):
        print('maybe lunar',t.iso)
    if c[2] > .9998 and c[1] == max(c):
        print('maybe solar',t.iso)
    

t = []
a = []

dt = 3600
for s in range(24*900):
    pos += vel*dt/2
    acc = 0*pos
    for i in range(N):
        for j in range(i):
            dr = pos[i] - pos[j]
            denom = np.sum(dr**2)**(3/2)
            acc[i] -= mass[j]*dr/denom
            acc[j] += mass[i]*dr/denom
    vel += acc*dt
    pos += vel*dt/2
    srel = pos[0] - pos[1]
    mrel = pos[2] - pos[1]
    t.append(s/24)
    d = dot(mrel,srel)
    a.append(d)
    if len(a) > 3:
        pred(t[-2],a[-3:])

ang = np.arccos(a)*180/np.pi

pl.plot(t,ang)
pl.xlabel('JD - 2460000')
pl.ylabel('Sun-Moon angle in degrees')
pl.show()

