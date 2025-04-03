import numpy as np
import matplotlib.pyplot as pl
from astropy.time import Time

data = np.loadtxt('inic.txt')

mass = data[:,0]
pos = data[:,1:4]
vel = data[:,4:]

N = len(mass)

def inner(x,y):
    return np.sum(x*y)/(np.sum(x*x)*np.sum(y*y))**(1/2)

def pred(dy,c):
    date = Time(dy+2460000,format='jd')
    if c[1] < -.9998 and c[1] == min(c):
        print('maybe lunar',date.iso)
    if c[2] > .9998 and c[1] == max(c):
        print('maybe solar',date.iso)
        
def drift():
    global pos
    pos += vel*dt/2

def kick():
    global pos,vel
    acc = 0*pos
    for i in range(N):
        for j in range(i):
            dr = pos[i] - pos[j]
            denom = np.sum(dr**2)**(3/2)
            acc[i] -= mass[j]*dr/denom
            acc[j] += mass[i]*dr/denom
    vel += acc*dt
    
    
t = 0
mjd = []
dot = []

dt = 3600
day = 86400

while t < day*1300:
#for s in range(24*900):
    drift()
    kick()
    drift()
    t += dt
    mjd.append(t/day)
    srel = pos[0] - pos[1]
    mrel = pos[2] - pos[1]
    d = inner(mrel,srel)
    dot.append(d)
    if len(mjd) > 3:
        pred(mjd[-2],dot[-3:])

ang = np.arccos(dot)*180/np.pi

pl.plot(mjd,ang)
pl.xlabel('JD - 2460000')
pl.ylabel('Sun-Moon angle in degrees')
pl.show()

