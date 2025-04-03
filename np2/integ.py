import numpy as np
import matplotlib.pyplot as pl
from astropy.time import Time

data = np.loadtxt('inic.txt')

sGM, eGM, mGM = data[:,0]
epos = data[1,1:4]
evel = data[1,4:]
mpos = data[2,1:4]
mvel = data[2,4:]

def inner(x,y):
    return np.sum(x*y)/(np.sum(x*x)*np.sum(y*y))**(1/2)

def pred(dy,c):
    date = Time(dy+2460000,format='jd')
    if c[1] < -.9998 and c[1] == min(c):
        print('maybe lunar',date.iso)
    if c[2] > .9998 and c[1] == max(c):
        print('maybe solar',date.iso)
        
def drift(dt):
    global epos,mpos
    epos += evel * dt/2
    mpos += mvel * dt/2
    
def accel(GM,dis):
    return -GM * dis/np.sum(dis**2)**(3/2)

def kick(dt):
    global epos,mpos,evel,mvel
    evel += accel(sGM,epos) * dt
    evel += accel(mGM,epos-mpos) * dt
    mvel += accel(sGM,mpos) * dt
    mvel += accel(eGM,mpos-epos) * dt
    
t = 0
mjd = []
dot = []

dt = 3600
day = 86400

while t < day*1300:
    drift(dt)
    kick(dt)
    drift(dt)
    t += dt
    mjd.append(t/day)
    d = inner(epos-mpos,epos)
    dot.append(d)
    if len(mjd) > 3:
        pred(mjd[-2],dot[-3:])

ang = np.arccos(dot)*180/np.pi

pl.plot(mjd,ang)
pl.xlabel('JD - 2460000')
pl.ylabel('Sun-Moon angle in degrees')
pl.show()

