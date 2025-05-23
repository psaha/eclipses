import math
import matplotlib.pyplot as pl
from astropy.time import Time

oem = False

sGM, eGM, mGM = 1.327124400e+20,  3.986004350e+14,  4.902800070e+12
ex,ey,ez = -1.344847290e+11,  6.189060650e+10, -2.459272850e+06
mx,my,mz = -1.341622790e+11,  6.208933150e+10, -4.759447620e+06
evx,evy,evz = -1.293204740e+04, -2.718492330e+04,  3.784017490e-01
mvx,mvy,mvz = -1.340964130e+04, -2.626952280e+04,  9.136169430e+01

def inner(ax,ay,az,bx,by,bz):
    num = ax*bx + ay*by + az*bz
    dena = ax*ax + ay*ay + az*az
    denb = bx*bx + by*by + bz*bz
    return num/(dena*denb)**(1/2)


def pred(dy,c):
    date = Time(dy+2460000,format='jd')
    if c[1] < -.9998 and c[1] == min(c):
        print('maybe lunar',date.iso,' Julian date %.2f' % (dy+2460000))
    if c[2] > .9998 and c[1] == max(c):
        print('maybe solar',date.iso,' Julian date %.2f' % (dy+2460000))
        
        
        
def drift(dt):
    global ex,ey,ez,mx,my,mz,evx,evy,evz,mvx,mvy,mvz
    ex += evx * dt/2
    ey += evy * dt/2
    ez += evz * dt/2
    mx += mvx * dt/2
    my += mvy * dt/2
    mz += mvz * dt/2
    
def accel(GM,x,y,z):
    f = -GM/(x*x + y*y + z*z)**(3/2)
    return x*f, y*f, z*f

def kick(dt):
    global ex,ey,ez,mx,my,mz,evx,evy,evz,mvx,mvy,mvz
    ax, ay, az = accel(sGM,ex,ey,ez)
    evx, evy, evz = evx+ax*dt, evy+ay*dt, evz+az*dt
    ax, ay, az = accel(mGM,ex-mx,ey-my,ez-mz)
    if not oem:
        evx, evy, evz = evx+ax*dt, evy+ay*dt, evz+az*dt
    ax, ay, az = accel(sGM,mx,my,mz)
    mvx, mvy, mvz = mvx+ax*dt, mvy+ay*dt, mvz+az*dt
    ax, ay, az = accel(eGM,mx-ex,my-ey,mz-ez)
    if not oem:
        mvx, mvy, mvz = mvx+ax*dt, mvy+ay*dt, mvz+az*dt
    
t = 0
mjd = []
dot = []

dt = 3600
day = 86400

yr2026 = 1405 * day
yr2100 = 28433 * day

while t < yr2026:
    # integration step
    drift(dt)
    kick(dt)
    drift(dt)
    t += dt
    # save and check for eclipses
    mjd.append(t/day)
    d = inner(ex-mx,ey-my,ez-mz,ex,ey,ez)
    dot.append(d)
    if len(mjd) > 3:
        pred(mjd[-2],dot[-3:])

    
ang = [math.acos(x)*180/math.pi for x in dot]

pl.plot(mjd,ang)
pl.xlabel('JD - 2460000')
pl.ylabel('Sun-Moon angle in degrees')

pl.show()

