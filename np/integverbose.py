import math
import matplotlib.pyplot as pl
from astropy.time import Time

with open('inic.txt','r') as fil:
    def getf():
        return [float(s) for s in fil.readline().split()]
    sGM = getf()[0]
    eGM, ex,ey,ez,evx,evy,evz = getf()
    mGM, mx,my,mz,mvx,mvy,mvz = getf()


def inner(ax,ay,az,bx,by,bz):
    num = ax*bx + ay*by + az*bz
    dena = ax*ax + ay*ay + az*az
    denb = bx*bx + by*by + bz*bz
    return num/(dena*denb)**(1/2)


def pred(dy,c):
    date = Time(dy+2460000,format='jd')
    if c[1] < -.9998 and c[1] == min(c):
        print('maybe lunar',date.iso)
    if c[2] > .9998 and c[1] == max(c):
        print('maybe solar',date.iso)
        
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
    evx, evy, evz = evx+ax*dt, evy+ay*dt, evz+az*dt
    ax, ay, az = accel(sGM,mx,my,mz)
    mvx, mvy, mvz = mvx+ax*dt, mvy+ay*dt, mvz+az*dt
    ax, ay, az = accel(eGM,mx-ex,my-ey,mz-ez)
    mvx, mvy, mvz = mvx+ax*dt, mvy+ay*dt, mvz+az*dt
    
t = 0
mjd = []
dot = []

dt = 3600
day = 86400

yr2026 = 1405 * day
yr2100 = 28433 * day

out = 2

def checkpoint():
    global out
    if out > 1:
        print(' Sonne_GM         Erde_GM          Mond_GM          m^3 / s^2')
        print('%16.9e %16.9e %16.9e' % (sGM,eGM,mGM))
        print()
    if out > 0:
        date = Time(2460000+(2-out)*dt/86400,format='jd')
        print(date.iso)
        print(' Erde_x           Erde_y           Erde_z           m')        
        print('%16.9e %16.9e %16.9e' % (ex,ey,ez))
        print(' Mond_x           Mond_y           Mond_z           m') 
        print('%16.9e %16.9e %16.9e' % (mx,my,mz))
        print(' Erde_vx          Erde_vy          Erde_vz          m / s')
        print('%16.9e %16.9e %16.9e' % (evx,evy,evz))
        print(' Mond_vx          Mond_vy          Mond_vz          m / s') 
        print('%16.9e %16.9e %16.9e' % (mvx,mvy,mvz))
        print()
    out -= 1

while t < yr2026:
    # integration step
    checkpoint()
    drift(dt)
    kick(dt)
    drift(dt)
    t += dt
    checkpoint()
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

