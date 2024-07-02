import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits import mplot3d

pl.style.use('dark_background')

def extract(fname):
    c = 299792458
    pos = np.zeros(3)
    vel = 0*pos
    fil = open(fname)
    print(fname)
    mread = False
    while True:
        lyne = fil.readline()[:-1]
        if 'km^3/s^2' in lyne and not mread:
            gms = lyne.split('km^3/s^2')[1].split()
            q = gms.index('=')
            M = float(gms[q+1]) * 1e9/c**3
            mread = True
        if lyne.startswith(' X ='):
            l = lyne.split('=')
            pos[0] = float(l[1].split()[0])
            pos[1] = float(l[2].split()[0])
            pos[2] = float(l[3])
        if lyne.startswith(' VX='):
            l = lyne.split('=')
            #print(l)
            vel[0] = float(l[1].split()[0])
            vel[1] = float(l[2].split()[0])
            vel[2] = float(l[3])
            #print(vx,vy,vz)
            break
    pos *= 1e3/c
    vel *= 1e3/c
    fil.close()
    print(M,pos,vel)
    return M, pos, vel
    

N = 3
mass = np.zeros(N)
pos = np.zeros(shape=(N,3))
vel = 0*pos
mass[0],pos[0],vel[0] = extract('sol.txt')
mass[1],pos[1],vel[1] = extract('earth.txt')
mass[2],pos[2],vel[2] = extract('luna.txt')
'''
mass[3],pos[3],vel[3] = extract('mercury.txt')
mass[4],pos[4],vel[4] = extract('venus.txt')
mass[5],pos[5],vel[5] = extract('mars.txt')
mass[6],pos[6],vel[6] = extract('jupiter.txt')
mass[7],pos[7],vel[7] = extract('saturn.txt')
mass[8],pos[8],vel[8] = extract('uranus.txt')
mass[9],pos[9],vel[9] = extract('neptune.txt')
'''


print()


cmass = 0
cpos = np.zeros(3)
cvel = 0*cpos
for n in range(N):
    cmass += mass[n]
    cpos += mass[n]*pos[n]
    cvel += mass[n]*vel[n]
cpos /= cmass
cvel /= cmass

'''
print(mass)
print(pos)
print(vel)
print()

print('cm pos and vel')
print(cpos)
print(cvel)
print()
'''

pos -= cpos
vel -= cvel

for n in range(N):
    if n==0:
        print('%9.6f e-6 s' % (1e6*mass[n]))
    else:
        print('%9.6f e-12 s' % (1e12*mass[n]))
for n in range(N):
    if n!=1:
        num = pos[n]-pos[1]
        print('%11.6f %11.6f %11.6f s' % (num[0],num[1],num[2]))
        num = (vel[n]-vel[1])*1e6
        print('%11.6f %11.6f %11.6f e-6' % (num[0],num[1],num[2]))


def virial():
    cmass = 0
    cpos = np.zeros(3)
    cvel = 0*cpos
    for n in range(N):
        cmass += mass[n]
        cpos += mass[n]*pos[n]
        cvel += mass[n]*vel[n]
    kin = pot = 0
    for i in range(N):
        kin += mass[i] * np.sum(vel[i]**2)/2
    for j in range(i):
        dr = pos[i]-pos[j]
        dr = np.sum(dr**2)**(1/2)
        pot -= mass[i]*mass[j]/dr
    print(pot+kin)#,2*kin,pot)
        
def dot(x,y):
    return np.sum(x*y)/(np.sum(x*x)*np.sum(y*y))**(1/2)

def latlon():
    pi = np.pi
    cos = np.cos
    sin = np.sin
    for lat in np.linspace(-pi/2,pi/2,7):
        lon = np.linspace(-pi,pi,100)
        x = cos(lat)*cos(lon)
        y = cos(lat)*sin(lon)
        z = sin(lat)
        ax.plot(x, y, z, color='0.6')
    for lon in np.linspace(-pi,pi,13):
        lat = np.linspace(-pi/2,pi/2,100)
        x = cos(lat)*cos(lon)
        y = cos(lat)*sin(lon)
        z = sin(lat)
        ax.plot(x, y, z, color='0.6')

t = []
a = []
amin = amax = 0
        
# virial()



ax = pl.axes(projection='3d')

# Hide grid lines
ax.grid(False)

# Hide axes ticks
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

pl.gca().set_aspect('equal')
pl.tight_layout()

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
    if d < -0.9998 or d > 0.9998:
        print('%4i %8.5f' % (round(t[-1]),d))
    a.append(d)
    srel /= np.linalg.norm(srel)
    mrel /= np.linalg.norm(mrel)
    ter = (pos[1] - pos[0])/600
    lun = ter + (pos[2] - pos[1])/8
    srel = ter
    mrel = lun
    if s==0:
        latlon()
        pl.plot(0,0,0,'o',color='orange',markersize=30)
        sp, = pl.plot(srel[0],srel[1],srel[2],'o',color='lightblue',markersize=10)
        mp, = pl.plot(mrel[0],mrel[1],mrel[2],'o',color='silver')
    sp.set_xdata(srel[0])
    sp.set_ydata(srel[1])
    sp.set_3d_properties(srel[2])
    mp.set_xdata(mrel[0])
    mp.set_ydata(mrel[1])
    mp.set_3d_properties(mrel[2])
    pl.pause(0.01)


    
pl.gca().set_aspect('equal')
pl.show()

pl.plot(t,np.arccos(a)*180/np.pi)
pl.xlabel('JD - 2460000')
pl.ylabel('Sun-Moon angle in degrees')
pl.text(44,-6,'solar eclipse')
pl.text(60,180,'lunar eclipse')
pl.savefig('eclipses')
pl.show()

