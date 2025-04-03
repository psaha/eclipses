import numpy as np

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
            vel[0] = float(l[1].split()[0])
            vel[1] = float(l[2].split()[0])
            vel[2] = float(l[3])
            break
    pos *= 1e3/c
    vel *= 1e3/c
    fil.close()
    return M, pos, vel
    

N = 3
mass = np.zeros(N)
pos = np.zeros(shape=(N,3))
vel = 0*pos
mass[0],pos[0],vel[0] = extract('../equator/sol.txt')
mass[1],pos[1],vel[1] = extract('../equator/earth.txt')
mass[2],pos[2],vel[2] = extract('../equator/luna.txt')
'''
mass[3],pos[3],vel[3] = extract('mercury.txt')
mass[4],pos[4],vel[4] = extract('venus.txt')
mass[5],pos[5],vel[5] = extract('mars.txt')
mass[6],pos[6],vel[6] = extract('jupiter.txt')
mass[7],pos[7],vel[7] = extract('saturn.txt')
mass[8],pos[8],vel[8] = extract('uranus.txt')
mass[9],pos[9],vel[9] = extract('neptune.txt')
'''

cmass = 0
cpos = np.zeros(3)
cvel = 0*cpos
for n in range(N):
    cmass += mass[n]
    cpos += mass[n]*pos[n]
    cvel += mass[n]*vel[n]
cpos /= cmass
cvel /= cmass

pos -= cpos
vel -= cvel

data = np.zeros(shape=(N,7))
data[:,0] = mass
data[:,1:4] = pos
data[:,4:] = vel

np.savetxt('inic.txt',data)

