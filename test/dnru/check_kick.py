import numpy as np
import matplotlib.pyplot as pl

f=open('lftr.out','r')

# parameters
print 'Value_1 = 1.2e-7'
print 'Value_2 = 1.7e-7'
print 'Value_3 = 0.3'
print 'Value_3 = 0.8'
print '=> 0.84e-7 < theta_x < 1.2e-7'
print '   0.34e-7 < theta_y < 1.7e-7'

coord0=[]
coord1=[]

pl.close('all')

with open('lftr.out','r') as f:
  while True:
    line=f.readline()
    if not line: break
    if 'before' in line:
      x,xp,y=line.split()[1:]
      yp,z,dp = f.readline().split()
      coord0.append((x,xp,y,yp,z,dp))
    if 'after' in line:
      x,xp,y=line.split()[1:]
      yp,z,dp = f.readline().split()
      coord1.append((x,xp,y,yp,z,dp))

coord0 = np.array(coord0,dtype=float)
coord1 = np.array(coord1,dtype=float)

pl.figure(figsize=(8,8))
for i in range(6):
  pl.subplot(3,2,i+1)
  # only take data which is not zero
  c = coord1[:,i]-coord0[:,i]
  pl.hist(c[c!=0]*1.e7,label='x(%s)'%i)
  pl.legend()
  pl.tight_layout()
  pl.xlabel('x(%s)*1.e7'%i)

pl.draw()
pl.show()
