# import modules
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
#filepath for output MATLAB file
filepath ='/home/gln09/Desktop/GN_Ubiquitination/distributive/distributivesampling'

# Define the terms of the equations
# Y = [protein]
# t = time

# knumber = rate of that reaction
k = []
for a in range(0,15):
    k.append(0)
for a in range(0,15):
    k[a]=1
#print "k = ",k, len(k)
ks = []

ks=np.arange(0,26,1).tolist()

#print ks    
L = 1
# L = [E3 ligase]

parameters = (k, L)
iX0 = 5
iLX0 = 0
Xu = 0
LXu = 0
Xuu = 0
LXuu = 0
Xuuu = 0
LXuuu = 0
Xuuuu = 0
LXuuuu = 0
Xuuuuu = 0
y = []
for a in range(0,11):
    y.append(0)

y[0] = iX0
y[1] = iLX0
y[2] = Xu
y[3] = LXu
y[4] = Xuu
y[5] = LXuu
y[6] = Xuuu
y[7] = LXuuu
y[8] = Xuuuu
y[9] = LXuuuu
y[10] = Xuuuuu

ys=[]
zs=[]
count = 0
for g in range(0,len(ks)):
    k[1]=ks[g]
    parameters=(k,L)
    zs.append(ks[g])
    def f(y, t, k, L):
        return (
            (k[1]*y[1]-k[0]*y[0]*L),#dX0/dt
            (k[0]*y[0]*L-k[1]*y[1]-k[2]*y[1]),#dLX0/dt
            (k[2]*y[1]+k[4]*y[3]-k[3]*y[2]*L),#dXu/dt
            (k[3]*y[2]*L-k[4]*y[3]-k[5]*y[3]),#dLXu/dt
            (k[5]*y[3]+k[7]*y[5]-k[6]*y[4]*L),#dXuu/dt
            (k[6]*y[4]*L-k[7]*y[5]-k[8]*y[5]),#dLXuu/dt
            (k[8]*y[5]+k[10]*y[7]-k[9]*y[6]*L),#dXuuu/dt
            (k[9]*y[6]*L-k[10]*y[7]-k[11]*y[7]),#dLXuuu/dt
            (k[11]*y[7]+k[13]*y[9]-k[12]*y[8]*L),#dXuuuu/dt
            (k[12]*y[8]*L-k[13]*y[9]-k[12]*y[9]),#dLXuuuu/dt
            (k[14]*y[9])#dXuuuuu/dt
            )

    iC = []
    for n in range(0,len(y)):
        iC.append(0)
    for n in range(0,len(y)):
        iC[n]=y[n]

    t = range(0, 101,1)

    res = odeint(f, iC, t, parameters)

    transres = np.transpose(res)
    toappend = []
    toappend.append(transres[0][-1])
    toappend.append(transres[2][-1])
    toappend.append(transres[4][-1])
    toappend.append(transres[6][-1])
    toappend.append(transres[8][-1])
    toappend.append(transres[10][-1])

    ys.append(toappend)
    toappend=[]

scipy.io.savemat(filepath,mdict={'ys':ys})
