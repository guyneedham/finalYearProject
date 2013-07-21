# import modules
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import numpy as np
plt.rcParams['font.size']=20
plt.rcParams['lines.linewidth']=2.5
# Define the terms of the equations
# X = [protein]
# t = time

# knumber = rate of that reaction
kon = 1
koff = []
for n in range(0,6):
    koff.append(1)
k = []
for i in range(0,5):
    k.append(0)
#formation of XuL
k[0] = 1
#formation of XuuL
k[1] = 1
#formation of XuuuL
k[2] = 1
#formation of XuuuuL
k[3] = 1
#formation of XuuuuuL
k[4] = 1

# L = [E3 ligase]
L = 1

# Parameters, in a tuple
parameters = (kon, koff, k, L)
start = 5
# iX0 = initial [unubiquitinated protein]
iX0 = start
# iLX0 = initial [complex LX0]
iLX0 = 0
# iX1uL = initial [complex XuL] 1L
iX1uL = 0
# X1u = initial [complex Xu] 1
iX1u = 0
# X2uL = initial [complex XuuL] 2L
iX2uL = 0
# X2u = initial [Xuu] 2
iX2u = 0
# iX3uL = initial [XuuuL] 3L
iX3uL = 0
# iX3u = initial [Xuuu] 3
iX3u = 0
# iX4uL = initial [XuuuuL] 4L
iX4uL = 0
# iX4u = initial [Xuuuu] 4
iX4u = 0
# iX5uL = initial [XuuuuuL] 5L
iX5uL = 0
# iX5u = initial [Xuuuuu] 5
iX5u = 0

y = []
for n in range(0,12):
    y.append(0)
y[0] = iX0
y[1] = iLX0
y[2] = iX1uL
y[3] = iX1u
y[4] = iX2uL
y[5] = iX2u
y[6] = iX3uL
y[7] = iX3u
y[8] = iX4uL
y[9] = iX4u
y[10] = iX5uL
y[11] = iX5u

# The derivative function:
def f(y, t, kon, koff, k, L):
    return (
        (koff[0]*y[1]-kon*y[0]*L),#dX0/dt
            (kon*y[0]*L-koff[0]*y[1]-k[0]*y[1]),#dLX0/dt
            (k[0]*y[1]-koff[1]*y[2]-k[1]*y[2]),#dXuL/dt
            (koff[1]*y[2]),#dXu/dt
            (k[1]*y[2]-koff[2]*y[4]-k[2]*y[4]),#dXuuL/dt
            (koff[2]*y[4]),#dXuu/dt
            (k[2]*y[4]-koff[3]*y[6]-k[3]*y[6]),#dXuuuL/dt
            (koff[3]*y[6]),#dXuuu/dt
            (k[3]*y[6]-koff[4]*y[8]-k[4]*y[8]),#dXuuuuL/dt
            (koff[4]*y[8]),#dXuuuu/dt
            (k[4]*y[8]-koff[5]*y[10]),#dXuuuuuL/dt
            (koff[5]*y[10])#dXuuuuu/dt
            )
 
# Define the initial concentrations of the protein states

iC = []
for i in range(0,len(y)):
    iC.append(0)
for i in range(0,len(y)):
    iC[i]=y[i]
# Times to evaluate the ODEs.
t = range(0, 41, 1)

# Compute the ODE:
res = odeint(f, iC, t, parameters)

# transpose the matrix
transres = np.transpose(res)
#print transres
'''
#Plot as a bar chart

ys=[]
for i in range(0,6):
    ys.append(0)
ys[0]=(transres[0][-1]/start)*100
ys[1]=(transres[3][-1]/start)*100
ys[2]=(transres[5][-1]/start)*100
ys[3]=(transres[7][-1]/start)*100
ys[4]=(transres[9][-1]/start)*100
ys[5]=(transres[11][-1]/start)*100
N=6
width=0.5
ind=np.arange(N)
fig=plt.figure()
ax=fig.add_subplot(111)
rects1=ax.bar(ind,ys,width)
ax.set_ylabel('Percentage of Total Protein')
ax.set_title(r'Predicted Protein Percentages from a Model of Processive Polyubiquitination where k$_{off}$ = 1')
ax.set_xticks(ind+(width/2))
ax.set_xticklabels((r'X$_0$','X$_u$','X$_{uu}$','X$_{uuu}$','X$_{uuuu}$','X$_{uuuuu}$'))
ax.set_ylim(0,100)
ax.legend()

'''
#plot as a line graph

x = t
y =[]
for i in range(0,len(iC)):
             y.append(0)
for i in range(0,len(iC)):
             y[i]=transres[i]
plt.plot(x,y[0],label=r'X$_0$')
plt.plot(x,y[3],label=r'X$_u$')
plt.plot(x,y[5],label=r'X$_{uu}$')
plt.plot(x,y[7],label=r'X$_{uuu}$')
plt.plot(x,y[9],label=r'X$_{uuuu}$')
plt.plot(x,y[11],label=r'X$_{uuuuu}$',color='black')
plt.legend()
#plt.title("Mass Action model of Processive Polyubiquitination with Michaelis Menten Kinetics")
plt.ylabel("Protein Concentration, A.U.")
plt.xlabel("Time, A.U.")
plt.ylim(0,6)

#'''
plt.show()
