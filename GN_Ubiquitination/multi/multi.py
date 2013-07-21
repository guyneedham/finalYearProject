# import modules
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import numpy as np
plt.rcParams['font.size']=20
plt.rcParams['lines.linewidth']=2.5
# Define the terms of the equations
# Y = [protein]
# t = time

# knumber = rate of that reaction
k=[]
for i in range(0,10):
    k.append(0)
    k[i]=1
k[2]=0.5
k[6]=0.5

if k[2] + k[6] != 1:
    print "k3 + k9 != 1. k3 + k9 must sum to 1."
    exit()

# L = [E3 ligase]
L = 1

# Parameters in a tuple
parameters = (k, L)

start = 5
y=[]
y0=start
LY0=0
Y11=0
LY11=0
Y21=0
LY21=0
Y2=0
for i in range(0,7):
    y.append(0)
y[0]=y0
y[1]=LY0
y[2]=Y11
y[3]=LY11
y[4]=Y21
y[5]=LY21
y[6]=Y2

# The derivative function:
def f(y, t, k, L):
    return (
        (k[1]*y[1]-k[0]*y[0]*L),#dY0/dt
        (k[0]*y[0]*L-k[1]*y[1]-k[2]*y[1]-k[6]*y[1]),#dLY0/dt
        (k[2]*y[1]+k[4]*y[3]-k[3]*y[2]*L),#dY11/dt
        (k[3]*L*y[2]-k[4]*y[3]-k[5]*y[3]),#dLy11/dt
        (k[6]*y[1]+k[8]*y[5]-k[7]*L*y[4]),#dY21/dt
        (k[7]*L*y[4]-k[8]*y[5]-k[9]*y[5]),#dLY21/dt
        (k[5]*y[3]+k[9]*y[5])#dY2/dt
            )

# Define the initial concentrations of the protein states
iC = []
for i in range(0,len(y)):
    iC.append(0)
    iC[i]=y[i]

# Times to evaluate the ODEs.
t = range(0, 41, 1)

# Compute the ODE:
res = odeint(f, iC, t, parameters)

# transpose the matrix
transres = np.transpose(res)
print "Plotting now."

# Plot the results
x = t
#for i in range(0,11):
#plt.plot(x, transres[i], label=str(i))
y1 = transres[0]
y2 = transres[2]
y3 = transres[4]
y4 = transres[6]
plt.plot(x, y1, label='Y$_0$', linestyle='-')
plt.plot(x, y2, label=r'Y$^{1}$$_{1}$',linestyle='--')
plt.plot(x, y3, label=r'Y$^{2}$$_{1}$',linestyle='-.')
plt.plot(x, y4, label='Y$_2$', linestyle='-',color='black')
plt.ylabel("Protein Concentration, A.U.")
plt.xlabel("Time, A.U.")
plt.ylim(0,6)
plt.xlim()
#plt.title("Mass Action Model of Multi-ubiquitination with Michaelis-Menten Kinetics")
plt.legend(loc=4)


plt.show()

