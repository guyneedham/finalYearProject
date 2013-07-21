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
k = []
for a in range(0,15):
    k.append(0)
for a in range(0,15):
    k[a]=1
k[1]=1
k[3]=1
k[6]=1
k[9]=1
k[13]=1
#print "k = ",k, len(k)

L = 1
# L = [E3 ligase]
start=5
parameters = (k, L)
iX0 = start
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

def f(y, t, k, L):
     return (
         (k[1]*y[1]-k[0]*y[0]*L),#dX0/dt
         (k[0]*y[0]*L-k[1]*y[1]-k[2]*y[1]),#dLX0/dt
         (k[2]*y[1]+k[3]*y[3]-k[4]*y[2]*L),#dXu/dt
         (k[4]*y[2]*L-k[3]*y[3]-k[5]*y[3]),#dLXu/dt
         (k[5]*y[3]+k[6]*y[5]-k[7]*y[4]*L),#dXuu/dt
         (k[7]*y[4]*L-k[6]*y[5]-k[8]*y[5]),#dLXuu/dt
         (k[8]*y[5]+k[9]*y[7]-k[10]*y[6]*L),#dXuuu/dt
         (k[10]*y[6]*L-k[9]*y[7]-k[11]*y[7]),#dLXuuu/dt
         (k[11]*y[7]+k[13]*y[9]-k[12]*y[8]*L),#dXuuuu/dt
         (k[12]*y[8]*L-k[13]*y[9]-k[14]*y[9]),#dLXuuuu/dt
         (k[14]*y[9])#dXuuuuu/dt
         )

iC = []
for n in range(0,len(y)):
    iC.append(0)
for n in range(0,len(y)):
    iC[n]=y[n]

t = range(0, 41,1)

res = odeint(f, iC, t, parameters)

transres = np.transpose(res)
'''
#Plot as a bar chart

ys=[]
for i in range(0,6):
    ys.append(0)
ys[0]=(transres[0][-1]/start)*100
ys[1]=(transres[2][-1]/start)*100
ys[2]=(transres[4][-1]/start)*100
ys[3]=(transres[6][-1]/start)*100
ys[4]=(transres[8][-1]/start)*100
ys[5]=(transres[10][-1]/start)*100
N=6
width=0.5
ind=np.arange(N)
fig=plt.figure()
ax=fig.add_subplot(111)
rects1=ax.bar(ind,ys,width)
ax.set_ylabel('Percentage of Total Protein')
ax.set_title(r'Predicted Protein Percentages from a Model of Distributive Polyubiquitination where k$_1$ = 3000')
ax.set_xticks(ind+(width/2))
ax.set_xticklabels((r'X$_0$','X$_u$','X$_{uu}$','X$_{uuu}$','X$_{uuuu}$','X$_{uuuuu}$'))
ax.legend()
'''

#Plot as a line graph
x=t

y0 = transres[0]
y2 = transres[2]
y4 = transres[4]
y6 = transres[6]
y8 = transres[8]
y10 = transres[10]

plt.plot(x, y0, label=r'X$_0$')
plt.plot(x, y2, label=r'X$_u$')
plt.plot(x, y4, label=r'X$_{uu}$')
plt.plot(x, y6, label=r'X$_{uuu}$')
plt.plot(x, y8, label=r'X$_{uuuu}$')
plt.plot(x, y10, label=r'X$_{uuuuu}$',color='black')

plt.legend(loc=4)
#plt.title("Mass Action model of Distributive Polyubiquitination with Michaelis Menten Kinetics")
plt.ylabel("Protein Concentration, A.U.")
plt.xlabel("Time, A.U.")
plt.ylim(0,6)


plt.show()
