#import statements
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import odeint
import random

xData = np.loadtxt ("beta")
def virus(x,t,g,k):
    D = x[0]
    A = x[1]
    S = x[5]
    F1 = x[2]
    F2 = x[3]
    dDdt = (-g*D*A)
    dAdt = (-g*D*A) - (g*S*A)
    dF1dt = ((2*g*D*A)+(g*S*A)-(k*F1))
    dF2dt = (k*F1)-(k*F2)
    dSdt = (k*F2)
    return [dDdt, dAdt, dF1dt, dF2dt, dSdt]

#comparing to real data   
def func(p):
    t = xData[:, 0]
    g = p[0]
    k = p[2]
    ym = odeint(virus,[p[1], (100-p[1]), 0, 0, 0], t, args = (g, k))
    global yi
    yi = xData[:,1]
    SSR = sum((ym[:,-1]-yi)**2)
    print(SSR)
    return SSR  
full_output = 1
results = sp.optimize.minimize(func, [.0001, 50, 0.025])

# plots the true data
plt.scatter(xData[:,0], xData[:,1], s=1.5)
g = results.x[0]
t = xData[:, 0]
y = odeint(virus,[results.x[1], (100-results.x[1]), 0, 0, 0], t, args = (g, results.x[2]))
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.plot(t, y[:, -1])

#comparing to bootstrapping    
def fitguas(p):  
    ym = odeint(virus, [p[1], (100-p[1]), 0, 0 , 0], t, args = (p[0], p[2]))
    #yi = xData[:,-1]
    SSR = sum((ym[:,-1] - ynew)**2)
    return SSR

file=open("BS_beta", "a")
data = np.zeros([1000, 4])
r = (y[:,-1]-yi)
new_guess = results.x

for i in range(0,1000):
    random.shuffle(r)
    ynew = y[:,-1] + (r)
    results = sp.optimize.minimize(fitguas, new_guess, method='Nelder-Mead')
    print(results)
    data[i,:] = [results.fun, results.x[0], results.x[1], results.x[2]]
    np.savetxt("BS_beta", data)
file.close()

plt.xlabel("Time Post Transfection (h)", fontsize=18)
plt.ylabel("Fusion (%GFP + Pixels)", fontsize=17)
plt.title("Beta", fontsize=18)

#SSR Value: 
#x:
