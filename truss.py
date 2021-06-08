import numpy as np
import time
import math
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt
from truss4 import truss
from math import sin, cos, sqrt, pi


B = 10
#make array containing stress max's acording to problem
stressmax = np.ones(B)*25000
stressmax[8] = 75000

#outer most function so we can decouple the objective and constrains from the same function
def runopt(stressmax):
    #objective function history
    objhist = []

    #objective and constrain functino (decouple them)
    def objcon(x):
        nonlocal objhist
        mass, stress = truss(x)
        f = mass/25000 #scallin
        g = np.zeros(20)

        g[0:10] = (stressmax + stress) / 25000 #decompose absolute value constraint and scall
        g[10:20] = (stressmax - stress) / 25000 #decompose absolute value constraint and scall

        objhist.append(mass)
        return f, g #return mass and constraints. both should be of order 1 (.1 to 1)

#---------------------------- SAME
    xlast = []
    flast = []
    glast = []
    #obective functino
    def obj(x):
        nonlocal xlast, flast, glast
        if not np.array_equal(x, xlast):
            flast, glast = objcon(x)
            xlast = x
        return flast

    #constrain function
    def con(x):
        nonlocal xlast, flast, glast
        if not np.array_equal(x, xlast):
            flast, glast = objcon(x)
            xlast = x
        return glast

#----------------

    #array of initial guesses
    x0 = np.ones(B)*1

    #bounds, constraints, options for minimize function
    bnds = [(0.1, 20)] * 10
    cons = {'type': 'ineq', 'fun': con}
    options = {'disp': True}

    res = minimize(obj, x0, constraints=cons, bounds=bnds, options = options)
    res.fun = res.fun*100000 #scall solution

    #print results
    print("Stress chack", truss(res.x)[1])
    print("optimized cross sectional area matrix =", res.x, "inches squared")
    print("optimized function =", truss(res.x)[0], "pounds")



    return res.x, res.fun, objhist

if __name__ == '__main__':

    xstar, fstar, objhist = runopt(stressmax)
    plt.figure()
    plt.plot(objhist)
    plt.show()