from binascii import a2b_hex
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings

# solve F(r) = (r-a)/c

class gwmodel():

    def __init__(self, rho, k, a, c):
        self.rho = rho      # threshold fraction
        self.k = k      # average degree
        self.a = a      # share of active
        self.c = c      # share of contingently active
        self.roots = []

    # distribution threshold
    # r is overall share of active nodes''
    def f(self, r):
        k = self.k
        rho = self.rho
        total_summation = 0

        for y in range(0, 2*k):
            inner_max = rho * y / (1 - rho) #float
            inner_max = int(inner_max)
            inner_sum = 0

            for x in range(0, inner_max):
                inner_sum = inner_sum + ((k*r)**x)/math.factorial(x)
        
            if(isinstance(inner_sum,np.ndarray)):
                inner_sum = inner_sum[0]

            inner_sum = float(inner_sum)
            outer_term = ((k-k*r)**y)/math.factorial(y)
            #print("uh:", outer_term)
            if(isinstance(outer_term,np.ndarray)):
                outer_term = outer_term[0]
            
            outer_term = float(outer_term)
            #temp = np.multiply(outer_term, inner_sum)
            temp = float(outer_term * inner_sum)
            #print("temp:", type(temp))
            #print("total:", type(total_summation))
            total_summation = total_summation + temp
        return 1 - (math.exp(-k) * total_summation)

    #line from (A,0) to (P,1)
    def line(self, r):
        return (r-self.a)/self.c

    def changeAC(self, a, c):
        self.a = a
        self.c = c

    # find fixed points
    def find_intersection(self, x):
        return self.f(x) - self.line(x)

    def plot_functions(self, root_x, root_y):

        fig = plt.figure()
        plt.grid(True)
        a = self.a
        c = self.c
        intervals = 5000
        x = np.linspace(0,1,intervals)
        a = [a] * intervals
        c = [c] * intervals

        distribution = list(map(self.f,x))
        line_ap = list(map(self.line, x))
        #difference = list(map(find_intersection, x, a, c))
        
        plt.plot(x, distribution, color = 'red')
        plt.plot(x, line_ap, color = 'blue')
        #plt.plot(x, difference, color = 'black')

        for i in range(0, len(root_x)):
            plt.plot(root_x[i],root_y[i], 'go')

        plt.title("F(t)")
        plt.xlabel('R(t)')
        plt.ylabel('R(t+1) = A + CF(R(t))')
        plt.show()

        return fig

    def solve_root(self, start_guess):
        root_x = fsolve(self.find_intersection, start_guess)
        return root_x

#   give the last iteration guess as the starting value 
    def guess(self, lastRoot):
        # deal with issues if only 1 fixed point
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            start_guess = [lastRoot, 0.5, 1-lastRoot]   # ************** fix
            try:
                root_x = self.solve_root(start_guess)
                if(root_x[0] == start_guess[0] or root_x[0] < 0):
                    start_guess = [0.2, 0.5, 0.8]
                    root_x = self.solve_root(start_guess)
                    if(root_x[0] == start_guess[0] or root_x[0] < 0):
                        start_guess = [0.3, 0.5, 0.7]
                        root_x = self.solve_root(start_guess)
            except:         # ********************* print out the error
                try:
                    start_guess = [0.9]
                    root_x = self.solve_root(start_guess)
                except:
                    start_guess = [0.1]
                    root_x = self.solve_root(start_guess)

            # sometimes duplicates when only 1 root
        if (len(root_x)>1):
            if(len(root_x) == 2):
                raise ValueError("not possible to have exactly 2 intersections")
            if(np.abs(root_x[0] - root_x[1]) < 0.00001):
                root_x = [np.mean(root_x)]

        # find correct y values, plot
        root_y = [self.line(x) for x in root_x]
        self.root_x = root_x
        #print("root(s): ", root_x)
        return root_x, root_y

    def roots(self):
        return self.root_x
