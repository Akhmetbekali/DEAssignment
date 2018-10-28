import numpy as np
from matplotlib import pyplot as plt

# Reading User input
xi = float(input("Enter x0: "))
yi = float(input("Enter yi: "))
Xf = float(input("Enter X: "))
num = int(input("Enter n: "))

# Numerical Methods
def euler(x0, y0, X, n):
    h = (X - x0) / n # Calculating step
    x = np.linspace(x0, X, n + 1) # Creating an array of x
    y = np.zeros(n + 1)
    y[0] = y0
    for i in range(1, n + 1):
        y[i] = h * (np.cos(x[i-1]) - (y[i-1])) + y[i-1] # Calculating value of y for each x[i]
    plt.subplot(211)
    plt.title("Numerical Methods")
    plt.plot(x,y, label = 'Euler Method')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 5.5}, borderaxespad=0.)
    return x,y

def improved_euler(x0, y0, X, n): # 1, 1, 10, 10
    h = (X - x0) / n
    x = np.linspace(x0, X, n + 1)
    y = np.zeros(n + 1) # Create an array with zeros in each position. Zero, because we need to calculate y
    y[0] = y0
    for i in range(1, n + 1): # Calculating y[i] according to formula
        k1 = np.cos(x[i-1]) - (y[i-1])
        k2 = np.cos(x[i]) - (y[i-1] + h*k1)
        y[i] = h * ((k1+k2)/2) + y[i-1]
    plt.subplot(211)
    plt.plot(x, y, label='Improved Euler Method')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size':5.5}, borderaxespad=0.)
    return x,y

def runge_kutta(x0, y0, X, n):
    h = (X - x0) / n
    x = np.linspace(x0, X, n + 1)
    y = np.zeros(n + 1)
    y[0] = y0
    for i in range(1, n + 1): # Calculating each k according to formula and calculating y[i]
        k1 = np.cos(x[i - 1]) - (y[i - 1])
        k2 = np.cos(x[i - 1] + h/2) - (y[i-1] + (h/2)*k1)
        k3 = np.cos(x[i - 1] + h / 2) - (y[i - 1] + (h / 2) * k2)
        k4 = np.cos(x[i - 1] + h) - (y[i - 1] + h * k3)
        y[i] = h/6 * (k1 + 2*k2 + 2*k3 + k4) + y[i-1]
    plt.subplot(211)
    plt.plot(x, y, label='Runge Kutta Method')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 5.5}, borderaxespad=0.)
    return x,y

def IVP(x0,y0, X, n):
    c = ((np.exp(x0))*(2*y0 - np.cos(x0) - np.sin(x0)))/2
    x = np.linspace(x0, X, n + 1)
    y = np.zeros(n + 1)
    y[0] = y0
    for i in range(1, n + 1):
        y[i] = ((np.cos(x[i]) + np.sin(x[i]))/2) + (c / np.exp(x[i]))
    plt.subplot(211)
    plt.plot(x, y, label = 'Exact solution of the Initial Value Problem')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 5.5}, borderaxespad=0.)
    return x,y

# Computation of errors
def error_euler(euler_x, euler_y, exact_x, exact_y):
    if euler_x.all() == exact_x.all(): # Checking for the same x coordinates
        error_x = exact_x
    error_y = exact_y - euler_y # Calculating error in y

    plt.subplot(212)
    plt.title("Graph of errors")
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.plot(error_x, error_y, label = 'Euler Method Error')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 5.5}, borderaxespad=0.)

def error_improved_euler(ieuler_x, ieuler_y, exact_x, exact_y):
    if ieuler_x.all() == exact_x.all():
        error_x = exact_x
    error_y = exact_y - ieuler_y
    plt.subplot(212)
    plt.plot(error_x, error_y, label = 'Improved Euler Method Error')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 5.5}, borderaxespad=0.)

def error_rk(rk_x, rk_y, exact_x, exact_y):
    if rk_x.all() == exact_x.all():
        error_x = exact_x
    error_y = exact_y - rk_y
    plt.subplot(212)
    plt.plot(error_x, error_y, label = 'Runge Kutta Method Error')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 5.5}, borderaxespad=0.)

# Displaying of Methods and Errors
ivp_x, ivp_y = IVP(xi, yi, Xf, num)
e_x, e_y = euler(xi, yi, Xf, num)
ie_x, ie_y = improved_euler(xi, yi, Xf, num)
rk_x, rk_y = runge_kutta(xi, yi, Xf, num)

error_euler(e_x, e_y, ivp_x, ivp_y)
error_improved_euler(ie_x, ie_y, ivp_x, ivp_y)
error_rk(rk_x, rk_y, ivp_x, ivp_y)

plt.show()
