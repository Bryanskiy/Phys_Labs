import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

def mid(npArr):
    return npArr.sum() / len(npArr)

def sigma_coss(val, sigsArr, parmsArr, powerArr):
    sum = 0
    for i in range(len(sigsArr)):
        sum += (powerArr[i] * sigsArr[i] / parmsArr[i])**2
    return val * (sum)**0.5

def sigma_rand(valsArr):
    N = len(valsArr)
    valMid = sum(valsArr) / N
    msum = 0
    for val in valsArr:
        msum += (val - valMid)**2
        
    return (msum / (N * (N - 1))) ** 0.5

def k_only(x, y):
    return mid(x * y) / mid(x * x)

def k_only_sigma(x, y):
    return ((mid(x * x) * mid(y * y) - mid(x * y)**2) / (len(x) * mid(x * x)**2))**0.5

def k_coef(x, y):
    return (mid(x * y) - mid(x) * mid(y)) / (mid(x * x) - mid(x)**2)

def b_coef(x, y):
    return mid(y) - k_coef(x, y) * mid(x)

def k_sigma(x, y):
    return ((mid(y * y) - mid(y)**2) / (mid(x * x) - mid(x)**2) - k_coef(x, y)**2)**0.5  / len(x)**0.5

def plot(x, y, plot_name, x_approx = [], y_approx = [], x_cross = [], y_cross = []) :
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))

    plt.ylabel("T exp, mC")
    plt.xlabel("T teor, mC")

    plt.plot(x, y, 'o', markersize=8)
    plt.plot(x_approx, y_approx)
    plt.errorbar(x, y, yerr = y_cross, xerr = x_cross, fmt = '.')

    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.savefig(plot_name)
    return 


def main():
    x = np.array([0.4, 0.63, 0.79, 0.89, 1.09, 1.26, 1.54, 1.78])
    y = np.array([0.44, 0.68, 0.83, 0.98, 1.11, 1.25, 1.6, 2])

    x_approx = np.linspace(x[0], x[x.size - 1], x.size)
    y_approx = k_only(x, y) * x_approx

    delta_L = sigma_rand([203.5, 199.0, 199.5])
    x_cross = np.array([])
    for elem in x:
        tmp = sigma_coss(elem, [delta_L], [200], [0.5])
        x_cross = np.append(x_cross, tmp)
    
    y_cross = np.array([])
    coefs = np.array([1.1, 3.4, 5, 3.9, 5, 5, 4, 4.2])
    for i in xrange(y.size):
        tmp = sigma_coss(y[i], [0.2] , [coefs[i]], [1])
        y_cross = np.append(y_cross, tmp)

    plot(x, y, "TT", x_approx, y_approx, x_cross, y_cross)
    return;

main()