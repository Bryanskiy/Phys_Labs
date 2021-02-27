import numpy as np

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def mid(npArr):
    return npArr.sum() / len(npArr)

def sigma_coss(val, sigsArr, parmsArr, powerArr):
    sum = 0
    for i in range(len(sigsArr)):
        sum += (powerArr[i] * sigsArr[i] / parmsArr[i])**2
    return val * (sum)**0.5

def sigma_log(val, sigma_val):
    return sigma_val / val

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


def main():
    # ----------------FRENEL------------------------------------------ 
    _lambda = 546.1 # nano
    D = 216 ## mikro
    m = np.array([1, 2, 3, 4, 5, 6])
    a = np.array([68, 66, 63, 60, 55, 48]) # mili 
    _2zm = 2 * np.sqrt(m * a * _lambda * 10**(-6)) #mili
    sigma_2zm = np.array([])
    for i in range(_2zm.size):
        sigma_2zm = np.append(sigma_2zm, sigma_coss(_2zm[i], [2], [a[i]], [0.5]))

    print (sigma_2zm)

    ## PLOT
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))
    plt.xlabel("m")
    plt.plot(m, _2zm, 'o' ,label = "width of Frenel zones, mm", markersize=8)
    D_arr = [D/(10**(3))] * 6
    plt.plot(m, D_arr, label = "width of gap, mm")
    plt.errorbar(m, _2zm, yerr = sigma_2zm, fmt = '.')
    plt.minorticks_on()
    plt.grid(which='major', axis='both')
    plt.grid(which='minor', axis='both', linestyle = ':')
    plt.legend()
    plt.savefig("frenel")

    # ----------------FRAINH------------------------------------------ 
    F2 = 12.8 # santi
    m = np.array([-4, -3, -2, -1, 1, 2, 3, 4])
    xm = np.array([-0.83, -0.56, -0.35, -0.20, 0.22, 0.38, 0.55, 0.80])
    xm_error = np.array([0.01] * 8)

    k = k_coef(m, xm)
    b = b_coef(m, xm)
    sigma_k = k_sigma(m, xm)

    x_approx = np.linspace(-4, 4, 1000)
    y_approx = x_approx * k + b

    ## PLOT
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))
    plt.xlabel("m")
    plt.ylabel("x, mm")    
    plt.plot(m, xm, 'o', markersize=8)
    plt.plot(x_approx, y_approx, 'o', markersize=1)
    plt.errorbar(m, xm, yerr = xm_error, fmt = '.')
    plt.minorticks_on()
    plt.grid(which='major', axis='both')
    plt.grid(which='minor', axis='both', linestyle = ':')
    plt.savefig("frainh")


main()    