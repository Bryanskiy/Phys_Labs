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

def plot(x, y, plot_name, x_approx = [], y_approx = [], x_cross = [], y_cross = []):
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))

    for i in range(len(y)):
        plt.plot(x, y[i], 'o', markersize=6)

    for i in range(len(y_approx)):
        plt.plot(x_approx, y_approx[i])

    for i in range(len(y_cross)):
        plt.errorbar(x, y[i], yerr = y_cross[i], xerr = x_cross, fmt = '.')  

    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.savefig(plot_name)
    return 


def grad_electormagnet():
    I = np.array([0.0, 0.23, 0.53, 0.71, 0.90, 1.11, 1.38]) #A
    Wb = np.array([0.1, 2.3, 4.4, 6.0, 7.4, 8.9, 10.9]) # mWb
    B = Wb / 75 * 10**4 # mTl

    x_approx = np.linspace(0.0, 1.4, 20)
    y_approx = k_only(I, B) * x_approx

    x_crosses = np.array(7 * [0.01])
    y_crosses = np.array([])
    for i in range(B.size):
        tmp = sigma_coss(B[i], [0.1], [Wb[i]], [1])
        y_crosses = np.append(y_crosses, tmp)

    plot(I, [B], "grad_electromagnet", x_approx, [y_approx], x_crosses, [y_crosses])

    print k_only(I, B)
    print k_only_sigma(I, B)

    return k_only(I, B)
    
def Holl():
    I = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
    I_approx = np.linspace(0, 1.1, 100)
    I_crosses = np.array(6 * [0.01])
    U_crosses = np.array(6 * [0.008])

    # U_1
    #----------------------------------

    I_0_1 = 0.34 # mA
    U_0_1 = 0.050 # mB

    U_1 = np.array([0.038, 0.010, -0.015, -0.038, -0.057, -0.065])
    U_1 = U_1 - U_0_1

    U_1_approx = k_coef(I, U_1) * I_approx + b_coef(I, U_1)
    
    # U_2
    #----------------------------------

    I_0_2 = 0.43 #mA
    U_0_2 = 0.065 # mB

    U_2 = np.array([0.050, 0.013, -0.019, -0.047, -0.072, -0.082])
    U_2 = U_2 - U_0_2

    U_2_approx = k_coef(I, U_2) * I_approx + b_coef(I, U_2)

    # U_3
    #----------------------------------

    I_0_3 = 0.53
    U_0_3 = 0.082

    U_3 = np.array([0.062, 0.018, -0.022, -0.057, -0.088, -0.100])
    U_3 = U_3 - U_0_3

    U_3_approx = k_coef(I, U_3) * I_approx + b_coef(I, U_3)

    # U_4
    #----------------------------------

    I_0_4 = 0.69
    U_0_4 = 0.107

    U_4 = np.array([0.082, 0.026, -0.034, -0.074, -0.115, -0.131])
    U_4 = U_4 - U_0_4

    U_4_approx = k_coef(I, U_4) * I_approx + b_coef(I, U_4)

    # U_5
    #----------------------------------

    I_0_5 = 0.82
    U_0_5 = 0.127

    U_5 = np.array([0.098, 0.031, -0.034, -0.088, -0.133, -0.153])
    U_5 = U_5 - U_0_5

    U_5_approx = k_coef(I, U_5) * I_approx + b_coef(I, U_5)

    # U_6
    #----------------------------------

    I_0_6 = 0.96
    U_0_6 = 0.150

    U_6 = np.array([0.114, 0.035, -0.037, -0.103, -0.158, -0.182])
    U_6 = U_6 - U_0_6

    U_6_approx = k_coef(I, U_6) * I_approx + b_coef(I, U_6)

    # conclude 
    plot(I, [U_1, U_2, U_3, U_4, U_5, U_6], "holl", I_approx, [U_1_approx, U_2_approx, U_3_approx, U_4_approx, U_5_approx, U_6_approx], 
    I_crosses, 6 * [U_crosses])

    return [np.array([k_coef(I, U_1), k_coef(I, U_2), k_coef(I, U_3), k_coef(I, U_4), k_coef(I, U_5), k_coef(I, U_6)]), 
    [k_sigma(I, U_1), k_sigma(I, U_2), k_sigma(I, U_3), k_sigma(I, U_4), k_sigma(I, U_5), k_sigma(I, U_6)]]

def general(holl_coefs, holl_sigmas):
    I = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
    holl_coefs_np = np.array(holl_coefs)
    x_approx = np.linspace(0.0, 1.1, 100)
    y_approx = k_coef(I, holl_coefs_np) * x_approx + b_coef(I, holl_coefs_np)

    I_crosses = np.array(6 * [0.01])
    coefs_crosses = np.array(holl_sigmas)

    plot(I, [holl_coefs_np], "general", x_approx, [y_approx], I_crosses, [coefs_crosses])

    return k_coef(I, holl_coefs_np)

def main():
    k_calibration = grad_electormagnet()
    tmp = Holl()
    holl_coefs = tmp[0] / 1.1
    holl_sigmas = tmp[1]
    general_coef = general(holl_coefs, holl_sigmas)

    Rx = -general_coef * 2.2 * 10.0**(-3)
    print Rx
main()