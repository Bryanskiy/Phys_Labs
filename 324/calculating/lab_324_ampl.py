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

def plot(x, y, plot_name, x_approx = np.array([]), y_approx = np.array([]), x_cross = np.array([]), y_cross = np.array([])) :
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))

    plt.ylabel("1/teta^2")
    plt.xlabel("1/R^2, 1 / Om^2")

    plt.plot(x, y, 'o', markersize=8)

    if(x_approx.size > 0 and y_approx.size > 0):
        plt.plot(x_approx, y_approx)

    if(x_cross.size > 0 and y_cross.size > 0):    
        plt.errorbar(x, y, yerr = y_cross, xerr = x_cross, fmt = '.')

    plt.legend()
    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.savefig(plot_name)
    return 

def main():
    x = np.array([1.0000, 0.6944, 0.4444, 0.3086, 0.2500, 0.1600, 0.1111])
    y = np.array([3.7351, 2.6166, 1.8794, 1.2272, 1.0564, 0.4294, 0.1539])

    x_approx = np.linspace(x[0], x[x.size - 1], x.size)
    y_approx = k_only(x, y) * x_approx

    sigma_R = 0.02
    R_arr = np.array([1, 1.2, 1.5, 1.8, 2, 2.5, 3])
    cross_x = np.array([])
    for i in xrange(x.size):
        tmp = sigma_coss(x[i], [sigma_R], [R_arr[i]], [2])
        cross_x = np.append(cross_x, tmp)


    u_k = np.array([8.5, 8.8, 3.7, 3, 2.8, 2.3, 3.2])
    u_kn = np.array([1.8, 0.4, 0.2, 0.2, 0.4, 0.5, 0.25])
    n = np.array([3, 5, 4, 3, 2, 1, 1])
    div = u_k / u_kn
    log = np.log(div)
    teta = log / n

    div_sigma = np.array([])
    for i in xrange(div.size):
        tmp = sigma_coss(div[i], [0.1, 0.1], [u_k[i], u_kn[i]], [1, 1])
        div_sigma = np.append(div_sigma, tmp)

    log_sigma = sigma_log(div, div_sigma)

    sigma_teta = np.array([])
    for i in xrange(teta.size):
        tmp = sigma_coss(teta[i], [log_sigma[i]], [log[i]], [1])
        sigma_teta = np.append(sigma_teta, tmp)

    cross_y = np.array([])
    for i in xrange(y.size):
        tmp = sigma_coss(y[i], [sigma_teta[i]], [teta[i]], [1])
        cross_y = np.append(cross_y , tmp)

    plot(x, y, "ampl", x_approx = x_approx, y_approx = y_approx, x_cross=cross_x, y_cross=cross_y )
    return 


main()