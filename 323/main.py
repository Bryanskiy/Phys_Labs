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

def plot(x, y, x_crosses, y_crosses, labels, plot_name):
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))

    plt.xlabel("x, mm")
    plt.ylabel("I, mA")

    for i in range(len(y)):
        plt.plot(x, y[i], 'o', markersize=8, label = labels[i])

    for i in range(len(y_crosses)):
        plt.errorbar(x, y[i], yerr = y_crosses[i], xerr = x_crosses, fmt = '.')    

    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.legend()
    plt.savefig(plot_name)
    return

def main():
    x = np.array([95, 90, 85, 80, 75, 71, 65, 60, 55, 52, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0]) # mm
    I = np.array([200, 150, 110 ,75, 50, 25, 45, 60, 75, 100, 130, 150, 167.5, 185, 200, 212, 230, 255, 275, 290]) #mA
    I_l = np.array([670, 620, 570, 520, 470, 425, 385, 360, 340, 320, 280, 260, 230, 215, 200, 170, 135, 100, 75, 50]) #mA
    I_c = np.array([410, 405, 410, 410, 405, 420, 405, 405, 410 ,410, 405, 410, 405, 400, 400, 390, 400, 410, 410, 400]) #mA

    midic = mid(I_c)
    print midic

    x_crosses = np.array(20 * [1])
    I_crosses = np.array(20 * [2.5])
    I_l_crosses = np.array(20 * [20])
    I_c_crosses = I_l_crosses

    plot(x, [I, I_l, I_c], x_crosses, [I_crosses, I_l_crosses, I_c_crosses], ["I", "Il", "Ic"],"Ampr")


main()   