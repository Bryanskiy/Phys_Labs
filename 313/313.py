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

    plt.plot(x, y, 'o', markersize=8)
    plt.plot(x_approx, y_approx)
    plt.errorbar(x, y, yerr = y_cross, xerr = x_cross, fmt = '.')

    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.savefig(plot_name)
    return 

def B_gor():
    N = 20
    n = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    t1 = np.array([8.97, 12.41, 19.15, 22.66, 29.72, 35.78, 49.25, 52.66, 57.78, 63.18])
    t2 = np.array([9.81, 13.21, 18.32, 21.31, 30.51, 37.32, 49.89, 52.84, 56.91, 62.87])
    T1 = t1 / N
    T2 = t2 / N

    T = (T1 + T2) / 2

    sigma_T_rand = np.array([])
    for i in range(10):
        tmp = sigma_rand([T1[i], T2[i]])
        sigma_T_rand = np.append(sigma_T_rand, tmp)

    B_gor_k = k_only(n, T)
    B_gor_k_sigma = k_only_sigma(n, T)

    x_approx = np.linspace(0, 12, 13)
    y_approx = B_gor_k * x_approx

    x_crosses = np.zeros(10)
    y_crosses = np.array([])
    for i in range(x_crosses.size):
        tmp = np.sqrt(sigma_T_rand[i]**2 + 0.26**2)
        y_crosses = np.append(y_crosses, tmp)
    
    plot(n, T, "T", x_approx, y_approx, x_crosses, y_crosses)

    return B_gor_k, B_gor_k_sigma

def B_vert():
    n = np.array([4, 6, 8, 10, 12])
    m = np.array([0.256,  0.198, 0.564, 0.378, 0.924])
    l_tmp = np.array([1., 2., 1., 2., 1.])
    l = 0.568 * l_tmp

    l_sigma = np.array([])
    for i in range(l.size):
        tmp = sigma_coss(l[i], [0.031], [0.568], [1])
        l_sigma = np.append(l_sigma, tmp)
 
    m_sigma = np.array([0.001, 0.001, 0.001, 0.001, 0.001])
    M = l * m * 9.81 * 100
    y_crosses = np.array([])
    for i in range(M.size):
        tmp = sigma_coss(M[i], [m_sigma[i], l_sigma[i]], [m[i], l[i]], [1, 1])
        y_crosses = np.append(y_crosses, tmp)

    x_crosses = np.zeros(5)

    B_vert_k = k_only(n, M)
    B_vert_k_sigma = k_only_sigma(n, M)

    x_approx = np.linspace(0, 12, 13)
    y_approx = B_vert_k * x_approx

    plot(n, M,"M", x_approx, y_approx, x_crosses, y_crosses)

    return B_vert_k, B_vert_k_sigma


def main():
    #mass, g
    tara_mass = 14.534
    masses = np.array([0.833, 0.835, 0.838, 0.846, 0.830])
    mid_mass = mid(masses)
    random_sigma_masses = sigma_rand(masses)
    sigma_mass = np.sqrt(random_sigma_masses**2 + 0.001**2)

    print("mid mass = {0} +- {1}, g".format(mid_mass, sigma_mass))

    #diams, sm
    diams = np.array([0.57, 0.57, 0.55, 0.57, 0.58])
    mid_diam = mid(diams)
    random_sigma_diams = sigma_rand(diams)
    sigma_diam = np.sqrt(random_sigma_diams**2 + 0.01**2)

    print("mid diam = {0} +- {1}, sm".format(mid_diam, sigma_diam))
    print("-------------------------------------------")

    #A method
    R_max = 2.8
    sigma_R_max = 0.1

    Pm_A = 10 * R_max**2 * np.sqrt(mid_mass * 9.81 / 6)
    sigma_Pm_A = sigma_coss(Pm_A, [sigma_R_max, sigma_mass], [R_max, mid_mass], [2, 0.5])

    print("Pm_A = {0} +- {1}".format(Pm_A, sigma_Pm_A))

    p_m = Pm_A * 6 / (3.14 * mid_diam**3)
    sigma_p_m = sigma_coss(p_m, [sigma_Pm_A, sigma_diam], [Pm_A, mid_diam], [1, 3])

    print("pm_A = {0} +- {1}".format(p_m, sigma_p_m))

    B_p_A = 2 * Pm_A / (mid_diam / 2)**3
    sigma_B_p_A = sigma_coss(B_p_A, [sigma_Pm_A, sigma_diam], [Pm_A, mid_diam], [1, 3])

    print("B_p_A = {0} +- {1}".format(B_p_A, sigma_B_p_A))

    B_r = 4 * 3.14 * p_m
    sigma_B_r = sigma_coss(B_r, [sigma_p_m], [p_m], [1])

    print("B_r = {0} +- {1}".format(B_r, sigma_B_r))

    #B method
    mass_sum = 500 + 35.523 - tara_mass;
    Pm_B = 10 * mid_diam**2 * np.sqrt((mass_sum) * 9.81 / 1.08 / 6)
    sigma_mass_sum = np.sqrt(10**2 + 0.001**2 + 0.001**2 + 0.836**2)
    sigma_Pm_B = sigma_coss(Pm_B, [sigma_diam, sigma_mass_sum], [mid_diam, mass_sum], [2, 0.5])

    print("-----------------------------------------")
    print("mass_sum = {0} +- {1}".format(mass_sum, sigma_mass_sum))
    print("Pm_B = {0} +- {1}".format(Pm_B, sigma_Pm_B))

    B_p_B = 16 * Pm_B / (mid_diam)**3
    sigma_B_p_B = sigma_coss(B_p_B, [sigma_Pm_B, sigma_diam], [Pm_B, mid_diam], [1, 3])

    print("B_p_B = {0} +- {1}".format(B_p_B, sigma_B_p_B))


    # B gor
    B_gor_k, B_gor_k_sigma = B_gor()
    B_gor_value = 3.14**2 * mid_mass * mid_diam**2 / 3 / B_gor_k**2 / Pm_B
    B_gor_sigma = sigma_coss(B_gor_value, [sigma_mass, sigma_diam, B_gor_k_sigma, sigma_Pm_B], [mid_mass, mid_diam, B_gor_k, Pm_B], [1, 2, 2, 1])

    print("-----------------------------------------")
    print("B_gor_k = {0} +- {1}".format(B_gor_k, B_gor_k_sigma))
    print("B_gor_value = {0} +- {1}".format(B_gor_value, B_gor_sigma))
    print("-----------------------------------------")

    # B vert
    B_vert_k, B_vert_k_sigma = B_vert()
    print("B_vert_k = {0} +- {1}".format(B_vert_k, B_vert_k_sigma))
    B_vert_value = B_vert_k / Pm_A
    B_vert_sigma = sigma_coss(B_vert_value, [B_vert_k_sigma, sigma_Pm_A], [B_vert_k, Pm_A], [1, 1])
    print("B_vert_value = {0} +- {1}".format(B_vert_value, B_vert_sigma))

    print("-----------------------------------------")
    B = np.sqrt(B_gor_value**2 + B_vert_value**2)
    B_sigma = 1 * (2 * B_gor_value * B_gor_sigma + 2 * B_vert_value * B_vert_sigma) / (2 * B)
    print("B = {0} +- {1}".format(B, B_sigma))

    tg_b = B_vert_value / B_gor_value
    tg_b_sigma = sigma_coss(tg_b, [B_vert_sigma, B_gor_sigma], [B_vert_value, B_gor_value], [1, 1])
    print("tg_b = {0} +- {1}".format(tg_b, tg_b_sigma))

    return 

main()