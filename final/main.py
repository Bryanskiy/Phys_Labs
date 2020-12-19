import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def mid(npArr):
    return npArr.sum() / len(npArr)

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

def b_sigma(x, y):
    return k_sigma(x, y) * (mid(x*x))**0.5

def r_line_point(x, y, a, b, c):
    return np.abs(a * x + b * y + c) / np.sqrt(a*a + b*b) 

def read_file(name):
    fst = np.loadtxt(name, dtype=np.float, usecols = 0)
    scd = np.loadtxt(name, dtype=np.float, usecols = 1)
    return fst, scd

def plot(x, y, x_label, y_label, plot_name):
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.plot(x, y, 'o', markersize=5)

    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.savefig(plot_name)
    return 

def lower_bound(A, key): 
    left = -1 
    right = len(A) 
    while right > left + 1: 
        middle = (left + right) // 2 
        if A[middle] >= key: 
            left = middle 
        else: 
            right = middle 
    return right

def callibrate(T, R, r):
    right = lower_bound(R, r)
    left = right - 1
    return T[left] + np.abs((T[right] - T[left]) / (R[right] - R[left]) * (r - R[left])) 
        

def find_points_in_temperature_range(T, R, left_border, right_border):
    T_ret = np.array([])
    R_ret = np.array([])

    for i in range(T.size):
        if (T[i] <= right_border) and (T[i] >= left_border):
            T_ret = np.append(T_ret, T[i])
            R_ret = np.append(R_ret, R[i])

    return T_ret, R_ret        

def linary_process(T, R, koeff):
    T_ret = np.array([])
    R_ret = np.array([])

    k = k_coef(T, R)
    b = b_coef(T, R)

    sigma_b = b_sigma(T, R)
    for i in range(T.size):
        if r_line_point(T[i], R[i], k, -1, b) < (sigma_b * koeff):
            T_ret = np.append(T_ret, T[i])
            R_ret = np.append(R_ret, R[i])
 
    return T_ret, R_ret           

def process_noize(T, R, normal_start, normal_end, noize_start, noize_end):
    T_ret = np.array([])
    R_ret = np.array([])

    T_tmp, R_tmp = find_points_in_temperature_range(T, R, normal_start, normal_end)

    k = k_coef(T_tmp, R_tmp)
    b = b_coef(T_tmp, R_tmp)
    sigma_b = b_sigma(T_tmp, R_tmp)

    T_full, R_full = find_points_in_temperature_range(T, R, np.minimum(normal_start, noize_start), np.maximum(normal_end, noize_end))

    for i in range(T_full.size):
        if r_line_point(T_full[i], R_full[i], k, -1, b) < (sigma_b * 0.5):
            T_ret = np.append(T_ret, T_full[i])
            R_ret = np.append(R_ret, R_full[i])

    return T_ret, R_ret 

def process(temperature, resistance):
    temperature_processed = np.array([])
    resistance_processed = np.array([])

    # tmp plot for illustretion
    T_tmp = np.array([])
    R_tmp = np.array([])
    for i in range(1700, 2300):    
        T_tmp = np.append(T_tmp, temperature[i])
        R_tmp = np.append(R_tmp, resistance[i])
    plot(T_tmp, R_tmp, "T, K", "R, Om", "processing_tmp.png")    
    #

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 0, 12.42)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 12.42, 12.56)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 12.56, 13.5)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 13.2, 16.5)
    for i in range(0, T0.size, 3):
        temperature_processed = np.append(temperature_processed, T0[i])
        resistance_processed = np.append(resistance_processed, R0[i])

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 13.5, 157)
    T0, R0 = linary_process(T0, R0, 15)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)    

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 170, 180)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 180, 190)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 190, 210)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0)

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 210, 241)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0) 

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 239, 300)
    T0, R0 = linary_process(T0, R0, 0.5)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0) 

    T0, R0 = find_points_in_temperature_range(temperature, resistance, 239, 243)
    for i in range(0, T0.size, 3):
        temperature_processed = np.append(temperature_processed, T0[i])
        resistance_processed = np.append(resistance_processed, R0[i])

    T0, R0 = process_noize(temperature, resistance, 100, 160, 160, 165)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0) 

    T0, R0 = process_noize(temperature, resistance, 170, 180, 165, 170)
    temperature_processed = np.append(temperature_processed, T0)
    resistance_processed = np.append(resistance_processed, R0) 

    return temperature_processed, resistance_processed

def main():
    I = 5 / 1000000

    # calibraton process
    T, R = read_file("calibration.dat")
    plot(T, R, "T, K", "R, Om", "calibration.png")

    # read main
    temperature, voltage = read_file("main.dat")
    for i in range(temperature.size):
        temperature[i] = callibrate(T, R, temperature[i])

    resistance = voltage / I    
    plot(temperature, resistance, "T, K", "R, Om", "without_processing.png")

    #processing
    temperature_processed ,resistance_processed = process(temperature, resistance)
    plot(temperature_processed, resistance_processed, "T, K", "R, Om", "processing.png")

main()    