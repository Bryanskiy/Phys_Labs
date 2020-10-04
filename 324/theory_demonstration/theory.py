import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plots(x, y_arr, x_label, y_labels_arr, plot_name):
    mpl.rcParams['font.size'] = 16
    plt.figure(figsize=(16,9))

    plt.xlabel(x_label)

    for i in xrange(len(y_arr)):
        plt.plot(x, y_arr[i], label = y_labels_arr[i])

    plt.legend()
    plt.grid(b=True, which='major', axis='both', alpha=5)

    plt.savefig(plot_name)
    return


def aperiodic_oscillations_voltage(t, C, L, R, voltage_0 = 1):
    gamma = R / (2.0 * L)
    omega_0 = np.sqrt(1 / (L * C))
    alpha = np.sqrt(gamma**2 - omega_0**2)
    return np.exp(-gamma * t) * voltage_0 * (np.cosh(alpha * t) + gamma / alpha * np.sinh(alpha * t))

def aperiodic_oscillations_amperage(t, C, L, R, voltage_0 = 1):
    gamma = R / (2.0 * L)
    omega_0 = np.sqrt(1 / (L * C))
    alpha = np.sqrt(gamma**2 - omega_0**2)
    R_crittical = 2.0 * np.sqrt(L / C)
    return -voltage_0 / R_crittical / 2.0 * omega_0 / alpha * np.exp(-gamma * t) * np.sinh(alpha * t)

def critical_oscillations_voltage(t, C, L, voltage_0 = 1):
    omega_0 = np.sqrt(1 / (L * C))
    return voltage_0 * (1 + omega_0 * t) * np.exp(-omega_0 * t)


def critical_oscillations_amperage(t, C, L, voltage_0 = 1):
    R = 2.0 * np.sqrt(L / C)
    omega_0 = np.sqrt(1.0 / (L * C))
    return -voltage_0 / R / 2.0 * omega_0 * t * np.exp(-omega_0 * t)

def attenuation_oscillations_voltage(t, C, L, R, voltage_0 = 1):
    gamma = R / (2.0 * L)
    omega_0 = np.sqrt(1.0 / (L * C))
    omega_1 = np.sqrt(omega_0**2 - gamma**2)
    U_0 = voltage_0 * omega_0 / omega_1
    fi_0 = -np.arctan(gamma / omega_1)
    return U_0 * np.exp(-gamma * t) * np.cos(omega_1 * t + fi_0)

def attenuation_oscillations_amperage(t, C, L, R, voltage_0 = 1):
    gamma = R / (2.0 * L)
    omega_0 = np.sqrt(1.0 / (L * C))
    omega_1 = np.sqrt(omega_0**2 - gamma**2)
    U_0 = voltage_0 * omega_0 / omega_1
    R_critical = 2.0 * np.sqrt(L / C)
    I_0 = U_0 / R_critical / 2
    return I_0 * np.exp(-gamma * t) * np.sin(omega_1 * t)

def main():
    # case 1: aperiodic oscillations

    C = 0.05 # mkFr
    L = 2000.0 # Gn
    R = 800.0 # Om

    t = np.linspace(0, 200, 200)
    voltage = aperiodic_oscillations_voltage(t, C, L, R)

    R_critical = 2.0 * np.sqrt(L / C)
    amperage = R_critical * aperiodic_oscillations_amperage(t, C, L, R)

    plots(t, [voltage, amperage], "t", ["votage", "amperage * R_critical"], "aperiodic_t")
    plots(voltage, [amperage], "voltage", ["amperage * R_kritical"], "aperiodic_u")

    #case 2: critical
    
    C = 0.05 # mkFr
    L = 2000.0 # Gn

    R_critical = 2.0 * np.sqrt(L / C)
    t = np.linspace(0, 200, 200)
    voltage = critical_oscillations_voltage(t, C, L)
    amperage = R_critical * critical_oscillations_amperage(t, C, L)

    plots(t, [voltage, amperage], "t", ["votage", "amperage * R_critical"], "critical_t")
    plots(voltage, [amperage], "voltage", ["amperage * R_kritical"], "critical_u")

    #case 3 : attenuation

    C = 0.2 # mkFr
    L = 200.0 # Gn
    R = 10.0 # Om

    R_critical = 2.0 * np.sqrt(L / C)
    t = np.linspace(0, 200, 200)
    voltage = attenuation_oscillations_voltage(t, C, L, R)
    amperage = R_critical * attenuation_oscillations_amperage(t, C, L, R)

    plots(t, [voltage, amperage], "t", ["votage", "amperage * R_critical"], "attenuation_t")
    plots(voltage, [amperage], "voltage", ["amperage * R_kritical"], "attenuation_u")

    return

main()    
