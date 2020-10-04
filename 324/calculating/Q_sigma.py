import numpy as np

def sigma_coss(val, sigsArr, parmsArr, powerArr):
    sum = 0
    for i in range(len(sigsArr)):
        sum += (powerArr[i] * sigsArr[i] / parmsArr[i])**2
    return val * (sum)**0.5

def sigma_log(val, sigma_val):
    return sigma_val / val

def main():
    x = np.array([1.0000, 0.6944, 0.4444, 0.3086, 0.2500, 0.1600, 0.1111])
    y = np.array([3.7351, 2.6166, 1.8794, 1.2272, 1.0564, 0.4294, 0.1539])

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

    print sigma_coss(3.2, [sigma_teta[4]], [teta[4]], [1])

main()