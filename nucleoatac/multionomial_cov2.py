from numba import jit

@jit("f8(f8[:], f8[:], f8)", nopython=True)
def calculateCov(p, v, r):
    if p.shape[0]!= v.shape[0]:
        raise ValueError("p and v must be same shape")
    value = 0.0
    for i in range(p.shape[0]):
        for j in range(i,p.shape[0]):
            if i==j:
                value += p[i]* (1-p[i]) * v[i]**2
            else:
                value += p[i] * p[j] * -2 * v[i] * v[j]
    return(value * r)