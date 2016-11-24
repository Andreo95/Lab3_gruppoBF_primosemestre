import numpy as np, scipy.special as sps, scipy.stats as sts


def Phi(x):                     #CDF of normal distribution
    return (sps.erf(x / np.sqrt(2)) + 1) / 2


def normProb(x, y):             # probability of finding a normally distributed variable between x and y
    if max(y, x) > 0: return sts.norm.cdf(-x) - sts.norm.cdf(-y)
    return sts.norm.cdf(y) - sts.norm.cdf(x)


def normalizeVar(x, m, s):      # transform normally distributed variable x (mean m, stdev s) into standard form (mean 0, stdev 1)
    return (x - m)/s


def intervalProb(x, m, s):      # given array x, outputs the array of the probabilities that a normal variable (mean m, stdev s) is found in the intervals delimited by x's values
    localZ = (x - m)/s
    localOut = np.zeros(len(x)+1)
    localOut[0] = normProb(-np.inf, localZ[0])
    for i in range(1, len(x)):
        localOut[i] = normProb(localZ[i-1], localZ[i])
    localOut[-1] = normProb(localZ[-1], np.inf)
    return localOut


def chiSquare(x, y, s):         # Returns the chi-square of data x (with error s)
    return np.sum( ((x-y)/s)**2 )


def chiAvg(x, s=1):              # error-weighted average of data x
    if (type(s) != np.ndarray):
        s = np.ones(len(x)) * s
    avg = np.sum(x/s**2) / np.sum(1/s**2)
    sigma = 1 / np.sqrt( np.sum(1/s**2) )
    return np.array([avg, sigma])
    

def chiLineFit(x, y, s=1):      # finds m and q for y = mx + q, with errors s on y
    if (type(s) != np.ndarray):
        s = np.ones(len(x)) * s
    m = ( sum(1/s**2)*sum(x*y/s**2) - sum(x/s**2)*sum(y/s**2) ) / ( sum(1/s**2)*sum(x**2/s**2) - sum(x/s**2)**2 )
    q = ( sum(x**2/s**2)*sum(y/s**2) - sum(x/s**2)*sum(x*y/s**2) ) / ( sum(1/s**2)*sum(x**2/s**2) - sum(x/s**2)**2 )
    sigma_q = np.sqrt( sum((x/s)**2) / ( sum(1/s**2)*sum(x**2/s**2) - sum(x/s**2)**2 ) )
    sigma_m = np.sqrt( sum(1/s**2) / ( sum(1/s**2)*sum(x**2/s**2) - sum(x/s**2)**2 ) )
    return np.array([[m, q], [sigma_m, sigma_q]])


def chiRatio(x, y, s=1):         # line fit with q=0
    if (type(s) != np.ndarray):
        s = np.ones(len(x)) * s
    m = np.sum((x/s)**2) / np.sum(x*y/s**2)
    return m