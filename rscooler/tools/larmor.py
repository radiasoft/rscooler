import numpy as np
from scipy import optimize


def findgyrostuff(array):
    """
    Find the Larmor radius and gyrocenter for a particle executing many revolutions.
    Assumes the circle will be well filled.
    Parameters
    ----------
    array: Array particle positions in (x,xp,y,yp,z,zp)

    Returns
    -------
    (grho, gcen): Tuple with Larmor radius and gyrocenter

    """
    xmin = np.argmin(array[:, 0])
    xmax = np.argmax(array[:, 0])

    ymin = np.argmin(array[:, 2])
    ymax = np.argmax(array[:, 2])

    grho = (abs(array[ymax, 2] - array[ymin, 2]) + abs(array[xmax, 0] - array[xmin, 0])) / 4.
    gcen = ((array[ymax, 0] + array[ymin, 0]) / 2., (array[xmax, 2] + array[xmin, 2]) / 2.)

    return grho, gcen


def _findR(x, y, xc, yc):
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def _f(c, x, y):
    R = _findR(x, y, *c)
    return R - np.average(R)

def opt_circle(x, y):
    """
    Preferred method for calculating Larmor radius and gyrocenter.
     Least squares fit to coordinates of a circle.
     Based on example from Scipy cookbook
    Parameters
    ----------
    x: Particle x coordinates
    y: Particle y coordinates

    Returns
    -------
    xc: Coordinate in x of circle center
    yc: Coordinate in y of circle center
    R: Circle radius
    residual: Residual of the fit (sum(Ri - R)**2)

    """
    x_m = np.average(x)
    y_m = np.average(y)
    ce = x_m, y_m
    center, ier = optimize.leastsq(_f, ce, args=(x, y))
    xc, yc = center
    Ri = _findR(x, y, *center)
    R = np.average(Ri)
    residual = np.sum((Ri - R)**2)
    return xc, yc, R, residual

def nrGR(B, v):
    """
    Calculate the nonrelativistic Larmor radius.
    Parameters
    ----------
    B: Magnetic field in Tesla
    v: Transverse velocity in m/s

    Returns
    -------
    r: Larmor radius in m

    """
    me = 9.10938356e-31

    p = me * v
    e = 1.6e-19
    r = p / (e * B)
    return r

def relGR(B, v, vl):
    """
    Calculate the relativistically corrected Larmor radius
    Parameters
    ----------
    B: Magnetic field in Tesla
    v: Transverse velocity in m/s
    vl: Total velocity in m/s

    Returns
    -------

    """
    c = 299792458
    me = 9.10938356e-31

    vtot = np.sqrt(v**2 + vl**2)
    gamma = 1/np.sqrt(1-vtot**2/c**2)
    p = gamma * me * v
    r = p / (1.6e-19 * B)
    return r