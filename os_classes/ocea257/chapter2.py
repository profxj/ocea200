import numpy as np


'''
Problem 2.2
'''

def euler_inertial(n, u0, f, Dt):
    """

    Parameters
    ----------
    n : int
        Number of time steps
    u0
    f
    Dt

    Returns
    -------
    x, y, u, v : np.ndarray

    """
    # Initial conditions
    u = [u0]
    v = [0.]
    x = [0.]
    y = [0.]

    # Proceed
    for ii in range(n):
        # Velocity
        u += [u[-1] + f*Dt*v[-1]]
        v += [v[-1] - f*Dt*u[-2]] # Note we had already updated u!
        # Position
        x += [x[-1] + Dt*u[-1]]
        y += [y[-1] + Dt*v[-1]]
    # Return
    return np.array(x), np.array(y), np.array(u), np.array(v)

def semi_implicit_inertial(n, u0, f, Dt):
    """

    Parameters
    ----------
    n : int
        Number of time steps
    u0
    f
    Dt

    Returns
    -------
    u, v : np.ndarray

    """
    # Initial conditions
    u = [u0]
    v = [0.]

    # Proceed
    for ii in range(n):
        # Velocity
        u += [ (f*Dt*v[-1] + u[-1]*(1-0.25*(f*Dt*f*Dt)))/(1+0.25*f*Dt*f*Dt) ]
        v += [(-f * Dt * u[-1] + v[-1] * (1 + 0.25 * (f * Dt * f * Dt))) / (1 - 0.25 * f * Dt * f * Dt)]
    # Return
    return np.array(u), np.array(v)
