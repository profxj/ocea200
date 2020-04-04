import numpy as np


'''
Problem 1.2
'''

def expo_func(omega, t):
    return np.exp(-1*omega*(t + np.abs(t)))

def first_forward(t, wdt, omega=1.):
    # For convenience, set omega = 1
    dt = wdt / omega
    # Evaluate
    tdt = t + dt
    expf = expo_func(omega, t)
    expf_plus = expo_func(omega, tdt)
    # Finite difference, forward
    dexpdt_forward = (expf_plus-expf)/dt
    return dexpdt_forward

def first_backward(t, wdt, omega=1.):
    # For convenience, set omega = 1
    dt = wdt / omega
    # Evaluate
    tdt = t - dt
    expf = expo_func(omega, t)
    expf_neg = expo_func(omega, tdt)
    # Finite difference, forward
    dexpdt_back = (expf-expf_neg)/dt
    return dexpdt_back

def second_centered(t, wdt, omega=1.):
    # For convenience, set omega = 1
    dt = wdt / omega
    # Evaluate
    tdtp = t + dt
    tdtn = t - dt
    expf_pos = expo_func(omega, tdtp)
    expf_neg = expo_func(omega, tdtn)
    # Finite difference, forward
    dexpdt_second = (expf_pos-expf_neg)/(2*dt)
    return dexpdt_second

'''
Problem 1.2
'''

def sinh_func(k, x):
    return np.sinh(k*x)

def sinh_first_forward(x, kdx, k=1.):
    dx = kdx / k
    # Evaluate
    xdx = x + dx
    expf = sinh_func(k, x)
    expf_plus = sinh_func(k, xdx)
    # Finite difference, forward
    dexpdt_forward = (expf_plus-expf)/dx
    return dexpdt_forward

def sinh_first_backward(x, kdx,k=1.):
    dx = kdx / k
    # Evaluate
    xdx = x - dx
    expf = sinh_func(k, x)
    expf_neg = sinh_func(k, xdx)
    # Finite difference, forward
    dexpdt_forward = (expf-expf_neg)/dx
    return dexpdt_forward

def sinh_second_center(x, kdx,k=1.):
    dx = kdx / k
    # Evaluate
    xdxn = x - dx
    xdxp = x + dx
    sinh_pos = sinh_func(k, xdxp)
    sinh_neg = sinh_func(k, xdxn)
    # Finite difference, forward
    dfdx_second = (sinh_pos-sinh_neg)/(2*dx)
    return dfdx_second

def sinh_fourth(x, kdx,k=1.):
    dx = kdx / k
    # Evaluate
    xdxn = x - dx
    xdxp = x + dx
    xdxn2 = x - 2*dx
    xdxp2 = x + 2*dx
    sinh_pos = sinh_func(k, xdxp)
    sinh_neg = sinh_func(k, xdxn)
    sinh_pos2 = sinh_func(k, xdxp2)
    sinh_neg2 = sinh_func(k, xdxn2)
    # Finite difference, forward
    dfdx_fourth = (4./3) * (sinh_pos-sinh_neg)/2/dx - (1/3)*(sinh_pos2-sinh_neg2)/4/dx
    return dfdx_fourth

