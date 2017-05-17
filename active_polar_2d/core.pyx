from libc.math cimport sqrt, cos, sin
import numpy as np
cimport numpy as np
from threefry cimport rng

cdef class cyfunc_2d:
    cdef void force(self, double x, double y, double theta, double *fx, double *fy, double *torque):
        fx[0] = 0
        fy[0] = 0
        torque[0] = 0
    def __init__(self):
        pass

cdef class pyfunc_2d(cyfunc_2d):
    cdef void force(self, double x, double y, double theta, double *fx, double *fy, double *torque):
        cdef double tmp_fx, tmp_fy, tmp_torque
        tmp_fx, tmp_fy, tmp_torque = self.py_force(x, y, theta)
        fx[0] = tmp_fx
        fy[0] = tmp_fy
        torque[0] = tmp_torque
    def __init__(self, force):
        self.py_force = force

def integrate_OD_2d_theta(double x, double y, double th, double gamma_par, double gamma_perp, double gamma_rot, double T, double v0, double dt, int interval, int steps, f=None, seed=None):
    cdef cyfunc_2d cy_f
    cdef int t_idx1, t_idx2, j, n_dims
    cdef double noise, th_noise

    r = rng(seed)

    if f is None:
        cy_f = cyfunc_2d()
    elif isinstance(f, cyfunc_2d):
        cy_f = f
    elif callable(f):
        cy_f = pyfunc_2d(f)
    else:
        raise ValueError("f should be a callable or a cyfunc_2d")


    cdef double D_par = T/gamma_par
    cdef double mu_par = 1/gamma_par
    cdef double noise_par = sqrt(2*D_par/dt)
    cdef double D_perp = T/gamma_perp
    cdef double mu_perp = 1/gamma_perp
    cdef double noise_perp = sqrt(2*D_perp*dt)
    cdef double D_rot = T/gamma_rot
    cdef double mu_rot = 1/gamma_rot
    cdef double noise_rot = sqrt(2*D_rot/dt)

    cdef double fx, fy, torque
    cdef double f_par, f_perp
    cdef double c_th, s_th
    cdef double xi_par, xi_perp

    cdef double[::1] x_out = np.empty((steps,), dtype=float)
    cdef double[::1] y_out = np.empty((steps,), dtype=float)
    cdef double[::1] th_out = np.empty((steps,), dtype=float)

    for t_idx1 in range(steps):
        for t_idx2 in range(interval):
            cy_f.force(x, y, th, &fx, &fy, &torque)
            c_th = cos(th)
            s_th = sin(th)
            f_par = c_th*fx + s_th*fy
            f_perp = -s_th*fx + c_th*fy
            xi_par = r.random_normal()
            xi_perp = r.random_normal()

            x = x + dt * (c_th * (mu_par*f_par + noise_par*xi_par + v0) - s_th * (mu_perp*f_perp+noise_perp*xi_perp) )
            y = y + dt * (s_th * (mu_par*f_par + noise_par*xi_par + v0) + c_th * (mu_perp*f_perp+noise_perp*xi_perp) )
            th = th + dt * (mu_rot*torque + noise_rot*r.random_normal())
            

        x_out[t_idx1] = x
        y_out[t_idx1] = y
        th_out[t_idx1] = th

    return np.asarray(x_out), np.asarray(y_out), np.asarray(th_out)
