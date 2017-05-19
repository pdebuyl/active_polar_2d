from active_polar_2d.core cimport cyfunc_2d
from libc.math cimport cos, sin

cdef class wall(cyfunc_2d):
    cdef double r_cut
    cdef double L
    cdef double bond_l
    cdef double grav

    cdef double wall_f(self, double x):
        cdef int sign
        cdef double r
        if x < self.r_cut:
            r = x
            sign = 1
        elif x > self.L-self.r_cut:
            r = self.L-x
            sign = -1
        else:
            return 0
        return sign*3*(3*r**-10 - r**-4)

    cdef void force(self, double x, double y, double theta, double *fx, double *fy, double *torque):
        cdef double c_th, s_th, x1, x2, f1, f2
        c_th = cos(theta)
        s_th = sin(theta)
        x1 = x + c_th*self.bond_l
        f1 = self.wall_f(x1)
        x2 = x - c_th*self.bond_l
        f2 = self.wall_f(x2)
        fx[0] = f1 + f2 - self.grav
        fy[0] = 0
        torque[0] = s_th*self.bond_l/2*(f2-f1)

    def __init__(self, L, bond_l, grav):
        self.r_cut = 3.**(1./6.)
        self.L = L
        self.bond_l = bond_l
        self.grav = grav

