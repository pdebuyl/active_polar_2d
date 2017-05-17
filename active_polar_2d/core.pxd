

cdef class cyfunc_2d:
    cdef void force(self, double x, double y, double theta, double *fx, double *fy, double *torque)


cdef class pyfunc_2d(cyfunc_2d):
    cdef object py_force

