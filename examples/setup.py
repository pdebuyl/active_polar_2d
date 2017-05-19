from setuptools import setup, Extension
from Cython.Build import cythonize

import active_polar_2d
import os.path

include = os.path.dirname(active_polar_2d.__file__)

setup(
    ext_modules=cythonize(Extension('wall_force_torque', ["wall_force_torque.pyx"], include_dirs=[include]))
)
