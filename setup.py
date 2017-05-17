from setuptools import setup, Extension
from Cython.Build import cythonize

try:
    import threefry
    threefry_include = threefry.get_include()
except:
    threefry_include = ''

active_polar_2d_ext = cythonize(Extension("active_polar_2d.core",
                     sources=["active_polar_2d/core.pyx"],
                     include_dirs=["active_polar_2d", threefry_include],
                     ))

setup(name='active_polar_2d',
      version='0.1.0.dev0',
      description='Numerical simulation of the Langevin equation for 2D active particles',
      author='Pierre de Buyl',
      license='BSD',
      packages=['active_polar_2d'],
      ext_modules = active_polar_2d_ext,
      setup_requires=['cython', 'threefry'],
      package_data={'active_polar_2d': ['core.pxd']},
      )
