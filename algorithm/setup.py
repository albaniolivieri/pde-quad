from distutils.core import setup
from Cython.Build import cythonize
setup(
    ext_modules = cythonize(["branch_and_bound.pyx",
                             "nearest_neighbor.pyx",
                             "quadratization.pyx",
                             "PolySys.pyx"])
)