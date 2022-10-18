from setuptools import setup
from Cython.Build import cythonize

setup(
    name='make_pair',
    ext_modules=cythonize("make_pair.pyx"),
    zip_safe=False,
)
