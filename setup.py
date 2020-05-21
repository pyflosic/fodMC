import setuptools
from numpy.distutils.core import Extension, setup

setup(
   name="fodMC",
   version="1.0",
   author="Kai Trepte",
   author_email="kai.trepte1987@gmail.com",
   description="Fermi-orbital descriptor generator",
   packages = setuptools.find_packages(),
   ext_modules=[Extension(name='fodmc', sources=['fodMC/lib/fodmc.f90'], f2py_options=['--quiet'])]
)
