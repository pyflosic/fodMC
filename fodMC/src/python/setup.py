import setuptools
from numpy.distutils.core import Extension, setup

setup(
   name="pyfodmc",
   version="0.1",
   author="Sebastian Schwalbe",
   author_email="theonov13@gmail.com",
   description="Example package to demonstrate wheel issue",
   packages = setuptools.find_packages(),
   ext_modules=[Extension(name='fodmc', sources=['pyfodmc/lib/fodmc.f90'], f2py_options=['--quiet'])]
)
