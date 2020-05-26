import setuptools
from numpy.distutils.core import Extension, setup

with open("README.md", "r") as fh:
     name="fodMC",                                                                          
     long_description = fh.read()  

setup(
   name="fodMC",
   version="1.0.1",
   author="Kai Trepte",
   author_email="kai.trepte1987@gmail.com",
   description="Fermi-orbital descriptor Monte-Carlo",
   long_description=long_description,
   long_description_content_type="text/markdown",
   include_package_data=True,
   packages = setuptools.find_packages(),
   ext_modules=[Extension(name='fodmc', sources=['fodMC/lib/fodmc.f90'], f2py_options=['--quiet'])]
)
