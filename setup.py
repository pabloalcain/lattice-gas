#!/usr/bin/env python
"""
Script for installation of the python module.

It also builds the c library
"""
from distutils.core import setup, Extension
#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.

files = ["latgas/*"]

lat_module = Extension('liblatgas', 
                       sources = ['src/evolution.c'],
                       libraries = ['m'],
                       extra_compile_args = ['-O3', '-g'],
                      )

setup(name = "Lattice Gas in 2D/3D",
      version = "0.1",
      description = "Lattice Gas implementation in Python/c",
      author = "Pablo Alcain",
      author_email = "pabloalcain@gmail.com",
      url = "none yet",
      packages = ['latgas'],
      package_data = {'package' : files},
      scripts = ["examples/test.py", "examples/lengths.py"],
      long_description = 
      """
Python implementation of a Lattice Gas.

The idea is to be able to set interactions as desired and work
*directly* from a Python interpreter!

So far everything is implemented in Python; obviously, soon enough,
calculation-heavy routines will be ported to c/c++
""",
      ext_modules = [lat_module]) 
