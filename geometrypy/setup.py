from distutils.core import setup, Extension
import numpy

extended_numpy = Extension('geometrypy',
                    sources = ['geometrypy.c'])

setup (name = 'GeometryPy',
       version = '1.0',
       description = 'This is a library to facilitate the calculation of geometric degrees-of-freedom using Numpy.',
       ext_modules = [extended_numpy])
