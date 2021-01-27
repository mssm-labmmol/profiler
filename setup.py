from setuptools import setup, Extension, find_packages

import numpy

extended_numpy = Extension(
    'geometrypy',
    sources = [
        'geometrypy/geometrypy.c',
        ]
)

setup (
    name='ProfilerOpt',
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
    ],
    entry_points={
        'console_scripts': [
            'ProfilerOpt=profilerOpt.profilerOpt:main',
        ]
    }
)

setup (
    name='ProfilereGen',
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
    ],
    entry_points={
        'console_scripts': [
            'ProfilerGen=profilerGen.profilerGen:main',
        ]
    }
)
