from setuptools import setup, Extension, find_packages

import numpy

extended_numpy = Extension(
    'geometrypy',
    sources = [
        'geometrypy/geometrypy.c',
        ]
)

setup (
    name='profilerOpt',
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
    ],
    entry_points={
        'console_scripts': [
            'profilerOpt=profilerOpt.profilerOpt:main',
        ]
    }
)

setup (
    name='profilereGen',
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
    ],
    entry_points={
        'console_scripts': [
            'profilerGen=profilerGen.profilerGen:main',
        ]
    }
)
