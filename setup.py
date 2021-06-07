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
    },
    install_requires=['numpy', 'deap']
)

setup (
    name='profilerGen',
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
    ],
    entry_points={
        'console_scripts': [
            'profilerGen=profilerGen.profilerGen:main',
        ]
    },
    install_requires=['numpy', 'deap']
)
