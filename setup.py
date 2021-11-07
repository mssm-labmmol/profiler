from setuptools import setup, Extension, find_packages

with open("VERSION") as fp:
    package_version = fp.read().strip()

extended_numpy = Extension(
    'geometrypy',
    sources = [
        'profilerTools/geometrypy/geometrypy.c',
        ]
)

ext_gsl_ridge = Extension(
    'gslpyridge',
    sources = [
        'profilerTools/gslpy_ridge/gslpy_ridge.c',
    ],
    libraries = [
        'gsl', 'gslcblas'
    ],
)

setup (
    name='profilerOpt',
    version=package_version,
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
        ext_gsl_ridge,
    ],
    entry_points={
        'console_scripts': [
            'profilerOpt=profilerTools.profilerOpt.profilerOpt:main',
        ]
    },
    install_requires=['numpy', 'deap', 'scikit-learn']
)

setup (
    name='profilerGen',
    version=package_version,
    packages=find_packages(),
    ext_modules=[
        extended_numpy,
    ],
    entry_points={
        'console_scripts': [
            'profilerGen=profilerTools.profilerGen.profilerGen:main',
        ]
    },
    install_requires=['numpy', 'deap', 'scikit-learn']
)
