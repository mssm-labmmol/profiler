from setuptools import setup, Extension, find_packages
import sys

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

opt_ext = []
gen_ext = []

if '--without-geometrypy' in sys.argv:
    sys.argv.remove('--without-geometrypy')
else:
    opt_ext.append(extended_numpy)
    gen_ext.append(extended_numpy)

if '--without-gsl' in sys.argv:
    sys.argv.remove('--without-gsl')
else:
    opt_ext.append(ext_gsl_ridge)


setup (
    name='profilerOpt',
    version=package_version,
    python_requires='>=3.6',
    packages=find_packages(),
    ext_modules=opt_ext,
    entry_points={
        'console_scripts': [
            'profilerOpt=profilerTools.profilerOpt.profilerOpt:main',
        ]
    },
    install_requires=['numpy', 'deap', 'scikit-learn'],
)

setup (
    name='profilerGen',
    version=package_version,
    python_requires='>=3.6',
    packages=find_packages(),
    ext_modules=gen_ext,
    entry_points={
        'console_scripts': [
            'profilerGen=profilerTools.profilerGen.profilerGen:main',
        ]
    },
    install_requires=['numpy', 'deap', 'scikit-learn'],
)
