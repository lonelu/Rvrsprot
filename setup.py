#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='Rvrsprot',
    version='0.0.0',
    author='Lei Lu',
    author_email='lonelur@gmail.com',
    url='https://github.com/lonelu//Rvrsprot',
    description='protein design',
    license='MIT',
    packages=find_packages(exclude=['tests*']),
    install_requires=[
        'numpy',
        'matplotlib',
        'prody',
        'sklearn',
        'numba',
        'scipy',
        #'dataclasses',
    ],
    # extras_require = {
    #     '':  [],
    # },
    entry_points={
        'console_scripts': [ ],
    },
    
    #long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
    ],
)
