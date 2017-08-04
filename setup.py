#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(
    name="sisterbayes",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=["sisterbayes"],
    scripts=[
        "bin/sisterbayes-simulate.py",
        "bin/sisterbayes-normalize.py",
        "bin/sisterbayes-sumstats.py",
        "bin/sisterbayes-reject.py",
        "bin/sisterbayes-summarize.py",
        ],
    url="http://pypi.python.org/pypi/sisterbayes/",
    test_suite = "sisterbayes.test",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.txt").read(),
    # install_requires=[ ],
)
