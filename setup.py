#! /usr/bin/env python


#from distutils.core import setup
from setuptools import setup
import sys
from sisterbayes import __version__, __project__
sys.stderr.write("-setup.py: {} {}\n".format(__project__, __version__))


setup(
    name="sisterbayes",
    version=__version__,
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
