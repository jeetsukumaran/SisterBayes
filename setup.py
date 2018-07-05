#! /usr/bin/env python

from setuptools import setup
import os
import sys
import re

def _read(names, **kwargs):
    path = os.path.join(os.path.dirname(__file__), *names)
    if sys.version_info.major < 3:
        return open(path, "rU").read()
    else:
        with open(path, encoding=kwargs.get('encoding', 'utf8')) as src:
            s = src.read()
        return s

__project__ = "SisterBayes"
__version__ = re.match(r".*^__version__\s*=\s*['\"](.*?)['\"]\s*$.*", _read(["src", "sisterbayes", "__init__.py"]), re.S | re.M).group(1)
sys.stderr.write("-setup.py: {} {}\n".format(__project__, __version__))

setup(
    name="sisterbayes",
    version=__version__,
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=["sisterbayes"],
    package_dir={"":"src"},
    scripts=[
        "bin/sisterbayes-simulate.py",
        "bin/sisterbayes-normalize.py",
        "bin/sisterbayes-sumstats.py",
        "bin/sisterbayes-reject.py",
        "bin/sisterbayes-summarize.py",
        ],
    url="http://pypi.python.org/pypi/sisterbayes/",
    test_suite = "tests",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.txt").read(),
    zip_safe = True,
    # install_requires=[ ],
)
