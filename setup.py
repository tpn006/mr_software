import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 1
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'mrsoftware/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )


setup(
    name='mrsoftware',
    version=VERSION,
    description='CSE185 Demo Project',
    author='Jonathan Narita and Tasha Nguyen',
    author_email='tpn006@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "mrsoftware=mrsoftware.mrsoftware:main"
        ],
    }
)