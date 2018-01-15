import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "philsol",
    version = "0.20",
    author = "Phil Main, Joel Collins",
    author_email = "philmain28@gmail.com",
	install_requires=['numpy', 'scipy'],
    description = ("A fully vectorial finite difference waveguide mode solver. Based on the algorithm of Zhu and Brown"),
    license = "GNU GPL",
    keywords = "vector finite difference waveguide mode solver",
    url = "https://github.com/philmain28/philsol",
    packages=['philsol'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)