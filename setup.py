# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup_args  :  {
    'name' : "backscatter",
    'version' : "2016.08",
    'description' : "A Python package of analysis tools for SuperDARN data",
    'long_description' : read('README.rst')
    'url' : "",
    'author' : "SuperDARN Canada",
    'license' : "GNU",
    'packages' : ["backscatter","tests"]
}

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(**setup_args)