# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

from setuptools import setup, find_packages
from setuptools.command.install import install

import errno
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


class CustomInstall(install):
    """Customized setuptools install command - creates configuration files."""
    def run(self):
        self.create_config_files()
        install.run(self)

    @staticmethod
    def create_config_files():
        # Creating default configuration files

        import configparser as cp
        import os

        config = cp.RawConfigParser()

        config.add_section('fitacf')
        config.set('fitacf', 'w_max', '90.0')
        config.set('fitacf', 'v_max', '30.0')
        config.set('fitacf', 'fitacf_revision_minor', '0')
        config.set('fitacf', 'fitacf_revision_major', '3')
        config.set('fitacf', 'fluctuation_cutoff_coeff', '2')
        config.set('fitacf', 'alpha_cutoff', '2.0')
        config.set('fitacf', 'acf_snr_cutoff', '1.0')
        config.set('fitacf', 'minimum_lags', '3')

        config.add_section('core')
        config.set('core', 'hdw_files_path', '/usr/local/hdw')

        for loc in os.path.expanduser("~"), "/etc/backscatter":
            file_path = os.path.join(loc, "backscatter.ini")

            if not os.path.exists(os.path.dirname(file_path)):
                try:
                    os.makedirs(loc)
                except OSError as exc:      # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        if exc.errno == errno.EACCES:
                            err_msg = """Could not create installation folder
                            at {0}. Please manually create or reinstall as
                            root.""".format(loc)
                            print(err_msg)
                        else:
                            raise

            try:
                with open(file_path, 'w') as cfg:
                    config.write(cfg)
            except OSError:
                err_msg = """Could not install configuration file
                at {0}. Please manually copy or reinstall as root.""".format(loc)
                print(err_msg)


# Set up arguments for installation
setup_args = {
    'name': "backscatter",
    'version': "2016.08",
    'description': "A Python package of analysis tools for SuperDARN data",
    'url': "",
    'author': "SuperDARN Canada",
    'license': "GNU",
    'packages': find_packages(exclude=['contrib', 'docs', 'tests']),
    'setup_requires': ['configparser'],
    'install_requires': ['numpy>=1.8'],
    'cmdclass': {'install': CustomInstall},
}

setup(**setup_args)
