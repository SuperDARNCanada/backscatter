# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


try:
    from setuptools import setup
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install


class CustomInstall(install):
    """Customized setuptools install command - creates configuration files."""
    def run(self):
        self.create_config_files()
        install.run(self)
    
    def create_config_files(self):
        #Creating default configuation files

        import ConfigParser as cp
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
        config.set('core', 'hdw_files_path', '/etc/backscatter')
        
        for loc in os.curdir, "/etc/backscatter":
        
            file_path = os.path.join(loc,"backscatter.ini")
        
            if not os.path.exists(os.path.dirname(file_name)):
                try:
                    os.makedirs(loc)
                except OSError as exc: # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            
                with open(file_path,'wb') as cfg:
                    config.write(cfg)






#Set up arguments for installation
setup_args  =  {
    'name' : "backscatter",
    'version' : "2016.08",
    'description' : "A Python package of analysis tools for SuperDARN data",
    'url' : "",
    'author' : "SuperDARN Canada",
    'license' : "GNU",
    'packages' : ["backscatter","tests"],
    'cmdclass' : {'install' : CustomInstall,},
}

setup(**setup_args)




    