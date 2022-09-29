"""
This file is for when processing a SuperDARN RAWACF file to FITACF from the command line.
The syntax for this is 'python -m backscatter.fitacf rawacf_file fitacf_file'
"""
from .fitacf import main
main()
