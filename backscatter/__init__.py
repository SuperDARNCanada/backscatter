import ConfigParser as cp # REVIEW #7 License at the top of the file?
import os
import core
# REVIEW #40 Need two lines here
#parsing config files # REVIEW #1 space after '#' required, same goes for below comment
config = core.parse_config_file()

hdw_files_path = config.get("core","hdw_files_path") # REVIEW #40 missing whitespace after comma

hdw_info = core.parse_hdw_files(hdw_files_path)

#adding modules that will be seen at import                
__all__ = ["fitacf","dmap"] # REVIEW #40 missing whitespace after comma

from . import dmap
from . import fitacf # REVIEW #40 pycharm complains about  module level import not at top of file






# REVIEW # 40 too many blank lines