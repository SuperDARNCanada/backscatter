import backscatter.core.core as core
from .dmap import dmap
from .fitacf import fitacf

# parsing config files
config = core.parse_config_file()
hdw_files_path = config.get("core", "hdw_files_path")
hdw_info = core.parse_hdw_files(hdw_files_path)

# adding modules that will be seen at import
__all__ = ["fitacf", "dmap"]








