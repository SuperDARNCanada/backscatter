import os
import configparser as cp

def parse_hdw_files(hdw_files_path):

    hdw_info = {}
    
    loop_ran = False
    try:
        for filename in os.listdir(hdw_files_path):
            loop_ran = True

            # Most code to read new hdw format taken from code used in Borealis experimentoptions.py
            if "hdw.dat" in filename:
                with open(os.path.join(hdw_files_path, filename), 'r') as f:
                    lines = f.readlines()
                    # Remove comments, blank lines
                    lines[:] = [line for line in lines if line[0] != "#"]
                    lines[:] = [line for line in lines if len(line.split()) != 0]

                    try:
                        hdw = lines[-1]
                    except IndexError:
                        errmsg = 'Cannot find any valid lines in the hardware file: ' \
                                 '{}'.format(filename)
                        raise IndexError(errmsg)

                    # we now have the correct line of data.
                    params = hdw.split()
                    if len(params) != 22:
                        errmsg = 'Found {} parameters in hardware file, expected 22'.format(
                            len(params))
                        raise IndexError(errmsg)

                    hdw_dict = {}

                    hdw_dict["stid"] = int(params[0])
                    hdw_dict["status"] = int(params[1])  # 1 operational, -1 offline

                    # Formats: YYYYMMDD, and HH:MM:SS
                    hdw_dict["date_valid"] = params[2]
                    hdw_dict["time_valid"] = params[3]

                    hdw_dict["lat"] = float(params[4])  # decimal degrees, S = negative
                    hdw_dict["lon"] = float(params[5])  # decimal degrees, W = negative
                    hdw_dict["altitude"] = float(params[6])  # metres
                    hdw_dict["boresight"] = float(params[7])  # degrees from geographic north, CCW = neg
                    hdw_dict["boresight_shift"] = float(params[8])  # degrees from physical boresight. nominal 0.0 degrees
                    hdw_dict["beamsep"] = float(params[9])  # degrees, nominal 3.24 degrees
                    hdw_dict["velsign"] = float(params[10])  # +1.0 or -1.0
                    hdw_dict["phasesign"] = float(params[11])  # +1 indicates correct int phase, -1 otherwise
                    hdw_dict["tdiff_a"] = float(params[12])  # ns for channel A
                    hdw_dict["tdiff_b"] = float(params[13])  # ns for channel B

                    # Meters, int offset from midpoint of main. [x, y, z] where x is along line of
                    # antennas, y is along array normal, z is altitude difference
                    hdw_dict["interoffx"] = float(params[14])
                    hdw_dict["interoffy"] = float(params[15])
                    hdw_dict["interoffz"] = float(params[16])

                    hdw_dict["rxrisetime"] = float(params[17])  # us
                    hdw_dict["rxatten"] = float(params[18])  # dB
                    hdw_dict["attenstages"] = float(params[19])  # number of stages

                    hdw_dict["maxgates"] = int(params[20])
                    hdw_dict["maxbeams"] = int(params[21])  # so beamnum points in a certain dir
                    break
    except OSError as e:
        error_msg = "No directory {0} found when locating hdw files!".format(hdw_files_path)      
        raise ImportError(error_msg)


    if loop_ran == False or len(hdw_info) == 0:
        error_msg = "No hardware files found at {0}".format(hdw_files_path)
        raise ImportError(error_msg)

    return hdw_info


def parse_config_file():
    config = cp.ConfigParser()

    valid_sections = ["core","fitacf"]
    for loc in os.curdir, os.path.expanduser("~"), "/etc/backscatter":
        
        file_path = os.path.join(loc,"backscatter.ini")
        config.read(file_path)

        if config.sections():
            if set(config.sections()).issubset(set(valid_sections)):
                return config
            else:
                raise ValueError("Config file is missing required sections")
    else:
        raise IOError("Config file missing!")
