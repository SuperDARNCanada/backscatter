import os
import configparser as cp


def parse_hdw_files(hdw_files_path):

    hdw_info = {}

    loop_ran = False
    try:
        for filename in os.listdir(hdw_files_path):
            loop_ran = True

            if "hdw.dat" in filename:
                with open(os.path.join(hdw_files_path,filename),'r') as f:
                    backwards_lines = reversed(f.readlines())

                    for line in backwards_lines:
                        if '#' in line or not line.strip():
                            continue
                        else:
                            fields = line.split()
                            if len(fields) != 19:
                                raise ValueError("Invalid number of hdw.dat fields!")

                            hdw_dict = {}

                            hdw_dict["stid"] = int(fields[0])
                            hdw_dict["year"] = int(fields[1])
                            hdw_dict["sc"] = int(fields[2])
                            hdw_dict["lat"] = float(fields[3])
                            hdw_dict["lon"] = float(fields[4])
                            hdw_dict["altitude"] = float(fields[5])
                            hdw_dict["boresight"] = float(fields[6])
                            hdw_dict["beamsep"] = float(fields[7])
                            hdw_dict["velsign"] = float(fields[8])
                            hdw_dict["rxattenstep"] = float(fields[9])
                            hdw_dict["tdiff"] = float(fields[10])
                            hdw_dict["phasesign"] = float(fields[11])
                            hdw_dict["interoffx"] = float(fields[12])
                            hdw_dict["interoffy"] = float(fields[13])
                            hdw_dict["interoffz"] = float(fields[14])
                            hdw_dict["rxrisetime"] = float(fields[15])
                            hdw_dict["attenstages"] = float(fields[16])
                            hdw_dict["maxgates"] = int(fields[17])
                            hdw_dict["maxbeams"] = int(fields[18])

                            hdw_info[hdw_dict["stid"]] = hdw_dict
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
