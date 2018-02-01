import os  # REVIEW #7 License
import ConfigParser as cp


# REVIEW #1 docstring
def parse_hdw_files(hdw_files_path):

    hdw_info = {}

    loop_ran = False
    try:  # REVIEW #37 This is a huuuuuuge try statement for checking only the directory. Should be refactored
        for filename in os.listdir(hdw_files_path):  # REVIEW You can put the filenames in a list and then iterate over them
            loop_ran = True

            if "hdw.dat" in filename:
                with open(os.path.join(hdw_files_path, filename), 'r') as f:  # REVIEW #26 can name 'f' something more meaningful, like 'hdw_file'
                    backwards_lines = reversed(f.readlines())  # REVIEW #3 just a comment about why you start from the end of the file? OR there may be a better way, where you just read the one line you want? (Not necessarily the last line)

                    for line in backwards_lines:
                        if '#' in line or not line.strip():
                            continue
                        else:
                            fields = line.split()
                            if len(fields) != 19:  # REVIEW #29 maagic number, put constant at top of file for number of fields?
                                raise ValueError(" {} is an invalid number "
                                                 "of hdw.dat fields!".format(len(fields)))

                            hdw_dict = {}

                            hdw_dict["stid"] = int(fields[0])
                            hdw_dict["year"] = int(fields[1])
                            hdw_dict["sc"] = int(fields[2])
                            hdw_dict["lat"] = float(fields[3])
                            hdw_dict["lon"] = float(fields[4])
                            hdw_dict["altitude"] = float(fields[5])
                            hdw_dict["boresight"] = float(fields[6])
                            hdw_dict["beamsep"] = float(fields[7])
                            hdw_dict["velsign"] = float(fields[8])  # REVIEW #28 should be integer I think? always +/- 1
                            hdw_dict["rxattenstep"] = float(fields[9])
                            hdw_dict["tdiff"] = float(fields[10])
                            hdw_dict["phasesign"] = float(fields[11])  # REVIEW #28 should be integer I think? always +/- 1
                            hdw_dict["interoffx"] = float(fields[12])
                            hdw_dict["interoffy"] = float(fields[13])
                            hdw_dict["interoffz"] = float(fields[14])
                            hdw_dict["rxrisetime"] = float(fields[15])
                            hdw_dict["attenstages"] = float(fields[16])  # REVIEW #28 I don't think it's possible to have a non-integer number of atten stages
                            hdw_dict["maxgates"] = int(fields[17])
                            hdw_dict["maxbeams"] = int(fields[18])

                            hdw_info[hdw_dict["stid"]] = hdw_dict  # REVIEW #0 Some radars have multiple lines in the hdw.dat files (example: ade) this will only grab the last line, which is the most up-to-date line. What if data from prior to the update needs to be parsed? We will have the incorrect hdw info.
                            break  # REVIEW #0 it might be good enough to just append to the dictionary, and don't have the break here. Then whenever the hdw info is used, check the dates of each entry for that radar to make sure that you're using the correct hdw info
    except OSError as e:  # REVIEW #15 There are two failure modes we found 1) Permission denied or 2) no such file or directory. Yoou should wrap only the call to listdir in the try catch block
        error_msg = "No directory {0} found when locating hdw files!".format(hdw_files_path)
        raise ImportError(error_msg)

    if loop_ran is False or len(hdw_info) == 0:  # REVIEW #37 the for loop above can have an else statement for this probably? Instead of keeping track of a boolean
        error_msg = "No hardware files found at {0}".format(hdw_files_path)
        raise ImportError(error_msg)

    return hdw_info

# REVIEW #1 docstring
def parse_config_file():
    config = cp.ConfigParser()

    valid_sections = ["core", "fitacf"]
    for loc in os.curdir, os.path.expanduser("~"), "/etc/backscatter":  # REVIEW #29 magic string, place at top? REVIEW #26 'loc' can be renamed to 'location' or similar

        file_path = os.path.join(loc, "backscatter.ini")  # REVIEW #29 magic string
        config.read(file_path)

        if config.sections():
            if set(config.sections()).issubset(set(valid_sections)):  # REVIEW #3 This line is difficult to understand at first so it could use a small comment
                return config
            else:
                raise ValueError("Config file is missing required sections")
    else:  # REVIEW #0 This else statement feels like it should be part of the 'if config.sections()' statement...
        raise IOError("Config file missing!")
