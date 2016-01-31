import os

def load_param(param_file=None):
    """ Method to load in the param file for a given run
        It will try and then ask for where the file is. if it doent know
    """
# Add a try catch statment incase you cant find the file

    param = {'file': param_file}

    if param['file'] is None:
        param['file'] = _get_param_file()

    fname = os.path.abspath(os.path.expandvars(param['file']))

    with open(fname) as f:
        content = f.readlines()

    for item in content:
        if '#define' in item and item[0] != '!':
            if len(item.split()) > 2:
                key = item.split()[1]
                val = item.split()[2]
                val = _convert(item.split()[2])
            else:
                key = item.split()[1]
                val = None

            param[key] = val

# An issue: In the param it is commen to say nchannels as pex. This
#           presents a problem in how we read the param file. So we
#           can run through the param dictionary and replace any value
#           with the coresponding key
    
    for key,val in param.iteritems():
        if param.has_key(val):
            param[key] = param[val]

    return param


def _get_param_file():
    fname = raw_input('Please Param File: ')
    fname = os.path.abspath(os.path.expandvars(fname))

    while not os.path.isfile(fname):
        error_text = '\nFile %s not found!\n' \
                     'Please Enter Param file: ' % fname

        fname = raw_input(error_text)
        fname = os.path.abspath(os.path.expandvars(fname))

    return fname


def _convert(val):
    constructors = [int, float, str]
    for c in constructors:
        try:
            return c(val)
        except ValueError:
            pass
