import numpy as np


class FlagWork(object):

    def __init__(self, flagValues,flagMeanings):
        self.flagValues = flagValues ##numpy array of values
        self.flagMeanings = flagMeanings ##list of meanings
        self.dType = str(flagValues.dtype)

    def code(self, maskList):
        #myCode = np.uint64(0)
        myCode = np.array(0).astype(self.dType)
        for flag in maskList:
            myCode |= self.flagValues[self.flagMeanings.index(flag)]
        return myCode

    def mask(self, flags, maskList):
        myCode = self.code(maskList)
        flags = np.array(flags).astype(self.dType)
        #print flags
        #print myCode
        return np.bitwise_and(flags, myCode)

    def decode(self, val):
        val = np.array(val).astype(self.dType)
        count = 0
        res = []
        mask = np.zeros(len(self.flagMeanings))
        for value in self.flagValues :
            if value & val :
                res.append(self.maskNames[count])
                mask[count] = 1
            count += 1
        return res, mask



def get_info_from_flag_variable(variable,key_error='flag_functions'):
    flag_meanings = None
    if 'flag_meanings' in variable.ncattrs():
        flag_meanings = [x.strip() for x in variable.flag_meanings.split(' ')]
    elif 'flag_list' in variable.ncattrs():
        flag_meanings = [x.strip() for x in variable.flag_list.split(',')]
    if flag_meanings is None:
        print(f'[ERROR][{key_error}] Flag list in not available for variable {variable.name}, attribute flag_meanings or flag_list is required')

    flag_values_var = None
    if 'flag_values' in variable.ncattrs():
        flag_values_var = variable.flag_values
    elif 'flag_masks' in variable.ncattrs():
        flag_values_var = variable.flag_masks
    elif 'flag_mask' in variable.ncattrs():
        flag_values_var = variable.flag_mask
    if isinstance(flag_values_var, str):
        try:
            flag_values_var = [x for x in flag_values_var.split(',')]
        except ValueError as ex:
            print(f'[ERROR][{key_error}] Value error {ex} while parsing flag_values_var from variable {variable.name}')
            flag_values_var = None
    if isinstance(flag_values_var, np.ndarray):
       flag_values_var = flag_values_var.tolist()
    if not isinstance(flag_values_var, list):
        print(f'[ERROR][{key_error}] Flag values list is not available for variable {variable.name}, attribute flag_values, flag_masks or flag_mask with a comma separated list of integer values is required')
        flag_values_var = None

    return flag_meanings, flag_values_var


def start_flag_work_from_variable(variable):
    flag_meanings, flag_values_var = get_info_from_flag_variable(variable)
    if flag_meanings is None or flag_values_var is None:
        return None
    if np.issubdtype(variable[:].dtype, np.floating):  ##floating points are not allowed
        flag_values = np.array(flag_values_var).astype('uint64')
    else:
        flag_values = np.array(flag_values_var).astype(variable[:].dtype)
    fw = FlagWork(flag_values,flag_meanings)
    return fw