#! /usr/bin/env python
'''
$Id$

Introduction:
'''
import numpy as np
import numpy.ma as ma


class Class_Flags_OLCI(object):

    def Code(self, maskList):
        myCode = np.uint64(0)
        for flag in maskList:
            myCode |= self.maskValues[self.maskNames.index(flag)]
        return myCode

    def Mask(self, flags, maskList):
        myCode = self.Code(maskList)
        flags = np.uint64(flags)
        # print(myCode,flags)
        # print myCode
        return np.bitwise_and(flags, myCode)

    def Decode(self, val):
        count = 0
        res = []
        mask = np.zeros(len(self.maskValues))
        for value in self.maskValues:
            if value & val:
                res.append(self.maskNames[count])
                mask[count] = 1
            count += 1
        return res, mask

    def __init__(self, flagMasks, flagMeanings):
        self.maskValues = flagMasks
        self.maskNames = flagMeanings.split(' ')

        # print(type(self.maskValues))
        # print(type(self.maskNames))


class Class_Flags_Idepix(object):
    def __init__(self, flagMasks, flagMeanings):
        self.maskValues = flagMasks
        flagMeanings = flagMeanings.replace('  ', ' ')
        self.maskNames = flagMeanings.split(' ')

    def Code(self, maskList):
        myCode = np.int16(0)
        for flag in maskList:
            myCode |= self.maskValues[self.maskNames.index(flag)]
        return myCode

    def Mask(self, flags, maskList):
        myCode = self.Code(maskList)
        flags = np.int16(flags)
        # print flags
        # print myCode
        return np.bitwise_and(flags, myCode)

    def Decode(self, val):
        count = 0
        res = []
        mask = np.zeros(len(self.maskValues))
        for value in self.maskValues:
            if value & val:
                res.append(self.maskNames[count])
                mask[count] = 1
            count += 1
        return (res, mask)


class Class_Flags_Polymer(object):

    def __init__(self, flagMasks, flagMeanings):
        self.maskValues = flagMasks
        flagMeanings = flagMeanings.replace('  ', ' ')
        self.maskNames = flagMeanings.split(' ')

    def MaskGeneral(self, flags):
        flags = np.int64(flags)
        code = np.empty(flags.shape, dtype=np.int64)
        code[:] = 1023
        # print flags
        # print myCode
        return np.bitwise_and(flags, code)

    def Code(self, maskList):
        myCode = np.int64(0)
        for flag in maskList:
            myCode |= self.maskValues[self.maskNames.index(flag)]
        return myCode

    def Mask(self, flags, maskList):
        myCode = self.Code(maskList)
        flags = np.int64(flags)
        # print flags
        # print myCode
        return np.bitwise_and(flags, myCode)

    def Decode(self, val):
        count = 0
        res = []
        mask = np.zeros(len(self.maskValues))
        for value in self.maskValues:
            if value & val:
                res.append(self.maskNames[count])
                mask[count] = 1
            count += 1
        return (res, mask)
