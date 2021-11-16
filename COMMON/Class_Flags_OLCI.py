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
        #print flags
        #print myCode
        return np.bitwise_and(flags, myCode)

    def Decode(self, val):
        count = 0
        res = []
        mask = np.zeros(len(self.maskValues))
        for value in self.maskValues :
            if value & val :
                res.append(self.maskNames[count])
                mask[count] = 1
            count += 1
        return (res, mask)

    def __init__(self, flagMasks,flagMeanings):
        self.maskValues = flagMasks
        self.maskNames = flagMeanings.split(' ')
        
        #print(type(self.maskValues))
        #print(type(self.maskNames))