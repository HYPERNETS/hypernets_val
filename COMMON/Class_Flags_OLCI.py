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

    def __init__(self, pars={}):

         self.maskValues = np.array([1, 2, 4, 8, 8388608, 16777216, 16, 32, 64, 128, 256, 512, 1024,2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 2097152, 33554432, 67108864, 134217728, 268435456, 4294967296, 8589934592, 17179869184, 34359738368, 68719476736, 137438953472, 274877906944, 549755813888, 1099511627776, 2199023255552, 4398046511104, 8796093022208, 17592186044416, 35184372088832, 70368744177664, 140737488355328, 281474976710656, 562949953421312, 1125899906842624, 2251799813685248, 4503599627370496, 9007199254740992, 18014398509481984, 36028797018963968], dtype='uint64')
  
         self.maskNames = ["INVALID", "WATER", "LAND", "CLOUD", "CLOUD_AMBIGUOUS", "CLOUD_MARGIN", "SNOW_ICE", "INLAND_WATER", "TIDAL", "COSMETIC", "SUSPECT", "HISOLZEN", "SATURATED", "MEGLINT", "HIGHGLINT", "WHITECAPS", "ADJAC", "WV_FAIL", "PAR_FAIL", "AC_FAIL", "OC4ME_FAIL", "OCNN_FAIL", "KDM_FAIL", "BPAC_ON", "WHITE_SCATT", "LOWRW", "HIGHRW", "ANNOT_ANGSTROM", "ANNOT_AERO_B", "ANNOT_ABSO_D", "ANNOT_ACLIM", "ANNOT_ABSOA", "ANNOT_MIXR1", "ANNOT_DROUT", "ANNOT_TAU06", "RWNEG_O1", "RWNEG_O2", "RWNEG_O3", "RWNEG_O4", "RWNEG_O5", "RWNEG_O6", "RWNEG_O7", "RWNEG_O8", "RWNEG_O9", "RWNEG_O10", "RWNEG_O11", "RWNEG_O12", "RWNEG_O16", "RWNEG_O17", "RWNEG_O18", "RWNEG_O21"]
