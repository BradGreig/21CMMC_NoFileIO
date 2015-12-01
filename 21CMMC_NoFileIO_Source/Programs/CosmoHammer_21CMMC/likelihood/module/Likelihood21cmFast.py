#!/usr/bin/env python

import os
import numpy as np
np.seterr(invalid='ignore', divide='ignore')
from decimal import *

from ctypes import cdll
from ctypes import c_float,c_int,c_uint,c_double

SharedLibraryLocation = "%s/libdrive_21cmMC_streamlined.so"%(os.getcwd())

Lib_21CMFAST = cdll.LoadLibrary(SharedLibraryLocation)
Function_21CMFAST = Lib_21CMFAST.drive_21CMMC

Function_21CMFAST.restype = c_float * 2
Function_21CMFAST.argtypes = [np.ctypeslib.ndpointer(np.complex64, ndim=1,flags='aligned, contiguous'), 
                                np.ctypeslib.ndpointer(np.float32, ndim=1,flags='aligned, contiguous'), 
                                np.ctypeslib.ndpointer(np.float64, ndim=1,flags='aligned, contiguous'),
                                np.ctypeslib.ndpointer(np.float64, ndim=1,flags='aligned, contiguous'),
                                c_float, c_float, c_float, c_float, c_float, c_float, c_float, c_int]

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

class Likelihood21cmFast_multiz(object):
    
    def __init__(self, DensityBoxes, Velocity_Boxes, PS_values, PS_Error, Redshift, Foreground_cut, Shot_Noise_cut, ModUncert):
        self.DensityBoxes = DensityBoxes
        self.Velocity_Boxes = Velocity_Boxes
        self.PS_values = PS_values
        self.PS_Error = PS_Error
        self.Redshift = Redshift
        self.Foreground_cut = Foreground_cut
        self.Shot_Noise_cut = Shot_Noise_cut
        self.ModUncert = ModUncert

    def Likelihood(self,ctx):
        params = ctx.getParams()

        nf_vals = np.zeros(len(self.Redshift))
        
        # determine the chi-squared linearly summed across redshift for all chosen redshift ranges        
        total_sum = 0;
        for i in range(len(self.Redshift)):

            # Run the C-library for 21cmFast
            output = Function_21CMFAST(self.DensityBoxes[i],self.Velocity_Boxes[i],self.PS_values[i],self.PS_Error[i],self.Foreground_cut,self.Shot_Noise_cut,
                                        self.ModUncert,self.Redshift[i],params[0],params[1],params[2],1)

            # Return and store the neutral fraction
            nf_vals[i] = output[1]

            # Return the chi-squared for this redshift and sum across multiple redshifts. If neutral fraction is zero (no 21cm PS is recovered), therefore return an arbitrarily 
            # large number to ensure the EMCEE sampler will not accept this parameter set as being an acceptable combination
            if nf_vals[i] == 0.:
                total_sum += 100.*sum(np.square((self.PS_values[i])/self.PS_Error[i]))
            else:                                
                total_sum += output[0]
       
        # Return the log-likelihood (summed over all redshifts) and the neutral fractions for each redshift selected
        return -0.5*total_sum,nf_vals

    def computeLikelihood(self, ctx):

        return self.Likelihood(ctx)

    def setup(self):
        print "Likelihood Fitting for 21cm Fast (Multi-z: 3 parameters, Zeta, Rmfp and Tvir)" 
