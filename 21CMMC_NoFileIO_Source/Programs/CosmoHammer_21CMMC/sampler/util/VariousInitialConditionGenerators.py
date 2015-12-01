
import numpy as np
import os
from decimal import *
import multiprocessing
import itertools

from ctypes import cdll
from ctypes import c_float,c_int,c_uint,c_double

SharedLibraryLocation = "%s/libdrive_21cmMC_streamlined.so"%(os.getcwd())

Lib_21CMFAST = cdll.LoadLibrary(SharedLibraryLocation)
Function_21CMFAST = Lib_21CMFAST.drive_21CMMC
Function_IntINIT = Lib_21CMFAST.IntFunction
Function_FloatINIT = Lib_21CMFAST.FloatFunction

Function_IntINIT.restype = c_int * 2
Function_FloatINIT.restype = c_float

Function_21CMFAST.restype = c_float * 2
Function_21CMFAST.argtypes = [np.ctypeslib.ndpointer(np.complex64, ndim=1,flags='aligned, contiguous'), 
                                np.ctypeslib.ndpointer(np.float32, ndim=1,flags='aligned, contiguous'), 
                                np.ctypeslib.ndpointer(np.float64, ndim=1,flags='aligned, contiguous'),
                                np.ctypeslib.ndpointer(np.float64, ndim=1,flags='aligned, contiguous'),
                                c_float, c_float, c_float, c_float, c_float, c_float, c_float, c_int]

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')

def ICposition(start_positions_individual,lowerbounds,upperbounds,paramValues,paramWidths,Redshift,DensityBoxes):

    # Function to generate suitable initial locations for the walkers

    # Random number generator for each CPU processor
    np.random.seed()

    # dummy arrays and float variables for unused arguments of the 21cmFast computation (i.e PS is not generated)
    dummy_array_float = np.zeros(10., dtype=np.float32)
    dummy_array = np.zeros(10.)
    dummy_float = 0.0

    # determine the new walker position is within the initial parameter bounds (yes, my python code looks like C!)
    for j in range(len(lowerbounds)):
        if (start_positions_individual[j] < lowerbounds[j] or start_positions_individual[j] > upperbounds[j]):
            new_start_parameter_logic = False
            while new_start_parameter_logic == False:
                new_start_parameter = paramValues[j]+3.*np.random.normal(size=1.0)*paramWidths[j]
                if (new_start_parameter > lowerbounds[j] and new_start_parameter < upperbounds[j]):
                    new_start_parameter_logic = True

            start_positions_individual[j] = new_start_parameter            

    # Once all walker positions are acceptable (within the parameter bounds) now check the neutral fraction. Accept walker positions with a non-zero neutral fraction
    # at the latest (lowest) redshift (this section does not need to be used)
    new_start_parameter_logic = False
    while new_start_parameter_logic == False:

        output = Function_21CMFAST(DensityBoxes,dummy_array_float,dummy_array,dummy_array,dummy_float,dummy_float,
                                        dummy_float,Redshift,start_positions_individual[0],start_positions_individual[1],start_positions_individual[2],0)

        nf_value = output[1]

        if nf_value == 0.:
            new_start_parameter_logic = False
        else:
            new_start_parameter_logic = True

        if new_start_parameter_logic == False:
            for j in range(len(lowerbounds)):
                start_parameter_logic_brandnew = False
                while start_parameter_logic_brandnew == False:
                    new_start_parameter = paramValues[j]+3.*np.random.normal(size=1.0)*paramWidths[j]
                    if (new_start_parameter > lowerbounds[j] and new_start_parameter < upperbounds[j]):
                        start_positions_individual[j] = new_start_parameter
                        start_parameter_logic_brandnew = True

    return start_positions_individual

def ICposition_star(all_arguments):

    return ICposition(*all_arguments)

class UniformPosition(object):
    """
        Generates samples in a very thight n-dimensional ball 
    """
    
    def __init__(self):
        """
            default constructor
        """
        pass

    def setup(self, sampler):
        """
            setup the generator
        """
        self.sampler = sampler
    
    def generate(self):
        """
            generates the positions
        """

        BOX_LEN = Function_FloatINIT()

        output = Function_IntINIT()

        HII_DIM = output[0]
        VelocityComponent = output[1]

        DensityBoxes = np.zeros(HII_DIM*HII_DIM*(HII_DIM/2+1),np.dtype('complex64'))
        IndividualBox = np.zeros(HII_DIM*HII_DIM*HII_DIM)

        f = open('../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc'%(self.sampler.Redshift,HII_DIM,BOX_LEN),'rb')  
        IndividualBox = np.fromfile(f, dtype = np.dtype('f'), count = HII_DIM*HII_DIM*HII_DIM)    
        f.close()

        for ii in range(HII_DIM):
            for jj in range(HII_DIM):
                for kk in range((HII_DIM/2)+1):
                    if (kk*2+1) < HII_DIM:
                        DensityBoxes[kk + ((HII_DIM/2)+1)*(jj+HII_DIM*ii)] = IndividualBox[(kk*2)+HII_DIM*(jj+HII_DIM*(ii))] + 1j*IndividualBox[(kk*2+1)+HII_DIM*((jj)+HII_DIM*(ii))]        

        # Note here, we do not need the velocity field box, just the density field. This is because we only calculate the neutral fraction of the box to check it is not already ionised
        # (we don't want an ionised box at a redshift where we expect a measurement. In principle this step is not needed as the walkers should move away from this location 
        # immediately during the sampling, but I include this just as an extra step. It has little impact)

        print('Generate Start Positions')
        start_positions = [self.sampler.paramValues+3.*np.random.normal(size=self.sampler.paramCount)*self.sampler.paramWidths for i in xrange(self.sampler.nwalkers)]

        pool = multiprocessing.Pool(self.sampler.threadCount)

        M = pool.map

        returned_list = list(M(ICposition_star,itertools.izip(start_positions, itertools.repeat(self.sampler.lowerbounds), itertools.repeat(self.sampler.upperbounds), itertools.repeat(self.sampler.paramValues), 
        itertools.repeat(self.sampler.paramWidths), itertools.repeat(self.sampler.Redshift), itertools.repeat(DensityBoxes))))

        print('Start Positions Generated')    
        return returned_list
    
    def __str__(self, *args, **kwargs):
        return "SampleBallPositionGenerator"
