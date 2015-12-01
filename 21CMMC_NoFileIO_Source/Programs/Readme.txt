This is a (likely, highly incomplete) read me file for 21CMMC.

21CMMC is not installed on one’s computer (currently). I have modified
CosmoHammer from the base code and so far have not changed (or tested) ‘building’
my version of CosmoHammer. Therefore clashes might occur if you already use 
CosmoHammer.

To run 21CMMC, you need to compile the C code. The folder includes a makefile for the
code. Most of the code is standard 21cmFAST C code, with the addition of the 
“drive_21cmMC_streamlined.c” file, which is what 21CMMC calls. It is a thinned out version 
of the EoR 21cmFAST code (mostly find_HII_bubbles.c and delta_T.c).

I have included in the folder “Boxes” z = 7,8,9 and 10 density and velocity field boxes from
21cmFAST. These were used in previous publications. Additionally, these have been used to 
generate the sensitivity curves using 21cmSense. Note: These boxes were generated using WMAP 
cosmology. To use Planck cosmology, you will need to generate new boxes using 21cmFAST (the 
same is true if you want to sample the PS at a redshift other than those above).

The boxes are 128^3 250 Mpc boxes. These were 768^3 boxes smoothed onto 128^3 boxes as performed
in INIT.H. These boxes are used for sampling within 21CMMC. To generate the sensitivity curves, 
I use the 21cmPS from a 1536^3 500 Mpc box, smoothed onto a 256^3 grid. Through testing, I found
these ratios to be the best. Therefore I recommend using these numbers for generating your own 
boxes for your own cosmology, redshifts or any other purposes.

*** Likelihood computation in 21CMMC ***

The likelihood computation (the main work in 21CMMC) is performed within Likelihood21cmFast.py, which is
located in “CosmoHammer/likelihood/module”. The likelihood computation is a chi^2 likelihood of the 21cm 
PS given some mock observation.

*** Generating Error (sensitivity curves) ***

Telescope sensitivities are not generated within 21CMMC, but rather from a separate code 
(21cmSense, https://github.com/jpober/21cmSense). The output noise files from 21cmSense are
read into 21CMMC.py and passed into the MCMC sampler. I have included some noise curves in 
“NoiseData/“. Note, these PS are binned exactly as those sampled for 250 Mpc^3 cubed boxes. At 
present there is no interpolation within 21CMMC of the 21cm PS, but this will change in a future
release.

The noise calibration files are also provided from an optimised SKA and 331 antennae HERA instrument, 
which are output from 21cmSense. These .npz files should only be used with calc_sense.py, to generate
sensitivity curves. To generate your own sensitivity files, I refer you to the 21cmSense readme files.

*** Notes on 21CMMC ***

Unfortunately, 21CMMC is not completely optimised, nor completely user friendly. Therefore, at 
present it likely will take extra time to get used to using/modifying 21CMMC.

Within 21CMMC.py, everything is set for running 21CMMC. Here you can set redshifts, telescope noise curves,
modelling uncertainties, parameter ranges and a few other things. Note: At present, 21CMMC is hard coded to
work only for 3 parameters (Zeta, R_mfp and Tvir) or 2 parameters (Zeta and Tvir).

This hard coding is a major issue for modifying the code for user defined purposes.

The public release will have all these removed

*** Potential issues with 21CMMC ***

I have noted one issue with 21CMMC that has arisen once or twice. The error occurs occasionally and causes
the code to crash. At a point in the MCMC sampling, because I am reading the PS files from file, occasionally
it cannot find one and hence terminates the program. In future, I will remove all file I/O from within the 
likelihood sampler, which will alleviate this issue. Additionally, in future all command line calls will be 
removed.

*** Modifying 21CMMC ***

Unfortunately due to the hard coded nature of 21CMMC and its unfinished state, modifying 21CMMC will take more
effort than the final release should require.

Most necessary changes to 21CMMC will need to be performed to drive_21cm_streamlined.c. This is the major component
of 21CMMC, where all the computations are done. Here, the likelihood statistic can be changed, additional parameters
can be added and many other things. Any changes to drive_21cm_streamlined.c will then need to be then incorporated into
Likelihood21cmFast.py and likely throughout 21CMMC.py, ensemble.py (CosmoHammer_21CMMC/emcee), 
CosmoHammer.py (CosmoHammer_21CMMC/sampler) and VariousInitialConditionGenerators.py (CosmoHammer_21CMMC/sampler/util)

Apologies for the completely opaque and lack of flexibility in the current version of 21CMMC!

*** Plotting output of 21CMMC ***

I have included the module triangle.py and Plotting_Triangle.py which should be able to provide a triangle plot for the 
MCMC data. Note the copyright within triangle.py. This module I downloaded from the web. This is not what was used to 
generate figures within the 21CMMC paper.