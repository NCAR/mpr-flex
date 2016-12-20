# Hydropt (to become MPR-flex)
Flexible hydrologic model parameter estimation driver. This program includes two main routines: 1) optimization and 2) mpr

# General overview
Three options to run this program (opt in namelist):
0-> Perform only MPR to estimate parameters listed in calPar.txt and restart file (for gamma parameter values), and output them in NetCDF.
1-> Perform model simulation (specified in idModel in namelist) with parameters listed in calPar.txt and restart file (for gamma parameter values) and default for the rest of parameters.
3-> Perform DDS calibration to estimate optimized parameters listed in calPar.txt. 

For opt=1, two options to calibrate 1) multiplier based calibration, 2) MPR calibration.
if model use spatially lumped mode, option 1) is equivalent to direct parameter adjustment 

Two inputs need to be prepared for all the run options, 1) calPar.txt 2)namelist.cal
1) list which parameters are calibrated (see input/calPar.txt)
2) configure/set run (see nml/namelist.cal)
   see comments for each nml variables for what kind of information need to be provided

To get parameter name in calPar.txt, See popMeta.f90 to see what kind of gamma and beta parameters are already defined.

if you include gamma parmaters in calPar, MPR use transfer function assoicated with calibrated gamma parameters to get associated beta parameter and output them in parameter file

Name convention of gamma parameters. 
<beta_parameter_name><integer1>gamma<integer2>
<beta_parameter_name> = beta parameter name assocated with this gamma parameter
<integer1> = <interger1>th transfer function form (most of transfer function have only one now) 
<integer2> = <interger2>th parameter in the transfer function (number of gamma parameters in transfer function vary) 
e.g.,ks1gamma1 is 1st gamma parameter for ks (saturated hydraulic conductivity) transfer function 1 (cosby equation) 

# Model specific I/O implementation
For run option 2 and 3, VIC and SAC are now implemented (as of 10/26/2016).

To use other model 
You will need to write model specific routines - use model_routines.f90.template
For VIC, vic_routines.f90.  For SAC, sac_routines.f90
Using this template, you need to write several subrutines
1. soil parameter reading routine
2. hru id reading routine
3. soil parameter adjustment routine (used for multiplier based calibration)
4. soil parameter replacement routine (used to replace a priori parameter with MPR derived parameter )
5. soil parameter writing routine
6. Other parameter adjustment routine e.g., snow, vege (for now, no MPR for this)

# Compiling the program
Compiled with pgi/15.7 and ifort/2015.2.164 (15.0.1)

# Terminologies
gamma parameters: parameter of transfer functions used in MPR
beta parameters:  model parameters
