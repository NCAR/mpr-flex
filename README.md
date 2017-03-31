# MPR-flex
Flexible hydrologic model parameter estimation driver. This program includes two main routines: 1) optimization and 2) mpr

## General overview
Three options to run this program (opt in namelist):
* 1-> Perform calibration to estimate optimized parameters listed in calPar.txt<sup>a</sup>. 
* 2-> Perform model simulation with calibrated beta multipliers and gamma parameters<sup>b</sup>.  
* 3-> Perform only MPR to estimate beta parameters listed in calPar.txt<sup>c</sup> and output them in NetCDF.

 <sup>a</sup> Model specified in "idModel" in ./nml/namelist.opt (See Model Options section). Use Dynamically Dimensional Search algorithm. Two options speficied in the 2nd column ("calMethod" column) in calPar.txt to calibrate I) multiplier of a priori beta parameter (calMethod=3), II) gamma parameter in MPR framework (calMethod=1). if model use spatially lumped mode, calMethod=3 is equivalent to direct parameter adjustment. a priori beta parameter file has to be read from model specific I/O. If only scaling parameter calibration desired, use calMethod=2 to avoild calibrating its gamma parameters.
 
 <sup>b</sup> This is supposed to be used after option 1. The parameter values listed in restart file ("restart_file" in ./nml/namelist.opt) and values read from a priori beta parameter file for the rest of parameters are used for the simulations. The list of parameter values in "restart_file" must be match up with a list of beta and gamma parameters and scaling parameter list expanded from calPar.txt. 

 <sup>c</sup> Only beta parameters with 1 or 2 in "calMethod" in calPar.txt. This option outputs parameter list in a text file specified in "mpr_param_file" in ./nml/namelist.opt. 

Two inputs need to be prepared for all the run options, 1) calPar.txt and 2) namelist.cal

1. list which parameters are calibrated (see input/calPar.txt for further option)
2. configure/set run (see nml/namelist.cal)
   see comments for each nml variables for what kind of information need to be provided

To get parameter name in calPar.txt, See popMeta.f90 to see what kind of gamma and beta parameters are already defined.

if you include gamma parmaters in calPar, MPR use transfer function assoicated with calibrated gamma parameters to get associated beta parameter and output them in parameter file

Name convention of gamma parameters. 

`<beta_parameter_name>``<integer1>`gamma`<integer2>`
* `<beta_parameter_name>` = beta parameter name assocated with this gamma parameter
* `<integer1>` = `<interger1>`th transfer function form (most of transfer function have only one now) 
* `<integer2>` = `<interger2>`th parameter in the transfer function (number of gamma parameters in transfer function vary)

e.g.,ks1gamma1 is 1st gamma parameter for ks (saturated hydraulic conductivity) transfer function type 1 (cosby equation) 

## Model options
0. No model (used for parameter estimation with gamma parameters specified in "restart_file")
1. VIC
2. SAC/SNOW17

## Model specific I/O implementation
For run option 1 and 2, VIC and SAC are now implemented (as of 10/26/2016).

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

## Compiling the program
pgi/15.7 or later
ifort/2015.2.164 (15.0.1) or later
gfortran 6.0.1 or later

## Terminologies
* gamma parameter: parameter of transfer functions used in MPR.
* beta parameters: model parameters.
* scaling parameter: coefficient of power mean. 
