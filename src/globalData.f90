module globalData

 use nrtype
 use public_var
 use data_type,  only: var_meta, par_meta, cpar_meta, beta_meta, scale_meta, input_meta
 use var_lookup, only: nBeta, nGamma, nVarSoilData, nVarVegData, nVarMapData, nPrpVeg

implicit none

private

! master meta data including all available parameters
type(par_meta),       save,            public   :: gammaMaster(nGamma)       ! meta data for all the avaialble gamma paramaeters
type(par_meta),       save,            public   :: betaMaster(nBeta)         ! meta data for all the avaialble beta paramaeters
type(beta_meta),      save,            public   :: beta(nBeta)               ! meta data for parameter dependency for all the beta parameters 
! Based on a list of beta parameters to be estimated from CalPar input
type(cpar_meta),      save,allocatable,public   :: parSubset(:)              ! meta data for beta and gamma parameters calibrated based on 'CalPar' input - direct calibration of beta
type(cpar_meta),      save,allocatable,public   :: gammaSubset(:)            ! subset of parSubset including only gamma parameters. 
type(input_meta),     save,allocatable,public   :: calParMeta(:)             ! meta data for speficification of beta parameteters in "CalPar" input
type(scale_meta),     save,allocatable,public   :: betaCalScale(:)           ! meta data for scaling operator of beta paramets 
integer(i2b),         save,            public   :: betaOrder(nBeta)          ! index of beta parameters sorted in computing order (only beta parameter estimated with MPR) 
character(len=strLen),save,allocatable,public   :: betaInGamma(:)            ! name of beta parameters to be estimated with MPR based on 'CalPar' input
character(len=strLen),save,allocatable,public   :: soilBetaInGamma(:)        ! subset of "betaInGamma"  including only soil parameters 
character(len=strLen),save,allocatable,public   :: vegBetaInGamma(:)         ! subset of "betaInGamma" including only vegetation parameters
real(dp),             save,allocatable,public   :: parArray(:,:)             ! calibrating parameter array - input for optimization routine
logical(lgc),         save,allocatable,public   :: parMask(:)                ! calibrating parameter mask vector - input for optimization routine
! Number of parameters 
integer(i2b),         save,            public   :: nBetaDirCal               ! number of beta parameters to be directly calibrated
integer(i2b),         save,            public   :: nGammaCal                 ! number of gamma parameters to be calibrated
integer(i2b),         save,            public   :: nBetaGammaCal             ! sum of beta parameter directly calibrated and gamma parameters calibrated 
integer(i4b),         save,            public   :: nParCalSum                ! Total number of calibrating parameters = n(gamma) + m1(beta) + m2(beta_per_layer) * nLyr (m1+m2=nBetaDirCal)
integer(i4b),         save,            public   :: nBetaNeed                 ! number of beta parameters computed in MPR including calibrating beta and their dependent beta 
integer(i2b),         save,            public   :: nSoilParModel             ! number of soil parameters to be estimated with MPR
integer(i2b),         save,            public   :: nVegParModel              ! number of vegetation parameters to be estimated with MPR
! meta data for input data 
type(var_meta),       save,            public   :: map_meta  (nVarMapData)   ! mapping data
type(var_meta),       save,            public   :: sdata_meta(nVarSoilData)  ! soil data
type(var_meta),       save,            public   :: vdata_meta(nVarVegData)   ! vegetation data
type(var_meta),       save,            public   :: vprp_meta (nPrpVeg)       ! mapping data

end module globalData
