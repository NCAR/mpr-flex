module globalData

 use nrtype
 use public_var
 use data_type,  only: var_meta, par_meta, cpar_meta, beta_meta, scale_meta, input_meta
 use var_lookup, only: nBeta, nGamma, nVarSoilData, nVarVegData, nVarMapData, nPrpVeg

implicit none

private

type(par_meta),       save,            public   :: gammaMaster(nGamma)           ! meta data for all the avaialble gamma paramaeters
type(par_meta),       save,            public   :: betaMaster(nBeta)             ! meta data for all the avaialble beta paramaeters
type(beta_meta),      save,            public   :: beta(nBeta)                   ! meta data for parameter dependency for all the beta parameters 
integer(i2b),         save,            public   :: betaOrder(nBeta)              ! number of vegetation parameters to be estimated with MPR
type(cpar_meta),      save,allocatable,public   :: parSubset(:)                  ! meta data for the parameters listed in 'CalPar' input
type(cpar_meta),      save,allocatable,public   :: gammaSubset(:)                ! meta data for gamma parameters associated with beta listed in 'CalPar' input  
character(len=strLen),save,allocatable,public   :: betaInGamma(:)                ! name of beta parameters to be estimated based on 'CalPar' input
character(len=strLen),save,allocatable,public   :: soilBetaInGamma(:)            ! name of soil beta parameters to be estimated based on 'CalPar' input
character(len=strLen),save,allocatable,public   :: vegBetaInGamma(:)             ! name of veg beta parameters to be estimated based on 'CalPar' input
character(len=strLen),save,allocatable,public   :: betaNeeded(:)                 ! name of beta parameters to be estimated and their dependent beta parameters 
type(input_meta),     save,allocatable,public   :: calParMeta(:)                 ! meta data for all the avaialble paramaeters for both gamma and beta
type(scale_meta),     save,allocatable,public   :: betaCalScale(:)               ! meta data for beta paramets whose scaling operator is calibrated

integer(i2b),         save,            public   :: nBetaGammaCal
integer(i2b),         save,            public   :: nBetaCal
integer(i2b),         save,            public   :: nGammaCal
integer(i2b),         save,            public   :: nSoilParModel                 ! number of soil parameters to be estimated with MPR
integer(i2b),         save,            public   :: nVegParModel                  ! number of vegetation parameters to be estimated with MPR

type(var_meta),       save,            public   :: map_meta  (nVarMapData)
type(var_meta),       save,            public   :: sdata_meta(nVarSoilData)
type(var_meta),       save,            public   :: vdata_meta(nVarVegData)
type(var_meta),       save,                  public   :: vprp_meta (nPrpVeg)

end module globalData
