module globalData

 use nrtype
 use public_var
 use data_type,  only: var_meta,par_meta,cpar_meta
 use var_lookup, only: nPar, nVarSoilData, nVarVegData, nVarMapData

implicit none
private

type(par_meta), save,              public   :: parMaster(nPar)               ! meta data for all the avaialble paramaeters for both gamma and beta
type(cpar_meta),      allocatable, public   :: parSubset(:)                  ! meta data for the parameters listed in 'CalPar' input
type(cpar_meta),      allocatable, public   :: gammaSubset(:)                ! meta data for gamma parameters listed in 'CalPar' input  
character(len=strLen),allocatable, public   :: betaInGamma(:)                ! name of beta parameters to be estimated based on 'CalPar' input
character(len=strLen),allocatable, public   :: betaNeeded(:)                 ! name of beta parameters to be estimated and their dependent beta parameters 

type(var_meta), save,              public   :: map_meta  (nVarMapData)
type(var_meta), save,              public   :: sdata_meta(nVarSoilData)
type(var_meta), save,              public   :: vdata_meta(nVarVegData)

end module globalData
