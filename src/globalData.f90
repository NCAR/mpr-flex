module globalData

 use nrtype
 use public_var
 use data_type,  only: var_meta,par_meta,cpar_meta
 use var_lookup, only: nPar, nVarSoilData, nVarVegData, nVarMapData

implicit none
private

type(par_meta), save,              public   :: parMaster(nPar)
type(cpar_meta),      allocatable, public   :: parSubset(:)
type(cpar_meta),      allocatable, public   :: gammaSubset(:)
character(len=strLen),allocatable, public   :: betaInGamma(:) 

type(var_meta), save,              public   :: map_meta  (nVarMapData)
type(var_meta), save,              public   :: sdata_meta(nVarSoilData)
type(var_meta), save,              public   :: vdata_meta(nVarVegData)

end module globalData
