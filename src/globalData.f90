module globalData

 use nrtype
 use public_var
 use data_type,  only: par_meta,cpar_meta
 use var_lookup, only: nPar

implicit none
private

type(par_meta), save,              public   :: parMaster(nPar)
type(cpar_meta),      allocatable, public   :: parSubset(:)
type(cpar_meta),      allocatable, public   :: gammaSubset(:)
character(len=strLen),allocatable, public   :: betaInGamma(:) 

end module globalData

