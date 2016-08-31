module globalData

 use nrtype
 use data_type,  only: par_meta,cpar_meta
 use var_lookup, only: nGamma,nBeta

implicit none

private

type(par_meta),  public                :: gammaMaster(nGamma)
type(par_meta),  public                :: betaMaster(nBeta)
type(cpar_meta), allocatable, public   :: parSubset(:)

end module globalData

