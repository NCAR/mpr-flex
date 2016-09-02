module globalData

 use nrtype
 use data_type,  only: par_meta,cpar_meta
 use var_lookup, only: nPar

implicit none

private

type(par_meta),  public                :: parMaster(nPar)
type(cpar_meta), allocatable, public   :: parSubset(:)

end module globalData

