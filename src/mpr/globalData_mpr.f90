module globalData_mpr

 use nrtype
 use data_type,  only: var_meta
 use var_lookup, only: nVarHru
 use var_lookup, only: nVarMapData
 use var_lookup, only: nVarSoilData
 use var_lookup, only: nVarVeglData
 use var_lookup, only: nPrpSoil

implicit none

private

type(var_meta), public      :: hru_meta  (nVarHru)
type(var_meta), public      :: map_meta  (nVarMapData)
type(var_meta), public      :: sdata_meta(nVarSoilData)
type(var_meta), public      :: vdata_meta(nVarVegData)
type(var_meta), public      :: sprp_meta (nPrpSoil)

end module globalData_mpr

