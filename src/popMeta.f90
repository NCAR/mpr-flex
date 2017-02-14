module popMeta

implicit none

private

public::popMprMeta

contains

subroutine popMprMeta(err,message)
  use nrtype
  use data_type,  only:var_meta
  use var_lookup, only:ixVarMapData
  use var_lookup, only:ixVarSoilData
  use var_lookup, only:ixVarVegData
  use globalData, only:map_meta
  use globalData, only:sdata_meta
  use globalData, only:vdata_meta
  
  implicit none
  !output variable
  integer(i4b),intent(out)      :: err     ! error code
  character(*),intent(out)      :: message ! error message
  
  ! initialize error control
  err=0; message='popMprMeta/'
  ! Mapping data meta
  map_meta(ixVarMapData%hru_id)          = var_meta('hru_id'       ,"hru id"                            , "-"              ,"1D", "integer")
  map_meta(ixVarMapData%weight)          = var_meta('weight'       ,"areal weight of intersect polygon" , "-"              ,"1D", "integer")
  map_meta(ixVarMapData%overlapPolyId)   = var_meta('overlapPolyId',"id of intersect polygo"            , "-"              ,"1D", "integer")
  ! Soil data variables
  sdata_meta(ixVarSoilData%polyid)       = var_meta('polyid'       , "soil polygon id"                  ,  "-"             ,"1D", "integer")
  sdata_meta(ixVarSoilData%soilclass)    = var_meta('soilclass'    , "USDA soil class"                  ,  "-"             ,"2D", "integer")
  sdata_meta(ixVarSoilData%hslyrs)       = var_meta('hslyrs'       , "soil layer thickness"             ,  "m"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%sand_frc)     = var_meta('sand_frc'     , "sand percentage"                  ,  "%"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%silt_frc)     = var_meta('silt_frc'     , "silt percentage"                  ,  "%"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%clay_frc)     = var_meta('clay_frc'     , "clay percentage"                  ,  "%"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%bulk_density) = var_meta('bulk_density' , "bulk density"                     ,  "kg/m^3"        ,"2D", "double" )
  sdata_meta(ixVarSoilData%ele_mean)     = var_meta('ele_mean'     , "mean elevation"                   ,  "m"             ,"1D", "double" )
  sdata_meta(ixVarSoilData%ele_std)      = var_meta('ele_std'      , "std elevation"                    ,  "m"             ,"1D", "double" )
  sdata_meta(ixVarSoilData%slp_mean)     = var_meta('slp_mean'     , "mean slope"                       ,  "-"             ,"1D", "double" )
  ! Vege data variables 
  vdata_meta(ixVarVegData%polyid)        = var_meta('polyid'       , "vege polygon id"                  ,  "-"             ,"1D", "integer")
  vdata_meta(ixVarVegData%vegclass)      = var_meta('vegclass'     , "vegetation class"                 ,  "-"             ,"1D", "double" )
  vdata_meta(ixVarVegData%grnfrc)        = var_meta('grnfrc'       , "green fraction"                   ,  "-"             ,"1D", "double" )
  vdata_meta(ixVarVegData%lai)           = var_meta('lai'          , "monthly lai"                      ,  "m^2/m^2"       ,"2D", "double" )

end subroutine

end module popMeta 
