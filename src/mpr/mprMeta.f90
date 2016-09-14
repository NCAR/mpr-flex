module mprMeta 

implicit none

private

public::popMeta

contains

  subroutine popMeta(err,message)
  
    use nrtype
    use data_type,  only:var_meta
    use var_lookup, only:ixVarHru
    use var_lookup, only:ixVarMap
    use var_lookup, only:ixVarSoil
    use var_lookup, only:ixPrpSoil
    use var_lookup, only:ixVarVeg
    use globalData, only:hru_meta
    use globalData, only:map_meta
    use globalData, only:sdata_meta
    use globalData, only:sprp_meta
    use globalData, only:vdata_meta
    
    implicit none
  
    !output variable
    integer(i4b),intent(out)      :: err     ! error code
    character(*),intent(out)      :: message ! error message
    
    ! initialize error control
    err=0; message='popMeta/'

    ! model hru variables 
    hru_meta(ixVarHru%lat)           = var_meta('lat'       ,"latitude of hru centroid"          , "decimal degree" ,"1D", "double" )
    hru_meta(ixVarHru%lon)           = var_meta('lon'       ,"longitude of hru centroid"         , "decimal degree" ,"1D", "double" )
    hru_meta(ixVarHru%ele)           = var_meta('ele'       ,"hru mean elevation "               , "m"              ,"1D", "double" )
    hru_meta(ixVarHru%ann_P)         = var_meta('ann_P'     ,"mean annual total precipitation"   , "mm"             ,"1D", "double" )
    hru_meta(ixVarHru%avg_T)         = var_meta('avg_T'     ,"mean annual air temperature"       , "Celcius degree" ,"1D", "double" )
    hru_meta(ixVarHru%july_T)        = var_meta('july_T'    ,"mean July air temperature"         , "Celcius degree" ,"1D", "double" )
    ! Mapping data meta
    map_meta(ixVarMap%hruid)         = var_meta('hruid'     ,"hru id"                            , "-"              ,"1D", "integer")
    map_meta(ixVarMap%weight)        = var_meta('weight'    ,"areal weight of intersect polygon" , "-"              ,"1D", "integer")
    map_meta(ixVarMap%polyid)        = var_meta('polyid'    ,"id of intersect polygo"            , "-"              ,"1D", "integer")
    ! Soil data variables
    sdata_meta(ixVarSoil%polyid)     = var_meta('polyid'    , "soil polygon id"                  ,  "-"             ,"1D", "integer")
    sdata_meta(ixVarSoil%soilclass)  = var_meta('soilclass' , "USDA soil class"                  ,  "-"             ,"2D", "integer")
    sdata_meta(ixVarSoil%h)          = var_meta('h'         , "soil layer thickness"             ,  "m"             ,"2D", "double" )
    sdata_meta(ixVarSoil%sand)       = var_meta('sand'      , "sand percentage"                  ,  "%"             ,"2D", "double" )
    sdata_meta(ixVarSoil%silt)       = var_meta('silt'      , "silt percentage"                  ,  "%"             ,"2D", "double" )
    sdata_meta(ixVarSoil%clay)       = var_meta('caly'      , "clay percentage"                  ,  "%"             ,"2D", "double" )
    sdata_meta(ixVarSoil%BD)         = var_meta('BD'        , "bulk density"                     ,  "kg/m^3"        ,"2D", "double" )
    sdata_meta(ixVarSoil%ele_mean)   = var_meta('ele_mean'  , "mean elevation"                   ,  "m"             ,"1D", "double" )
    sdata_meta(ixVarSoil%ele_std)    = var_meta('ele_std'   , "std elevation"                    ,  "m"             ,"1D", "double" )
    sdata_meta(ixVarSoil%slp_mean)   = var_meta('slp_mean'  , "mean slope"                       ,  "-"             ,"1D", "double" )
    ! soil physical properties
    sprp_meta(ixPrpSoil%bd)          = var_meta('bd'        , "bulk density"                     , "kg/m^3"         ,"2D", "double" )
    sprp_meta(ixPrpSoil%ks)          = var_meta('ks'        , "saturated hydraulic conductivity" , "mm/s"           ,"2D", "double" )
    sprp_meta(ixPrpSoil%phi)         = var_meta('phi'       , "porosity"                         , "frac"           ,"2D", "double" )
    sprp_meta(ixPrpSoil%b)           = var_meta('b'         , "slope of retention curve"         , "-"              ,"2D", "double" )
    sprp_meta(ixPrpSoil%psis)        = var_meta('psis'      , "saturated matric potential"       , "kPa"            ,"2D", "double" )
    sprp_meta(ixPrpSoil%fc)          = var_meta('fc'        , "field capacity"                   , "frac"           ,"2D", "double" )
    sprp_meta(ixPrpSoil%wc)          = var_meta('wc'        , "wilting point"                    , "frac"           ,"2D", "double" )
    sprp_meta(ixPrpSoil%myu)         = var_meta('myu'       , "specific yield"                   , "-"              ,"2D", "double" )
    sprp_meta(ixPrpSoil%z)           = var_meta('z'         , "layer interface depth"            , "m"              ,"2D", "double" )
    sprp_meta(ixPrpSoil%h)           = var_meta('h'         , "layer thickness"                  , "m"              ,"2D", "double" )
    ! Vege data variables 
    vdata_meta(ixVarVeg%polyid)      = var_meta('polyid'    , "vege polygon id"                  ,  "-"             ,"1D", "integer")
    vdata_meta(ixVarVeg%vegclass)    = var_meta('vegclass'  , "vegetation class"                 ,  "-"             ,"1D", "double" )
    vdata_meta(ixVarVeg%grnfrc)      = var_meta('grnfrc'    , "green fraction"                   ,  "-"             ,"1D", "double" )
    vdata_meta(ixVarVeg%lai)         = var_meta('lai'       , "monthly lai"                      ,  "m^2/m^2"       ,"2D", "double" )

  end subroutine popMeta

end module mprMeta 
