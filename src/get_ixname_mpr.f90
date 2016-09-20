module get_ixname

! Purpose
! Define functions to get the index of a named variable

USE nrtype                 ! variable types, etc.
USE public_var             ! Including common constant (physical constant, other e.g., missingVal, etc.)

implicit none

! the followings are public
public::get_ixDataMap       ! assign variable index to variables in mapping netCDF data
public::get_ixDataVeg       ! assign variable index to variables in veg netCDF data
public::get_ixDataSoil      ! assign variable index to variables in soil netCDF data
public::get_ixVarHru        ! assign variable index to model hru variables
public::get_ixVarTopo       ! assign variable index to topographic variables 
public::get_ixPrpVeg
! everything else
private

contains

! *******************************************************************************************************************
! function: get the index of the named variables for mapping data 
! *******************************************************************************************************************
 function get_ixDataMap(varName)
   USE var_lookup,only:ixVarMapData               ! indices of the named variables
   implicit none
   ! define dummy variables
   character(*), intent(in) :: varName            ! variable name
   integer(i4b)             :: get_ixDataMap      ! index of the named variable
   ! get the index of the named variables
   select case(trim(varName))
     case('hru_id');         get_ixDataMap = ixVarMapData%hru_id         ! hru ID 
     case('weight');         get_ixDataMap = ixVarMapData%weight         ! weight of overlapping polygon
     case('overlapPolyId');  get_ixDataMap = ixVarMapData%overlapPolyId  ! overlapping polygon ID 
     ! get to here if cannot find the variable
     case default;     get_ixDataMap = imiss
   endselect
 end function get_ixDataMap

! *******************************************************************************************************************
! function: get the index of the named variables for veg data 
! *******************************************************************************************************************
 function get_ixDataVeg(varName)
   USE var_lookup,only:ixVarVegData               ! indices of the named variables
   implicit none
   ! define dummy variables
   character(*), intent(in) :: varName            ! variable name
   integer(i4b)             :: get_ixDataVeg      ! index of the named variable
   ! get the index of the named variables
   select case(trim(varName))
     case('polyid');     get_ixDataVeg = ixVarVegData%polyid     ! veg polygon ID 
     case('lai');        get_ixDataVeg = ixVarVegData%lai        ! monthly lai [m2 m-2] 
     case('grnfrc');     get_ixDataVeg = ixVarVegData%grnfrc     ! monthly greeness fraction [-] 
     case('vegclass');   get_ixDataVeg = ixVarVegData%vegclass   ! veg class in veg polygon and layer
     ! get to here if cannot find the variable
     case default;     get_ixDataVeg = imiss
   endselect
 end function get_ixDataVeg

! *******************************************************************************************************************
! function: get the index of the named variables for soil data 
! *******************************************************************************************************************
 function get_ixDataSoil(varName)
  USE var_lookup,only:ixVarSoilData                 ! indices of the named variables
  implicit none
  ! define dummy variables
  character(*), intent(in) :: varName            ! variable name
  integer(i4b)             :: get_ixDataSoil     ! index of the named variable
  ! get the index of the named variables
  select case(trim(varName))
   case('polyid');       get_ixDataSoil = ixVarSoilData%polyid        ! soil polygon ID 
   case('hslyrs');       get_ixDataSoil = ixVarSoilData%hslyrs        ! soil layer thickness [m] 
   case('soilclass');    get_ixDataSoil = ixVarSoilData%soilclass     ! soil class in soil polygon and layer
   case('sand_frc');     get_ixDataSoil = ixVarSoilData%sand_frc      ! sand fraction in soil polygon and layer [-]
   case('silt_frc');     get_ixDataSoil = ixVarSoilData%silt_frc      ! silt fraction in soil polygon and layer [-]
   case('clay_frc');     get_ixDataSoil = ixVarSoilData%clay_frc      ! clay fraction in soil polygon and layer [-]
   case('bulk_density'); get_ixDataSoil = ixVarSoilData%bulk_density  ! bulk density in soil polygon and layer [-]
   case('ele_mean');     get_ixDataSoil = ixVarSoilData%ele_mean      ! average elev over soil polygon [m]
   case('ele_std');      get_ixDataSoil = ixVarSoilData%ele_STD       ! std elev over soil polygon [m]
   case('slp_mean');     get_ixDataSoil = ixVarSoilData%slp_mean      ! average slope over soil polygon [-]
   ! get to here if cannot find the variable
   case default;     get_ixdataSoil = imiss
  endselect
 end function get_ixDataSoil

! *******************************************************************************************************************
! Function: get the index of the named variables for the hru characteristics
! *******************************************************************************************************************
 function get_ixVarHru(varName)
  USE var_lookup,only:ixVarHru                  ! indices of the named variables
  implicit none
  ! define dummy variables
  character(*), intent(in) :: varName            ! variable name
  integer(i4b)             :: get_ixVarHru         ! index of the named variable
  ! get the index of the named variables
  select case(trim(varName))
   case('lat');    get_ixVarHru = ixVarHru%lat      ! latitude (degrees north)
   case('lon');    get_ixVarHru = ixVarHru%lon      ! longitude (degrees east)
   case('ele');    get_ixVarHru = ixVarHru%ele      ! elevation (m)
   case('ann_P');  get_ixVarHru = ixVarHru%ann_P    ! average annual precipitation (mm)
   case('avg_T');  get_ixVarHru = ixVarHru%avg_T    ! average annual temperatur e(C)
   case('july_T'); get_ixVarHru = ixVarHru%july_T   ! average July Temperature (C)
   ! get to here if cannot find the variable
   case default;   get_ixVarHru = imiss
  endselect
 end function get_ixVarHru

! *******************************************************************************************************************
! function: get the index of the named variables for vegetation properties 
! *******************************************************************************************************************
 function get_ixPrpVeg(varName)
  USE var_lookup,only:ixPrpVeg                   ! indices of the named variables
  implicit none
  ! define dummy variables
  character(*), intent(in) :: varName            ! variable name
  integer(i4b)             :: get_ixPrpVeg       ! index of the named variable
  ! get the index of the named variables
  select case(trim(varName))
   case('lai01');        get_ixPrpVeg = ixPrpVeg%lai01
   case('lai02');        get_ixPrpVeg = ixPrpVeg%lai02
   case('lai03');        get_ixPrpVeg = ixPrpVeg%lai03
   case('lai04');        get_ixPrpVeg = ixPrpVeg%lai04
   case('lai05');        get_ixPrpVeg = ixPrpVeg%lai05
   case('lai06');        get_ixPrpVeg = ixPrpVeg%lai06
   case('lai07');        get_ixPrpVeg = ixPrpVeg%lai07
   case('lai08');        get_ixPrpVeg = ixPrpVeg%lai08
   case('lai09');        get_ixPrpVeg = ixPrpVeg%lai09
   case('lai10');        get_ixPrpVeg = ixPrpVeg%lai10
   case('lai11');        get_ixPrpVeg = ixPrpVeg%lai11
   case('lai12');        get_ixPrpVeg = ixPrpVeg%lai12
   case('vegtype');      get_ixPrpVeg = ixPrpVeg%vegtype
   case('nroot');        get_ixPrpVeg = ixPrpVeg%nroot
   case('snup');         get_ixPrpVeg = ixPrpVeg%snup
   case('rs');           get_ixPrpVeg = ixPrpVeg%rs
   case('mrs');          get_ixPrpVeg = ixPrpVeg%mrs
   case('leafDim');      get_ixPrpVeg = ixPrpVeg%leafDim
   case('can_top_h');    get_ixPrpVeg = ixPrpVeg%can_top_h
   case('can_bot_h');    get_ixPrpVeg = ixPrpVeg%can_bot_h
   case('c_veg');        get_ixPrpVeg = ixPrpVeg%c_veg
   case('maxMassVeg');   get_ixPrpVeg = ixPrpVeg%maxMassVeg
   ! get to here if cannot find the variable
   case default;         get_ixPrpVeg = imiss
  endselect
 end function get_ixPrpVeg

end module get_ixname
