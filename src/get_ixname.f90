module get_ixname

! Purpose
! Define functions to get the index of a named variable

use nrtype                  ! variable types, etc.
use public_var

implicit none

! the followings are public
public::get_ixPar     ! assign variable index to parameter 
public::get_ixDataMap       ! assign variable index to variables in mapping netCDF data
public::get_ixDataVeg       ! assign variable index to variables in veg netCDF data
public::get_ixDataSoil      ! assign variable index to variables in soil netCDF data
public::get_ixVarHru        ! assign variable index to model hru variables
public::get_ixPrpVeg
! everything else
private

contains

! *******************************************************************************************************************
! function: get the index of the named variables for gamma parameters 
! *******************************************************************************************************************
 function get_ixPar(varName)
  use var_lookup,only:ixPar             ! indices of the named variables
  implicit none
  ! define dummy variables
  character(*), intent(in) :: varName     ! variable name
  integer(i4b)             :: get_ixPar   ! index of the named variable

  ! get the index of the named variables
  select case(trim(varName))
    case('ks1gamma1');        get_ixPar = ixPar%ks1gamma1       ! 
    case('ks1gamma2');        get_ixPar = ixPar%ks1gamma2       ! 
    case('ks1gamma3');        get_ixPar = ixPar%ks1gamma3       ! 
    case('ks2gamma1');        get_ixPar = ixPar%ks2gamma1       ! 
    case('ks2gamma2');        get_ixPar = ixPar%ks2gamma2       ! 
    case('ks2gamma3');        get_ixPar = ixPar%ks2gamma3       ! 
    case('phi1gamma1');       get_ixPar = ixPar%phi1gamma1      ! 
    case('phi1gamma2');       get_ixPar = ixPar%phi1gamma2      ! 
    case('phi1gamma3');       get_ixPar = ixPar%phi1gamma3      ! 
    case('phi2gamma1');       get_ixPar = ixPar%phi2gamma1      ! 
    case('phi2gamma2');       get_ixPar = ixPar%phi2gamma2      ! 
    case('phi2gamma3');       get_ixPar = ixPar%phi2gamma3      ! 
    case('phi2gamma4');       get_ixPar = ixPar%phi2gamma4      ! 
    case('phi2gamma5');       get_ixPar = ixPar%phi2gamma5      ! 
    case('phi2gamma6');       get_ixPar = ixPar%phi2gamma6      ! 
    case('fc1gamma1');        get_ixPar = ixPar%fc1gamma1       ! 
    case('wp1gamma1');        get_ixPar = ixPar%wp1gamma1       ! 
    case('b1gamma1');         get_ixPar = ixPar%b1gamma1        ! 
    case('b1gamma2');         get_ixPar = ixPar%b1gamma2        ! 
    case('b1gamma3');         get_ixPar = ixPar%b1gamma3        ! 
    case('psis1gamma1');      get_ixPar = ixPar%psis1gamma1     ! 
    case('psis1gamma2');      get_ixPar = ixPar%psis1gamma2     ! 
    case('psis1gamma3');      get_ixPar = ixPar%psis1gamma3     ! 
    case('myu1gamma1');       get_ixPar = ixPar%myu1gamma1      ! 
    case('myu1gamma2');       get_ixPar = ixPar%myu1gamma2      ! 
    case('z1gamma1');         get_ixPar = ixPar%z1gamma1        ! Layer bottom depth [m] 
    case('h1gamma1');         get_ixPar = ixPar%h1gamma1        ! top layer thickness [m] 
    case('h1gamma2');         get_ixPar = ixPar%h1gamma2        ! 2nd layer thickness [m] 
    case('binfilt1gamma1');   get_ixPar = ixPar%binfilt1gamma1  ! variable infilitration curve parameter 
    case('binfilt1gamma2');   get_ixPar = ixPar%binfilt1gamma2  ! variable infilitration curve parameter 
    case('D11gamma1');        get_ixPar = ixPar%D11gamma1       ! fraction of Dsmax where nonlinear baseflow begins 
    case('D21gamma1');        get_ixPar = ixPar%D21gamma1       ! fraction of maximum soil moisture where nonlinear baseflow begins 
    case('D31gamma1');        get_ixPar = ixPar%D31gamma1       ! fraction of maximum soil moisture where nonlinear baseflow begins 
    case('D41gamma1');        get_ixPar = ixPar%D41gamma1       ! fraction of maximum soil moisture where nonlinear baseflow begins 
    case('exp1gamma1');       get_ixPar = ixPar%exp1gamma1      ! exponent in Campbell equatin for Kh
    case('exp1gamma2');       get_ixPar = ixPar%exp1gamma2      ! exponent in Campbell equatin for Kh
    case('Ksat1gamma1');      get_ixPar = ixPar%ksat1gamma1     ! saturated hydrologic conductivity [mm/day]
    case('bbl1gamma1');       get_ixPar = ixPar%bbl1gamma1      ! bubbling pressure of soil [cm] 
    case('bbl1gamma2');       get_ixPar = ixPar%bbl1gamma2      ! bubbling pressure of soil [cm] 
    case('SD1gamma1');        get_ixPar = ixPar%SD1gamma1       ! soil particle density [kg/m^3]
    case('BD1gamma1');        get_ixPar = ixPar%BD1gamma1       ! soil particle density [kg/m^3]
    case('WcrFrac1gamma1');   get_ixPar = ixPar%WcrFrac1gamma1  ! Fractional soil moisture content at critical point [-] 
    case('WpwpFrac1gamma1');  get_ixPar = ixPar%WpwpFrac1gamma1 ! Fractional soil moisture content at wilting point [-]  

    case('uhshape');          get_ixPar = ixPar%uhshape           ! gamma pdf uh shape parameter [-]
    case('uhscale');          get_ixPar = ixPar%uhscale           ! gamma pdf uh scale parameter [-] 
    case('h1');               get_ixPar = ixPar%h1                ! top layer thickness [m] 
    case('h2');               get_ixPar = ixPar%h2                ! 2nd layer thickness [m] 
    case('h3');               get_ixPar = ixPar%h3                ! 3rd layer thickness [m] 
    case('h4');               get_ixPar = ixPar%h4                ! 4th layer thickness [m] 
    case('h5');               get_ixPar = ixPar%h5                ! 5th layer thickness [m] 
    case('binfilt');          get_ixPar = ixPar%binfilt           ! variable infilitration curve parameter 
    case('D1');               get_ixPar = ixPar%D1                ! fraction of Dsmax where nonlinear baseflow begins 
    case('D2');               get_ixPar = ixPar%D2                ! fraction of maximum soil moisture where nonlinear baseflow begins 
    case('D3');               get_ixPar = ixPar%D3                ! fraction of maximum soil moisture where nonlinear baseflow begins 
    case('D4');               get_ixPar = ixPar%D4                ! fraction of maximum soil moisture where nonlinear baseflow begins 
    case('expt');             get_ixPar = ixPar%expt              ! exponent in Campbell equatin for Kh
    case('ks');               get_ixPar = ixPar%ks                ! saturated hydrologic conductivity [mm/day]
    case('bbl');              get_ixPar = ixPar%bbl               ! bubbling pressure of soil [cm] 
    case('SD');               get_ixPar = ixPar%SD                ! soil particle density [kg/m^3]
    case('BD');               get_ixPar = ixPar%BD                ! soil particle density [kg/m^3]
    case('WcrFrac');          get_ixPar = ixPar%WcrFrac           ! Fractional soil moisture content at critical point [-] 
    case('WpwpFrac');         get_ixPar = ixPar%WpwpFrac          ! Fractional soil moisture content at wilting point [-]  
    case('phi');              get_ixPar = ixPar%phi               ! porosity [-]  
    case('fc');               get_ixPar = ixPar%fc                ! field capacity[-]  
    case('wp');               get_ixPar = ixPar%wp                ! wilting point [-]  
    case('psis');             get_ixPar = ixPar%psis              ! matric potential [kPa]  
    case('b');                get_ixPar = ixPar%b                 ! Fractional soil moisture content at wilting point [-]  
    case('myu');              get_ixPar = ixPar%myu               ! Fractional soil moisture content at wilting point [-]  
    case('rmin');             get_ixPar = ixPar%rmin              ! minimum stomatal resistance
    case('lai');              get_ixPar = ixPar%lai               ! monthly lai
    ! get to here if cannot find the variable
    case default;             get_ixPar = imiss
  endselect
 end function get_ixPar

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
     case default;           get_ixDataMap = imiss
   end select
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
  use var_lookup,only:ixPrpVeg                   ! indices of the named variables
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
