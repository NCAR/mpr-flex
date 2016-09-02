module get_ixname

! Purpose
! Define functions to get the index of a named variable

use nrtype                  ! variable types, etc.
use public_var

implicit none

! the followings are public
public::get_ixPar     ! assign variable index to parameter 
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
   case('Ks');               get_ixPar = ixPar%ks                ! saturated hydrologic conductivity [mm/day]
   case('bbl');              get_ixPar = ixPar%bbl               ! bubbling pressure of soil [cm] 
   case('SD');               get_ixPar = ixPar%SD                ! soil particle density [kg/m^3]
   case('BD');               get_ixPar = ixPar%BD                ! soil particle density [kg/m^3]
   case('WcrFrac');          get_ixPar = ixPar%WcrFrac           ! Fractional soil moisture content at critical point [-] 
   case('WpwpFrac');         get_ixPar = ixPar%WpwpFrac          ! Fractional soil moisture content at wilting point [-]  
   case('rmin');             get_ixPar = ixPar%rmin              ! Fractional soil moisture content at wilting point [-]  
   case('lai');              get_ixPar = ixPar%lai               ! Fractional soil moisture content at wilting point [-]  
   ! get to here if cannot find the variable
   case default;         get_ixPar = imiss
  endselect
 end function get_ixPar

end module get_ixname
