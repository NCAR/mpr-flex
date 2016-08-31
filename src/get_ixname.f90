module get_ixname

! Purpose
! Define functions to get the index of a named variable

use nrtype                  ! variable types, etc.
use public_var

implicit none

! the followings are public
public::get_ixGamma     ! assign variable index to gamma parameter 
public::get_ixBeta      ! assign variable index to beta parameter 
! everything else
private

contains

! *******************************************************************************************************************
! function: get the index of the named variables for gamma parameters 
! *******************************************************************************************************************
 function get_ixGamma(varName)
  use var_lookup,only:ixGamma             ! indices of the named variables
  implicit none
  ! define dummy variables
  character(*), intent(in) :: varName     ! variable name
  integer(i4b)             :: get_ixGamma ! index of the named variable

  ! get the index of the named variables
  select case(trim(varName))
   case('ks1gamma1');        get_ixGamma = ixGamma%ks1gamma1       ! 
   case('ks1gamma2');        get_ixGamma = ixGamma%ks1gamma2       ! 
   case('ks1gamma3');        get_ixGamma = ixGamma%ks1gamma3       ! 
   case('ks2gamma1');        get_ixGamma = ixGamma%ks2gamma1       ! 
   case('ks2gamma2');        get_ixGamma = ixGamma%ks2gamma2       ! 
   case('ks2gamma3');        get_ixGamma = ixGamma%ks2gamma3       ! 
   case('phi1gamma1');       get_ixGamma = ixGamma%phi1gamma1      ! 
   case('phi1gamma2');       get_ixGamma = ixGamma%phi1gamma2      ! 
   case('phi1gamma3');       get_ixGamma = ixGamma%phi1gamma3      ! 
   case('phi2gamma1');       get_ixGamma = ixGamma%phi2gamma1      ! 
   case('phi2gamma2');       get_ixGamma = ixGamma%phi2gamma2      ! 
   case('phi2gamma3');       get_ixGamma = ixGamma%phi2gamma3      ! 
   case('phi2gamma4');       get_ixGamma = ixGamma%phi2gamma4      ! 
   case('phi2gamma5');       get_ixGamma = ixGamma%phi2gamma5      ! 
   case('phi2gamma6');       get_ixGamma = ixGamma%phi2gamma6      ! 
   case('fc1gamma1');        get_ixGamma = ixGamma%fc1gamma1       ! 
   case('wp1gamma1');        get_ixGamma = ixGamma%wp1gamma1       ! 
   case('b1gamma1');         get_ixGamma = ixGamma%b1gamma1        ! 
   case('b1gamma2');         get_ixGamma = ixGamma%b1gamma2        ! 
   case('b1gamma3');         get_ixGamma = ixGamma%b1gamma3        ! 
   case('psis1gamma1');      get_ixGamma = ixGamma%psis1gamma1     ! 
   case('psis1gamma2');      get_ixGamma = ixGamma%psis1gamma2     ! 
   case('psis1gamma3');      get_ixGamma = ixGamma%psis1gamma3     ! 
   case('myu1gamma1');       get_ixGamma = ixGamma%myu1gamma1      ! 
   case('myu1gamma2');       get_ixGamma = ixGamma%myu1gamma2      ! 
   case('z1gamma1');         get_ixGamma = ixGamma%z1gamma1        ! Layer bottom depth [m] 
   case('h1gamma1');         get_ixGamma = ixGamma%h1gamma1        ! top layer thickness [m] 
   case('h1gamma2');         get_ixGamma = ixGamma%h1gamma2        ! 2nd layer thickness [m] 
   case('binfilt1gamma1');   get_ixGamma = ixGamma%binfilt1gamma1  ! variable infilitration curve parameter 
   case('binfilt1gamma2');   get_ixGamma = ixGamma%binfilt1gamma2  ! variable infilitration curve parameter 
   case('D11gamma1');        get_ixGamma = ixGamma%D11gamma1       ! fraction of Dsmax where nonlinear baseflow begins 
   case('D21gamma1');        get_ixGamma = ixGamma%D21gamma1       ! fraction of maximum soil moisture where nonlinear baseflow begins 
   case('D31gamma1');        get_ixGamma = ixGamma%D31gamma1       ! fraction of maximum soil moisture where nonlinear baseflow begins 
   case('D41gamma1');        get_ixGamma = ixGamma%D41gamma1       ! fraction of maximum soil moisture where nonlinear baseflow begins 
   case('exp1gamma1');       get_ixGamma = ixGamma%exp1gamma1      ! exponent in Campbell equatin for Kh
   case('exp1gamma2');       get_ixGamma = ixGamma%exp1gamma2      ! exponent in Campbell equatin for Kh
   case('Ksat1gamma1');      get_ixGamma = ixGamma%ksat1gamma1     ! saturated hydrologic conductivity [mm/day]
   case('bbl1gamma1');       get_ixGamma = ixGamma%bbl1gamma1      ! bubbling pressure of soil [cm] 
   case('bbl1gamma2');       get_ixGamma = ixGamma%bbl1gamma2      ! bubbling pressure of soil [cm] 
   case('SD1gamma1');        get_ixGamma = ixGamma%SD1gamma1       ! soil particle density [kg/m^3]
   case('BD1gamma1');        get_ixGamma = ixGamma%BD1gamma1       ! soil particle density [kg/m^3]
   case('WcrFrac1gamma1');   get_ixGamma = ixGamma%WcrFrac1gamma1  ! Fractional soil moisture content at critical point [-] 
   case('WpwpFrac1gamma1');  get_ixGamma = ixGamma%WpwpFrac1gamma1 ! Fractional soil moisture content at wilting point [-]  
   ! get to here if cannot find the variable
   case default;         get_ixGamma = imiss
  endselect
 end function get_ixGamma

! *******************************************************************************************************************
! function: get the index of the named variables for beta parameter 
! *******************************************************************************************************************
 function get_ixBeta(varName)
  use var_lookup,only:ixBeta             ! indices of the named variables
  implicit none
  ! define dummy variables
  character(*), intent(in) :: varName     ! variable name
  integer(i4b)             :: get_ixBeta ! index of the named variable

  ! get the index of the named variables
  select case(trim(varName))
   case('h1');        get_ixBeta = ixBeta%h1        ! top layer thickness [m] 
   case('h2');        get_ixBeta = ixBeta%h2        ! 2nd layer thickness [m] 
   case('h3');        get_ixBeta = ixBeta%h3        ! 3rd layer thickness [m] 
   case('h4');        get_ixBeta = ixBeta%h4        ! 4th layer thickness [m] 
   case('h5');        get_ixBeta = ixBeta%h5        ! 5th layer thickness [m] 
   case('binfilt');   get_ixBeta = ixBeta%binfilt   ! variable infilitration curve parameter 
   case('D1');        get_ixBeta = ixBeta%D1        ! fraction of Dsmax where nonlinear baseflow begins 
   case('D2');        get_ixBeta = ixBeta%D2        ! fraction of maximum soil moisture where nonlinear baseflow begins 
   case('D3');        get_ixBeta = ixBeta%D3        ! fraction of maximum soil moisture where nonlinear baseflow begins 
   case('D4');        get_ixBeta = ixBeta%D4        ! fraction of maximum soil moisture where nonlinear baseflow begins 
   case('expt');      get_ixBeta = ixBeta%expt      ! exponent in Campbell equatin for Kh
   case('Ks');        get_ixBeta = ixBeta%ks        ! saturated hydrologic conductivity [mm/day]
   case('bbl');       get_ixBeta = ixBeta%bbl       ! bubbling pressure of soil [cm] 
   case('SD');        get_ixBeta = ixBeta%SD        ! soil particle density [kg/m^3]
   case('BD');        get_ixBeta = ixBeta%BD        ! soil particle density [kg/m^3]
   case('WcrFrac');   get_ixBeta = ixBeta%WcrFrac   ! Fractional soil moisture content at critical point [-] 
   case('WpwpFrac');  get_ixBeta = ixBeta%WpwpFrac  ! Fractional soil moisture content at wilting point [-]  
   case('rmin');      get_ixBeta = ixBeta%rmin      ! Fractional soil moisture content at wilting point [-]  
   case('lai');       get_ixBeta = ixBeta%lai       ! Fractional soil moisture content at wilting point [-]  
   ! get to here if cannot find the variable
   case default;         get_ixBeta = imiss
  endselect
 end function get_ixBeta

end module get_ixname
