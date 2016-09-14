module data_type

! Here define custum data type    

  use nrtype

  implicit none

! ***********************************************************************************************************
! Define data structure of master parameter (both gamma and beta) metadata
! ***********************************************************************************************************
type,public  :: par_meta
  character(len=64)   :: pname=''        ! parameter name
  real(dp)            :: val  =-999.0_dp ! default bound 
  real(dp)            :: lwr  =-999.0_dp ! lower and upper bounds
  real(dp)            :: upr  =-999.0_dp ! lower and upper bounds
  character(len=64)   :: beta ='beta'    ! name of parent beta parameter - if parameter is beta parameter, use "beta"
  logical(lgt)        :: flag =.False.   ! flag to calibrate or not 
  character(len=64)   :: hups ='wamean'  ! name of parent beta parameter - if parameter is beta parameter, use "beta"
  character(len=64)   :: vups ='wamean'  ! name of parent beta parameter - if parameter is beta parameter, use "beta"
endtype par_meta

! extended parameter meta data for selected set 
type,extends(par_meta), public  :: cpar_meta
  integer(i4b)        :: ixMaster=-999   ! idex of master parameter list
endtype cpar_meta

! ***********************************************************************************************************
! Define data structure of beta parameter meta data 
! ***********************************************************************************************************
type,public  :: beta_meta
  character(len=64)                     :: pname=''             ! beta parameter name
  character(len=128)                    :: desc=''              ! description
  character(len=64)                     :: units=''             ! units
  character(len=32)                     :: dims=''              ! data dimension (scaler, vector, 2D, 3D)
  logical(lgt)                          :: v_agg=.FALSE.        ! flag to perform vertical scaling 
  character(len=32)                     :: vaggmethod           ! vertical upscaling operator to be used 
  logical(lgt)                          :: h_agg=.FALSE.        ! flag to perform horizontal scaling 
  character(len=32)                     :: haggmethod           ! horizontal upscaling operator to be used 
endtype beta_meta

! ***********************************************************************************************************
! Define data structure of variable metadata - soil propeties, topography, vege propeties, model hru propeties 
! ***********************************************************************************************************
type,public :: var_meta
  character(len=64)                     :: varname=''               ! name
  character(len=128)                    :: vardesc=''               ! description
  character(len=64)                     :: varunit=''               ! units
  character(len=32)                     :: vardims=''               ! dimension (scaler, vector, 2D, 3D)
  character(len=32)                     :: vartype=''               ! type (integer, double)
endtype var_meta

! *****
! data structure to hold polygon data only integer type
! ********************************************
! ipolydata(:)%var(:)%data(:)
! 1st level - horizontal
! 2nd level - variable index 
! 3rd level - data
type,public :: ivar 
  integer(i4b),allocatable       :: dat(:)
endtype ivar 
type,public :: ipolydata 
  type(ivar),allocatable         :: var(:)
endtype ipolydata 

! *****
! data structure to hold polygon data only real type
! ********************************************
! ipolydata(:)%var(:)%data(:)
! 1st level - horizontal
! 2nd level - variable index 
! 3rd level - data
type,public :: dvar 
  real(dp),allocatable      :: dat(:)
endtype dvar 
type,public :: dpolydata 
  type(dvar),allocatable    :: var(:)
endtype dpolydata 

! *****
! data structure to hold integer type data (both vector and 2D) 
! ********************************************
! poly(:)%layer(:)%weight(:)
!                 %ixSubLyr(:)
! 1st level - horizontal
! 2nd level - model layer 
! 3rd level - weight and index of sublayer 
type,public :: mapping 
  real(dp),    allocatable          :: weight(:)
  integer(i4b),allocatable          :: ixSubLyr(:)
endtype mapping 
type,public :: poly 
  type(mapping),allocatable    :: layer(:)
endtype poly

! *****
! Other data type (need to clean up)
! ********************************************
! Define derived types to hold data values in vector (soil properties or model parameters) given their indices 
! use "layer" for 2nd layer name
! ** double precision type 
 type,public :: lyr_d
   real(dp),allocatable              :: layer(:)
 endtype lyr_d
 ! ** integer type 
 type,public :: lyr_i
   integer(i4b),allocatable          :: layer(:)
 endtype lyr_i
! use "var" for 2nd layer name
! ** double precision type 
 type, public :: var_d
  real(dp),allocatable               :: var(:)
 endtype var_d
 ! ** integer type 
 type, public :: var_i
  integer(i4b),allocatable           :: var(:)
 endtype var_i
! use "dat" for 2nd layer name
! ** double precision type
 type,public :: dat_d1d
   real(dp),allocatable              :: dat(:)
 endtype dat_d1d
 ! ** integer type
 type,public :: dat_i1d
   integer(i4b),allocatable          :: dat(:)
 endtype dat_i1d

! Define derived types to hold data values in 2D array (soil properties or model parameters) given their indices 
! ** double precision type
 type,public :: dat_d2d
   real(dp),allocatable         :: dat(:,:)
 endtype dat_d2d
 ! ** integer type
 type,public :: dat_i2d
   integer(i4b),allocatable     :: dat(:,:)
 endtype dat_i2d

! Define derived types to hold data values in vector and 2D array (soil properties or model parameters) given their indices 
! ** double precision type
 type,public :: dat_d12d
   real(dp),allocatable         :: dat1d(:)
   real(dp),allocatable         :: dat2d(:,:)
 endtype dat_d12d
 ! ** integer type
 type,public :: dat_i12d
   integer(i4b),allocatable     :: dat1d(:)
   integer(i4b),allocatable     :: dat2d(:,:)
 endtype dat_i12d

! more layered data type
! var_dlength(:)%var(:)%layer(:)
! ** double precision type
 type,public :: var_dlength
   type(lyr_d),allocatable     :: var(:)
 endtype var_dlength

 ! ** integer type
 type,public :: var_ilength
   type(lyr_i),allocatable     :: var(:)
 endtype var_ilength

! data type containing a name and vector variable (double precision)
 type namedvar
  character(len=256)           :: varName
  real(dp),allocatable         :: varData(:)
 endtype namedvar

! data type containing a name and vector variable (integer)
 type nameivar
   character(len=256)          :: varName
   integer(i4b),allocatable    :: varData(:)
 endtype nameivar

! data type containing a name and vector variable (double precision)
 type namedvar2
  character(len=256)           :: varName
  real(dp),allocatable         :: varData(:,:)
 endtype namedvar2

! data type containing a name and vector variable (integer)
 type nameivar2
   character(len=256)          :: varName
   integer(i4b),allocatable    :: varData(:,:)
 endtype nameivar2

! data type containing a name and vector and 2D array variable for both integer and double precision
 type namevar
   character(len=256)          :: varName
   real(dp),    allocatable    :: dvar1(:)
   integer(i4b),allocatable    :: ivar1(:)
   real(dp),    allocatable    :: dvar2(:,:)
   integer(i4b),allocatable    :: ivar2(:,:)
 endtype namevar

 type mapvar
   type(namevar),allocatable         :: var(:)
 endtype mapvar

end module data_type
