module public_var

  use nrtype 

  implicit none

  save

  ! some common constant variables
  integer,     parameter,public   :: strLen=256             ! length of character string
  integer(i4b),parameter,public   :: nMonth=12              ! number of months 
  integer(i4b),parameter,public   :: imiss=-999             ! missing value for integer value
  real(dp),    parameter,public   :: dmiss=-999.0_dp        ! missing value for floating value
  real(dp),    parameter,public   :: verySmall=tiny(1.0_dp) ! a very small number 
  real(dp),    parameter,public   :: valMin=1.e-10_dp       ! minimum value for positive value

  ! Namelist variables
  ! &runconfig 
  integer(i4b),         public    :: opt
  integer(i4b),         public    :: opt_method
  character(len=strLen),public    :: mpr_param_file
  ! mprconfig
  character(len=strLen),public    :: mpr_input_dir 
  character(len=strLen),public    :: mpr_output_dir
  character(len=strLen),public    :: soil_param_nc 
  character(len=strLen),public    :: veg_param_nc 
  character(len=strLen),public    :: fname_soil
  character(len=strLen),public    :: fname_veg
  character(len=strLen),public    :: fname_smapping
  character(len=strLen),public    :: fname_vmapping
  character(len=strLen),public    :: dname_overSpoly
  character(len=strLen),public    :: dname_overVpoly
  character(len=strLen),public    :: dname_hru
  character(len=strLen),public    :: sclass_table
  integer(i4b)                    :: nSclass
  character(len=strLen),public    :: vclass_table
  integer(i4b)                    :: nVclass
  character(len=strLen),public    :: dname_spoly 
  character(len=strLen),public    :: dname_slyrs
  character(len=strLen),public    :: dname_vpoly 
  ! calconfig 
  character(len=strLen),public    :: filelist_name
  character(len=strLen),public    :: cellfrac_name
  character(len=strLen),public    :: origparam_name
  character(len=strLen),public    :: calibparam_name
  character(len=strLen),public    :: origvege_name
  character(len=strLen),public    :: calivege_name
  character(len=strLen),public    :: region_info
  character(len=strLen),public    :: sim_dir
  character(len=strLen),public    :: obs_name
  character(len=strLen),public    :: executable
  character(len=strLen),public    :: basin_objfun_weight_file
  integer(i4b),         public    :: agg
  integer(i4b),         public    :: dt
  integer(i4b),         public    :: sim_len 
  integer(i4b),         public    :: start_cal
  integer(i4b),         public    :: end_cal 
  integer(i4b),         public    :: nHru                       ! sum of hrus where model run
  integer(i4b),         public    :: nbasin                     ! number of basin to be jointly calibrated 
  logical(lgc),         public    :: isRoute                    ! T if sim runoff is routed w/ gamma UH implemented in mpr-flex, F if sim runoff is routed outside
  ! modelconfig
  integer(i4b),         public    :: idModel 
  integer(i4b),         public    :: TotNpar  != 54(VIC)
  integer(i4b),         public    :: nLyr     != 3 (VIC)
  character(len=strLen),public    :: inParList 
  ! optimization routine
  ! common
  integer(i8b),         public    :: maxn                       ! maximum number of trials before optimization is terminated
  integer(i8b),         public    :: nseed                      ! starting seed for random number generator
  character(len=strLen),public    :: state_file                 ! state file
  character(len=strLen),public    :: restrt_file                ! restart file
  logical(lgc),         public    :: isRestart                  ! ues restart option? 
  ! DDS 
  real(dp),             public    :: rpar                       ! search radius
  logical(lgc),         public    :: isMax                      ! maximization or minimization
  ! SCE 
  real(dp),             public    :: percen 
  integer(i4b),         public    :: numcpx                     ! 
  integer(i4b),         public    :: cpxstop                    ! 

end module public_var
