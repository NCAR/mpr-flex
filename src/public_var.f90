module public_var

  use nrtype 

  implicit none

  save

  ! some common variables
  integer,     parameter,public  :: strLen=256             ! length of character string
  ! some constant variables
  integer(i4b),parameter,public  :: imiss=-999             ! missing value for integer value
  real(dp),    parameter,public  :: dmiss=-999.0_dp        ! missing value for floating value
  real(dp),    parameter,public  :: verySmall=tiny(1.0_dp) ! a very small number 
  real(dp),    parameter,public  :: valMin=1.e-10_dp       ! minimum value for positive value

  ! Namelist variables
  ! runconfig 
  integer(i4b),public             :: opt
  ! mprconfig
  character(len=strLen),public    :: mpr_input_dir 
  character(len=strLen),public    :: mpr_output_dir
  character(len=strLen),public    :: param_nc 
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
  character(len=strLen),public    :: objfntype 
  integer(i4b),public             :: agg
  integer(i4b),public             :: dt
  integer(i4b),public             :: sim_len 
  integer(i4b),public             :: start_cal
  integer(i4b),public             :: end_cal 
  integer(i4b),public             :: nHru
  integer(i4b),public             :: nbasin
  ! modelconfig
  integer(i4b),public             :: idModel 
  integer(i4b),public             :: TotNpar  != 54(VIC)
  integer(i4b),public             :: nLyr     != 3 (VIC)
  character(len=strLen),public    :: calpar 
  !
  !DDS 
  !input parameters
  real(dp),public                 :: rpar           ! search  
  integer(i8b),public             :: nseed          ! starting seed for random number generator
  integer(i8b),public             :: maxn           ! maximum number of trials before optimization is terminated
  logical(lgc),public             :: isMax          ! maximization or minimization
  character(len=strLen),public    :: state_file     ! state file
  character(len=strLen),public    :: restrt_file    ! restart file
  logical,public                  :: isRestart      ! ues restart option? 

  integer(i8b)                    :: nParCalSum     ! Total number of calibratin parameters including per layer parameters 

end module public_var
