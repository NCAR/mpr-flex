module public_var

  use nrtype 
  implicit none

  save

  ! some common variables
  integer,parameter,public        :: strLen=256          ! length of character string
  ! some constant variables
  integer(i4b),parameter, public  :: imiss=-9999 

  ! Namelist variables
  !
  ! Control variables
  character(len=strLen),public    :: filelist_name
  character(len=strLen),public    :: cellfrac_name
  character(len=strLen),public    :: origparam_name
  character(len=strLen),public    :: calibparam_name
  character(len=strLen),public    :: region_info
  character(len=strLen),public    :: outputmod_name
  character(len=strLen),public    :: outputobs_name
  character(len=strLen),public    :: executable
  character(len=strLen),public    :: spin_exe
  character(len=strLen),public    :: basin_objfun_weight_file
  integer(i4b),public             :: dt
  integer(i4b),public             :: sim_len 
  integer(i4b),public             :: start_cal
  integer(i4b),public             :: end_cal 
  integer(i4b),public             :: opt
  integer(i4b),public             :: Ncells
  integer(i4b),public             :: ntot
  integer(i4b),public             :: nbasin
  real(dp),public                 :: BasinArea
  integer(i4b),public             :: Npro
  integer(i4b),public             :: region_flag = 0
  integer(i4b),public             :: upscale_flag = 0
  real,dimension(16),public       :: initcell,endcell
  character(len=strLen),public    :: eval_length
  character(len=strLen),public    :: calpar 
  !
  !DDS 
  !input parameters
  integer(i4b),public             :: NparCal        ! number of parameters in namelist
  real(dp),public                 :: r              ! search  
  integer(i8b),public             :: seed           ! starting seed for random number generator
  integer(i8b),public             :: maxiter        ! maximum number of trials before optimization is terminated
  logical(lgt),public             :: maxit          ! maximization or minimization
  character(len=strLen),public    :: state_file     ! state file
  character(len=strLen),public    :: restrt_file    ! restart file
  logical,public                  :: isRestart      ! ues restart option? 
  !
  ! VIC parameters
  integer, parameter,public       :: TotNparVic = 54

end module public_var
