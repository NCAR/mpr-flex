module mo_nml 

  use nrtype
  use public_var

  implicit none

  public :: read_nml
  private

! time control namelist for running models
  namelist / INIT_CONTROL /  filelist_name,           &
                             cellfrac_name,           &
                             origparam_name,          & 
                             calibparam_name,         &
                             region_info,             &
                             sim_dir,                 &
                             obs_name,                &
                             executable,              & 
                             basin_objfun_weight_file,&
                             dt,                      &
                             sim_len,                 & 
                             start_cal,               &
                             end_cal,                 &
                             opt,                     & 
                             Ncells,                  &
                             ntot,                    &
                             nbasin,                  &
                             upscale_flag,            &
                             Npro,                    & 
                             initcell,                & 
                             endcell,                 &
                             eval_length,             &
                             calpar

! DDS algorithm 
  namelist / DDS / NparCal,     & 
                   r,           &
                   isRestart,   &
                   seed,        & 
                   maxiter,     &
                   maxit,       &
                   restrt_file, &
                   state_file

contains

! --------------------------
subroutine read_nml(nmlfile, err, message)
  
  implicit none
  ! input 
  character(*), intent(in)  :: nmlfile
  ! output variables
  integer                   :: err
  character(len=256)        :: message    ! error message for downwind routine

  ! Start procedure here
  err=0; message="read_nml/"
  ! Open namelist file 
  open(UNIT=30, file=trim(nmlfile),status="old", action="read", iostat=err )
  if(err/=0)then; message=trim(message)//"Error:Open namelist"; return; endif
  ! read INIT_CONTROL group 
  read(unit=30, NML=INIT_CONTROL, iostat=err)
  if (err/=0)then; message=trim(message)//"Error:Read INIT_CONTROL"; return; endif
  ! read DDS group 
  read(UNIT=30, NML=DDS, iostat=err)
  if (err/=0)then; message=trim(message)//"Error:Read DDS"; return; endif

  close(UNIT=30)

  print *, 'Namelist file has been successfully processed'

end subroutine read_nml

end module mo_nml 
