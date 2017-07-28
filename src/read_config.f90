module read_config 
  use nrtype
  use public_var

  implicit none

  private

  public :: read_nml

! Main configuration 
  namelist / runconfig /  opt,                     &
                          opt_method,              &
                          mpr_param_file
! MPR configuration
  namelist / mprconfig /  mpr_input_dir,           & 
                          mpr_output_dir,          &
                          soil_param_nc,           &
                          veg_param_nc,            &
                          fname_soil,              &
                          fname_veg,               &
                          fname_smapping,          &
                          fname_vmapping,          &
                          dname_overSpoly,         &
                          dname_overVpoly,         &
                          dname_hru,               &
                          sclass_table,            &
                          vclass_table,            &
                          nVclass,                 &
                          dname_spoly,             &
                          dname_slyrs,             &              
                          dname_vpoly
! calibration run specification 
  namelist / calconfig /  filelist_name,           &
                          cellfrac_name,           &
                          origparam_name,          & 
                          calibparam_name,         &
                          origvege_name,           & 
                          calivege_name,           &
                          region_info,             &
                          sim_dir,                 &
                          obs_name,                &
                          executable,              & 
                          basin_objfun_weight_file,&
                          agg,                     &
                          dt,                      &
                          sim_len,                 & 
                          start_cal,               &
                          end_cal,                 &
                          nHru,                    &
                          nbasin,                  &
                          isRoute
! Model specification 
  namelist / modelconfig / idModel,                &
                           TotNpar,                &    
                           nLyr,                   &
                           inParList
! DDS algorithm 
  namelist / DDS / rpar,        &
                   isRestart,   &
                   nseed,       & 
                   maxn,        &
                   ismax,       &
                   restrt_file, &
                   state_file
! SCE algorithm 
  namelist / SCE / maxn,        & 
                   percen,      &
                   nseed,       &
                   numcpx,      &
                   cpxstop,     &
                   state_file

contains

! --------------------------
subroutine read_nml(nmlfile, err, message)
  implicit none
  ! input 
  character(*), intent(in)  :: nmlfile
  ! output variables
  integer(i4b)              :: err
  character(len=strLen)     :: message    ! error message for downwind routine

  ! Start procedure here
  err=0; message="read_nml/"
  ! Open namelist file 
  open(UNIT=30, file=trim(nmlfile),status="old", action="read", iostat=err )
  if(err/=0)then; message=trim(message)//"Error:Open namelist"; return; endif
  ! read "runconfig" group 
  read(unit=30, NML=runconfig, iostat=err)
  if (err/=0)then; message=trim(message)//"Error:Read runconfig"; return; endif
  ! read "mprconfig" group 
  read(unit=30, NML=mprconfig, iostat=err)
  if (err/=0)then; message=trim(message)//"Error:Read mprconfig"; return; endif
  ! read "calconfig" group 
  read(unit=30, NML=calconfig, iostat=err)
  if (err/=0)then; message=trim(message)//"Error:Read calconfig"; return; endif
  ! read "modelconfig" group 
  read(unit=30, NML=modelconfig, iostat=err)
  if (err/=0)then; message=trim(message)//"Error:Read modelconfig"; return; endif
  select case (opt_method)
    case (1)
    ! read DDS group 
    read(UNIT=30, NML=DDS, iostat=err)
    if (err/=0)then; message=trim(message)//"Error:Read DDS"; return; endif
    case (2)
    ! read SCE group 
    read(UNIT=30, NML=SCE, iostat=err)
    if (err/=0)then; message=trim(message)//"Error:Read SCE"; return; endif
  end select
  close(UNIT=30)
  print *, 'Namelist file has been successfully processed'
  return
end subroutine

end module read_config
