program main_calibration

  use nrtype 
  use public_var
  use read_config,     only: read_nml 
  use popMeta,         only: paramMaster
  use globaldata,      only: betaInGamma, parMask, parArray
  use process_meta,    only: read_calPar, get_parm_meta, param_setup, print_config
  use tf,              only: betaDependency, betaCompOrder
  use mo_dds,          only: dds
  use mo_opt_run,      only: opt_run
  use eval_model,      only: objfn
  use mpr_routine,     only: run_mpr
  use read_soildata,   only: check_polyID 

  implicit none
 
  character(len=strLen)             :: nmlfile         ! namelist containing configuration
  integer(i4b)                      :: ierr            ! error code 
  character(len=strLen)             :: cmessage        ! error message from suroutine

  nmlfile='namelist.dds.local'
  ! Read configuration namelists and save variables 
  call read_nml( trim(nmlfile), ierr, cmessage ); call handle_err(ierr,cmessage)
  ! Populate master parameter meta.  Saved data: betaMaster, gammaMaster 
  call paramMaster( ierr, cmessage ); call handle_err(ierr,cmessage)
  ! Read 'CalPar' input listing metadata of beta parameters to be estimated. Saved data: 'calParMeta'
  call read_calPar( trim(calpar), ierr,  cmessage ); call handle_err(ierr,cmessage)
  ! Process 'CalParMeta' along with master parameter meta data. 
  ! Saved data: 'parSubset','gammaSubset', 'betaInGamma', betaCalScale, nBetaDirCal, nBetaGammaCal, nGammaCal, nnParCalSum 
  call get_parm_meta(ierr,cmessage); call handle_err(ierr,cmessage)
  ! Compute beta parameter dependency. Saved data: beta 
  call betaDependency (ierr, cmessage); call handle_err(ierr,cmessage)
  if (size(betaInGamma)/=0)then
    ! Compute computing order of beta parameters including dependent parameters. Saved data: 'betaOrder', nBetaNeed
    call betaCompOrder (betaInGamma, ierr, cmessage); call handle_err(ierr,cmessage)
    call check_polyID(trim(mpr_input_dir)//trim(fname_soil), dname_spoly , ierr, cmessage); call handle_err(ierr, cmessage)
  endif
  ! initialize parameter and mask arrays 
  call param_setup( ierr, cmessage ); call handle_err(ierr,cmessage)
  ! Print out calibration configuration
  call print_config()
  ! main routine starts depending on option
  select case (opt)
    case (2)     ! just run model and output ascii of sim and obs series (parameter values use default or ones specified in restart file)
      call opt_run(objfn, restrt_file)
    case (1)     ! perform calibration with DDS
      call dds(objfn,                   & ! function to get object function
               parArray(:,1),           & ! initial parameter values
               parArray(:,2:3),         & ! lower and upper bounds of each parameters
               isRestart,               & ! .true. - use restart file to start, otherwise from beginning 
               restrt_file,             & ! restart file to write the most recent param values, etc
               r=rpar,                  & ! perturbation window (optional)
               mask=parMask,            & ! mask (optional)
               seed=nseed,              & ! seed for random number
               maxiter=maxn,            & ! maximum iteration
               maxit=isMax,             & ! minimzation (0) or maximization (1)
               tmp_file=state_file)       !
    case (3)     ! just perform MPR only and output parameters in netCDF
      call run_mpr( parArray(:,1), mpr_param_file , ierr, cmessage ); call handle_err(ierr,cmessage)
    case default
      print*, 'integer to specify optimization scheme is not valid' 
  end select 
  stop 

contains

  subroutine handle_err(err,message)
    ! handle error codes
    implicit none
    integer(i4b),intent(in)::err             ! error code
    character(*),intent(in)::message         ! error message
    if(err/=0)then
     print*,'FATAL ERROR: '//trim(message)
     stop
    endif
  end subroutine

end program main_calibration
