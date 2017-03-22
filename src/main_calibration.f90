program main_calibration

  use nrtype 
  use public_var
  use read_config,          only: read_nml 
  use popMeta,              only: paramMaster
  use globaldata,           only: betaMaster, parSubset, betaInGamma, gammaSubset, betaOrder, parMask, parArray, nBetaNeed
  use process_meta,         only: read_calPar, get_parm_meta, param_setup
  use tf,                   only: betaDependency, betaCompOrder
  use mo_dds,               only: dds
  use mo_opt_run,           only: opt_run
  use eval_model,           only: objfn
  use mpr_routine,          only: run_mpr

  implicit none
 
  character(len=strLen)             :: nmlfile         ! namelist containing configuration
  integer(i4b)                      :: i               ! loop index for writing
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
  print*,"!-- Beta and Gamma parameters ----"
  write(*,*) (trim(adjustl(parSubset(i)%pname)),new_line('a'), i=1,size(parSubset))
  print*,"!-- Beta parameters to be estimated with MPR excluding z and h----"
  if (size(betaInGamma)/=0)then
    write(*,*) (trim(adjustl(betaInGamma(i))),new_line('a'), i=1,size(betaInGamma))
    print*,"!-- List of gamma parameters calibrated----"
    write(*,*) (trim(adjustl(gammaSubset(i)%pname)),new_line('a'), i=1,size(gammaSubset))
    ! Compute computing order of beta parameters including dependent parameters. Saved data: 'betaOrder', nBetaNeed
    call betaCompOrder (betaInGamma, ierr, cmessage); call handle_err(ierr,cmessage)
    print*,"!-- All beta parameters to be computed with MPR including dependent beta parameters ----"
    write(*,*) (trim(adjustl(betaMaster(betaOrder(i))%pname)),new_line('a'), i=1,nBetaNeed)
  else
    write(*,*) "No beta parameters estimated with MPR" 
  endif
  ! initialize parameter and mask arrays 
  call param_setup( ierr, cmessage )
  do i=1,size(parArray,1)
    print*,parArray(i,1),parMask(i)
  enddo
  ! main routine starts depending on option
  select case (opt)
    case (0)     ! just run model and output ascii of sim and obs series (parameter values use default or ones specified in restart file)
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
    case (2)     ! just perform MPR only and output parameters in netCDF
      call run_mpr( parArray(:,1), restrt_file, ierr, cmessage ); call handle_err(ierr,cmessage)
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
