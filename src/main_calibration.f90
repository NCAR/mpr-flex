program main_calibration

  use nrtype 
  use public_var
  use read_config,          only: read_nml 
  use globaldata,           only: parSubset, betaInGamma,gammaSubset, betaNeeded
  use process_meta,         only: read_parm_master_meta, read_calPar, get_parm_meta, betaCollection, total_calParam, param_setup, check_gammaZ, check_gammaH
  use mo_dds,               only: dds
  use mo_opt_run,           only: opt_run
  use eval_model,           only: objfn
  use mpr_routine,          only: run_mpr

  implicit none
 
  character(len=strLen)             :: nmlfile         ! namelist containing configuration
  real(dp),             allocatable :: param(:,:)      ! initial value for parameter to be calibrated 
  logical(lgc),         allocatable :: parMask(:)      ! mask of parameter to be calibrated
  integer(i4b)                      :: i               ! loop index for writing
  integer(i4b)                      :: ierr            ! error code 
  character(len=strLen)             :: cmessage        ! error message from suroutine

  nmlfile='namelist.dds.local'
  ! Read configuration namelists and save variables 
  call read_nml( trim(nmlfile), ierr, cmessage ); call handle_err(ierr,cmessage)
  ! Read parameter master metadata and populate "parMaster" structure
  call read_parm_master_meta(trim(param_master_meta), ierr, cmessage); call handle_err(ierr,cmessage)
  ! Read 'CalPar' input listing metadata of beta parameters to be estimated.  Saved data structure: 'calParMeta'
  call read_calPar( trim(calpar), ierr,  cmessage ); call handle_err(ierr,cmessage)
  ! Process 'CalParMeta', save a subset of parameter meta from master. 'betaInGamma', Saved data strucutres: 'parSubset','gammaSubset'  
  call get_parm_meta(ierr,cmessage); call handle_err(ierr,cmessage)
  ! check if MPR is used, z and h gamma parameters are required
  call check_gammaZ( ierr, cmessage); call handle_err(ierr,cmessage)
  call check_gammaH( ierr, cmessage); call handle_err(ierr,cmessage)
  ! Identify all the beta parameters to be estimated - 'betaNeeded' incluging betaInGamma
  call betaCollection( ierr, cmessage); call handle_err(ierr,cmessage)
  ! print out list of gamma/beta parameters
  print*,"!-- Beta and Gamma parameters ----"
  write(*,*) (trim(adjustl(parSubset(i)%pname)),new_line('a'), i=1,size(parSubset))
  print*,"!-- Beta parameters to be estimated with MPR ----"
  write(*,*) (trim(adjustl(betaInGamma(i))),new_line('a'), i=1,size(betaInGamma))
  print*,"!-- List of gamma parameters calibrated----"
  write(*,*) (trim(adjustl(gammaSubset(i)%pname)),new_line('a'), i=1,size(gammaSubset))
  print*,"!-- All beta parameters to be computed with MPR including dependent beta parameters ----"
  write(*,*) (trim(adjustl(betaNeeded(i))),new_line('a'), i=1,size(betaNeeded))
  ! initialize parameter and mask arrays 
  call total_calParam()
  allocate(param(nParCalSum,3))
  allocate(parMask(nParCalSum))
  call param_setup(param, parMask)
  ! main routine starts depending on option
  select case (opt)
    case (0)     ! just run model and output ascii of sim and obs series (parameter values use default or ones specified in restart file)
      call opt_run(objfn, restrt_file)
    case (1)     ! perform calibration with DDS
      call dds(objfn,                   & ! function to get object function
               param(:,1),              & ! initial parameter values
               param(:,2:3),            & ! lower and upper bounds of each parameters
               isRestart,               & ! .true. - use restart file to start, otherwise from beginning 
               restrt_file,             & ! restart file to write the most recent param values, etc
               r=rpar,                  & ! perturbation window (optional)
               mask=parMask,            & ! mask (optional)
               seed=nseed,              & ! seed for random number
               maxiter=maxn,            & ! maximum iteration
               maxit=isMax,             & ! minimzation (0) or maximization (1)
               tmp_file=state_file)       !
    case (2)     ! just perform MPR only and output parameters in netCDF
      call run_mpr( param(:,1), restrt_file, ierr, cmessage ); call handle_err(ierr,cmessage)
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
