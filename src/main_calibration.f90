program main_calibration

  use nrtype 
  use public_var
  use mo_nml,                   only: read_nml 
  use paramMaster,              only: popMeta
  use subset_meta,              only: get_parm_meta, param_setup
  use mo_dds,                   only: dds
  use mo_opt_run,               only: opt_run
  use vic_routines,             only: vic_soil_param 
  use eval_model,               only: objfn

  implicit none
 
  character(len=strLen)            :: nmlfile         ! namelist containing configuration
  real(dp),allocatable             :: param(:,:)      ! initial value for parameter to be calibrated 
  logical,allocatable              :: mask(:)         ! mask of parameter to be calibrated
  integer(i4b)                     :: ierr            ! error code 
  character(len=strLen)            :: cmessage        ! error message of downwind routine

  ! read calibration namelists and save variables 
  nmlfile='namelist.dds.local'
  call read_nml( trim(nmlfile), ierr, cmessage ); call handle_err(ierr,cmessage)
  ! Populate master parameter meta  
  call popMeta( ierr, cmessage ); call handle_err(ierr,cmessage)
  ! read calibrating parameter list 
  call get_parm_meta( trim(calpar), ierr,cmessage); call handle_err(ierr,cmessage)
  ! initialize parameter and mask arrays 
  allocate(param(nParCal,3))
  allocate(mask(nParCal))
  call param_setup(param, mask)
  ! optimization start
  select case (opt)
    case (1)     ! DDS
      call dds(objfn,                   & ! function to get object function
               param(:,1),              & ! initial parameter values
               param(:,2:3),            & ! lower and upper bounds of each parameters
               isRestart,               & ! .true. - use restart file to start, otherwise from beginning 
               restrt_file,             & ! restart file to write the most recent param values, etc
               r=r,                     & ! perturbation window (optional)
               mask=mask,               & ! mask (optional)
               seed=seed,               & ! seed for random number
               maxiter=maxiter,         & ! maximum iteration
               maxit=maxit,             & ! minimzation (0) or maximization (1)
               tmp_file=state_file)       !
    case (2)     ! just output ascii of sim and obs series 
      call opt_run(objfn, restrt_file) 
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
  end subroutine handle_err

end program main_calibration
