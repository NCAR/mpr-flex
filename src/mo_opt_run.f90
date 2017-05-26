module mo_opt_run

  use nrtype,    only: i4b, i8b, dp
  use public_var

  implicit none

  private

  public :: opt_run    ! run model optimized parameter 

contains

subroutine opt_run( restartFile, err, message )
  use globalData, only:nCalParSum ! meta for beta parameter listed in 'inParList' input
  use eval_model, only:out_sim ! meta for beta parameter listed in 'inParList' input
  implicit none
  ! input
  character(len=strLen),      intent(in)  :: restartFile ! name of restart file including iteration, the most recent parameter values 
  !output variables
  integer(i4b),               intent(out) :: err         ! error id 
  character(len=strLen),      intent(out) :: message     ! error message
  ! Local variables
  real(dp),dimension(nCalParSum)          :: pval        ! inital value of decision (parameter) variables
  integer(i8b)                            :: i           ! loop index 
  integer(i8b)                            :: iDummy      ! dummy interger: fist line of restart file starting index of objective function evaluation 
  logical(lgc)                            :: isExistFile ! logical to check if the file exist or not
  character(len=strLen)                   :: cmessage    ! error message
  
  err=0; message='opt_run/' ! to initialize error control
  ! restart option
  print*, 'read restart file'
  inquire(file=trim(adjustl(restartFile)), exist=isExistFile)
  if ( isExistFile ) then !  if state file exists, update iStart and pval, otherwise iteration start with very beginning
    open(unit=70,file=trim(adjustl(restartFile)), action='read', status = 'unknown')
    read(70,*) iDummy
    read(70,*) (pval(i),i=1,nCalParSum)    
    close(70)
  else
    stop 'no restart file:do optimization first'
  endif
  call out_sim(pval, err, cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  return
end subroutine

end module mo_opt_run
