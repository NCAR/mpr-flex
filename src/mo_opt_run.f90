module mo_opt_run

  use nrtype,    only: i4b, i8b, dp
  use public_var

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: opt_run    ! run model optimized parameter 

CONTAINS
  ! ------------------------------------------------------------------
  !     NAME
  !>        \opt_run

  !     PURPOSE
  !>        \run optimized run using restart file

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: obj_func(p)"             Function on which to search the minimum
  !>        \param[in] "real(dp) :: pval(:)"                 inital value of decision variables
  !>        \logical,            :: restart                  logical to read state file and initialize, .false. -> start from begining

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     RETURN
  !>        None 

  !     RESTRICTIONS
  !         None.

  subroutine opt_run(obj_func, restartFile)

    implicit none

    INTERFACE
       function obj_func(pp)
         use nrtype, only: dp
         implicit none
         real(dp), dimension(:), intent(in) :: pp
         real(dp) :: obj_func
       end function obj_func
    END INTERFACE
    character(len=strLen),      intent(in)  :: restartFile ! name of restart file including iteration, the most recent parameter values 

    ! Local variables
    real(dp),   dimension(:), allocatable   :: pval        ! inital value of decision (parameter) variables
    integer(i8b)                            :: i           ! loop index 
    integer(i8b)                            :: iDummy      ! dummy interger: fist line of restart file starting index of objective function evaluation 
    real(dp)                                :: rDummy      ! dummy real: intermediate results for objective function values logical                                 
    logical                                 :: isExistFile ! logical to check if the file exist or not
    
    allocate ( pval(NparCal) )
    ! restart option
    print*, 'read restart file'
    inquire(file=trim(adjustl(restartFile)), exist=isExistFile)
    if ( isExistFile ) then !  if state file exists, update iStart and pval, otherwise iteration start with very beginning
      open(unit=70,file=trim(adjustl(restartFile)), action='read', status = 'unknown')
      read(70,*) iDummy
      read(70,*) (pval(i),i=1,NparCal)    
      close(70)
    endif
    ! Evaluate initial solution and return objective function value
    ! and Initialise the other variables (e.g. of_best)
    rDummy =  obj_func(pval)

  end subroutine

end module mo_opt_run
