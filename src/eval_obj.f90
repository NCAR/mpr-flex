module eval_obj

  use nrtype 
  use public_var
  use vic_subroutines_parallel, only: eval_objfn

  implicit none

  private

  public functn

contains
 
  function functn(pin)
    implicit none
    ! input variables
    real(dp),dimension(NparCal),intent(in) :: pin     ! parameter set 
    ! local variables 
    integer(i4b)                           :: err          ! error code
    character(len=strLen)                  :: message      ! error message
    real(dp)                               :: objfnc
    real(dp)                               :: functn
    
    ! initialize error control
    err=0; message='functn/'
    ! Evaluate objective function given parameter set
    call eval_objfn(pin,objfnc,err,message)
    if (err/=0) stop message 
    ! save objective function value
    functn = objfnc 
  end function functn 

end module eval_obj
