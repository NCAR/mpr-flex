module model_wrapper 

  use nrtype
  use public_var

  implicit none
  public :: adjust_param 
  public :: read_sim
  private

contains

subroutine adjust_param( idModel, param, err, message)
  use vic_routines, only: vic_soil_param,vic_vege_param
  implicit none
  
  ! input 
  integer(i4b),         intent(in)   :: idModel 
  real(dp),dimension(:),intent(in)   :: param        ! calibrating parameters
  ! output
  integer(i4b),         intent(out)  :: err          ! error code
  character(len=strLen),intent(out)  :: message      ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage     ! error message from downward subroutine

  ! Start procedure here
  err=0; message="adjust_param/"

  select case (idModel)
    case (1)
      ! Read/Adjust/Output soil model parameters 
      call vic_soil_param( param, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
      ! Read/Adjust/Output vege parameters 
      call vic_vege_param( param, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine adjust_param 
  

subroutine read_sim( idModel, sim, err, message)
  use vic_routines, only: read_vic_sim
  implicit none

  ! input 
  integer(i4b),         intent(in)   :: idModel 
  ! output
  real(dp),             intent(out)  :: sim(:,:)
  integer(i4b),         intent(out)  :: err          ! error code
  character(len=strLen),intent(out)  :: message      ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage     ! error message from downward subroutine

  ! Start procedure here
  err=0; message="read_sim/"

  select case (idModel)
    case (1)
      ! Read/Adjust/Output soil model parameters 
      call read_vic_sim( sim, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return

end subroutine read_sim 

end module model_wrapper 
