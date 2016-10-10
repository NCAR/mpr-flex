module model_wrapper 
  use nrtype
  use data_type                                        ! Including custum data structure definition
  use public_var

  implicit none
  
  private

  public :: adjust_param 
  public :: replace_param
  public :: read_sim
  public :: read_soil_param 
  public :: write_soil_param 
  public :: read_soil_lyr
  public :: read_hru_id

contains

subroutine adjust_param( idModel, param, multiplier, adjParam, err, message)
  use vic_routines, only: adj_soil_param_vic,adj_vege_param_vic
  implicit none
  ! input 
  integer(i4b),         intent(in)   :: idModel 
  real(dp),             intent(in)   :: param(:,:)    ! original soil parameters 
  real(dp),             intent(in)   :: multiplier(:) ! calibrating parameters
  ! output
  real(dp),             intent(out)  :: adjParam(:,:) ! adjusted soil parameter
  integer(i4b),         intent(out)  :: err           ! error code
  character(len=strLen),intent(out)  :: message       ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage      ! error message from downward subroutine

  ! Start procedure here
  err=0; message="adjust_param/"
  select case (idModel)
    case (1)
      ! Read/Adjust/Output soil model parameters 
      call adj_soil_param_vic( param, multiplier, adjParam, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
      ! Read/Adjust/Output vege parameters 
      call adj_vege_param_vic( multiplier, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine

subroutine replace_param( idModel, param, hModel, parMxyMz, adjParam, err, message)
  use vic_routines, only: replace_soil_param_vic
  implicit none
  
  ! input 
  integer(i4b),         intent(in)   :: idModel       ! model id
  real(dp),             intent(in)   :: param(:,:)    ! original soil parameters 
  real(dp),             intent(in)   :: hModel(:,:)   ! Model layer thickness for model hru
  type(namedvar2),      intent(in)   :: parMxyMz(:)   ! model soil parameter 
  ! output
  real(dp),             intent(out)  :: adjParam(:,:) ! adjusted soil parameter
  integer(i4b),         intent(out)  :: err           ! error code
  character(len=strLen),intent(out)  :: message       ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage      ! error message from downward subroutine

  ! Start procedure here
  err=0; message="replace_param/"
  select case (idModel)
    case (1)
      ! replace soil model parameters 
      call replace_soil_param_vic( param, hModel, parMxyMz, adjParam,  err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine
  
subroutine read_soil_param(idModel, param, err, message)
  use vic_routines, only: read_soil_param_vic
  implicit none
  ! input 
  integer(i4b),         intent(in)   :: idModel 
  ! output
  real(dp),             intent(out)  :: param(:,:)   !  
  integer(i4b),         intent(out)  :: err          ! error code
  character(len=strLen),intent(out)  :: message      ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage     ! error message from downward subroutine

  ! Start procedure here
  err=0; message="read_soil_param/"
  select case (idModel)
    case (1)
      call read_soil_param_vic( param, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine

subroutine write_soil_param(idModel, hruid, param, err, message)
  use vic_routines, only: write_soil_param_vic
  implicit none
  ! input 
  integer(i4b),         intent(in)   :: idModel 
  integer(i4b),         intent(in)   :: hruid(:)     ! hru ID
  real(dp),             intent(in)   :: param(:,:)   !  
  ! output
  integer(i4b),         intent(out)  :: err          ! error code
  character(len=strLen),intent(out)  :: message      ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage     ! error message from downward subroutine

  ! Start procedure here
  err=0; message="write_soil_param/"
  select case (idModel)
    case (1)
      call write_soil_param_vic( hruid, param, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine

subroutine read_hru_id(idModel, hruid, err, message)
  use vic_routines, only: vic_hru_id
  implicit none
  ! input 
  integer(i4b),         intent(in)   :: idModel 
  ! output
  integer(i4b),         intent(out)  :: hruid(:)     !  
  integer(i4b),         intent(out)  :: err          ! error code
  character(len=strLen),intent(out)  :: message      ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage     ! error message from downward subroutine

  err=0; message="read_hru_id/"
  select case (idModel)
    case (1)
      call vic_hru_id( hruid, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine

subroutine read_soil_lyr(idModel, hlyr, err, message)
  use vic_routines, only: vic_soil_layer
  implicit none
  ! input 
  integer(i4b),         intent(in)   :: idModel 
  ! output
  real(dp),             intent(out)  :: hlyr(:,:)    !  
  integer(i4b),         intent(out)  :: err          ! error code
  character(len=strLen),intent(out)  :: message      ! error message
  ! LOCAL VARIABLES
  character(len=strLen)              :: cmessage     ! error message from downward subroutine

  err=0; message="read_soil_lyr/"
  select case (idModel)
    case (1)
      call vic_soil_layer( hlyr, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine

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

  err=0; message="read_sim/"
  select case (idModel)
    case (1)
      call read_vic_sim( sim, err, cmessage)
      if (err/=0)then; message=message//cmessage; return; endif
    case default
      err=10; message=message//"model is not implemented"; return
  end select  
  return
end subroutine

end module model_wrapper 
