module model_wrapper 
  use nrtype
  use data_type                                        ! Including custum data structure definition
  use public_var

  implicit none
  
  private

  public :: adjust_param 
  public :: replace_param
  public :: read_sim
  public :: read_simRouted
  public :: read_soil_param 
  public :: write_soil_param 
  public :: read_soil_lyr
  public :: read_hru_id
contains

  subroutine adjust_param( idModel, param, multiplier, adjParam, err, message)
    use vic_routines, only: adj_soil_param_vic,adj_vege_param_vic
    use sac_routines, only: adj_soil_param_sac,adj_snow_param_sac
    implicit none
    ! input 
    integer(i4b),         intent(in)   :: idModel 
    real(dp),             intent(in)   :: param(:,:)    ! original soil parameters 
    type(var_d),          intent(in)   :: multiplier(:) ! calibrating parameters
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
        if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
        ! Read/Adjust/Output vege parameters 
        call adj_vege_param_vic( multiplier, err, cmessage)
        if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case (2)
        call adj_soil_param_sac( param, multiplier, adjParam, err, cmessage)
        if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
        call adj_snow_param_sac( multiplier, err, cmessage)
        if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case default; err=10; message=trim(message)//"model is not implemented"; return
    end select  
    return
  end subroutine
  
  subroutine replace_param( idModel, param, hModel, parMxyMz, adjParam, err, message)
    use vic_routines, only: replace_soil_param_vic
    use sac_routines, only: replace_soil_param_sac
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
  
    err=0; message="replace_param/"
    select case (idModel)
      case (1); call replace_soil_param_vic( param, hModel, parMxyMz, adjParam,  err, cmessage)
      case (2); call replace_soil_param_sac( param, parMxyMz, adjParam,  err, cmessage)
      case default; err=10; message=trim(message)//"model is not implemented"; return
    end select  
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
    return
  end subroutine
    
  subroutine read_soil_param(idModel, param, err, message)
    use vic_routines, only: read_soil_param_vic
    use sac_routines, only: read_soil_param_sac
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
      case (1); call read_soil_param_vic( param, err, cmessage)
      case (2); call read_soil_param_sac( param, err, cmessage)
      case default; err=10; message=trim(message)//"model is not implemented"; return
    end select  
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
    return
  end subroutine
  
  subroutine write_soil_param(idModel, hruid, param, err, message)
    use vic_routines, only: write_soil_param_vic
    use sac_routines, only: write_soil_param_sac
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
      case (1); call write_soil_param_vic( hruid, param, err, cmessage)
      case (2); call write_soil_param_sac( hruid, param, err, cmessage)
      case default; err=10; message=trim(message)//"model is not implemented"; return
    end select  
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
    return
  end subroutine
  
  subroutine read_hru_id(idModel, hruid, err, message)
    use vic_routines, only: vic_hru_id
    use sac_routines, only: sac_hru_id
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
      case (0); call generic_hru_id( hruid, err, cmessage) !if not specific model, read hru id list from netcdf
      case (1); call vic_hru_id( hruid, err, cmessage)
      case (2); call sac_hru_id( hruid, err, cmessage)
      case default
        err=10; message=trim(message)//"model is not implemented"; return
    end select  
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
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
      case (1); call vic_soil_layer( hlyr, err, cmessage)
      case default
        err=10; message=trim(message)//"model is not implemented"; return
    end select  
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
    return
  end subroutine
  
  subroutine read_sim( idModel, sim, err, message)
    use vic_routines, only: read_vic_sim
    use sac_routines, only: read_sac_sim
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
      case (1); call read_vic_sim( sim, err, cmessage)
      case (2); call read_sac_sim( sim, err, cmessage)
      case default
        err=10; message=trim(message)//"model is not implemented"; return
    end select  
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
    return
  end subroutine

  subroutine read_simRouted( sim, err, message)
    implicit none
    ! input 
    ! output
    real(dp),             intent(out)  :: sim(:,:)
    integer(i4b),         intent(out)  :: err          ! error code
    character(len=strLen),intent(out)  :: message      ! error message
    ! LOCAL VARIABLES
    character(len=strLen)              :: cmessage     ! error message from downward subroutine
    err=0; message="read_simRouted/"
    call read_lohmann_sim( sim, err, cmessage)
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
    return

    contains

    subroutine read_lohmann_sim( sim, err, message )
      !output variables
      real(dp),              intent(out) :: sim(:,:)
      integer(i4b),          intent(out) :: err            ! error code
      character(*),          intent(out) :: message        ! error message
      !Local variables
      real(dp)                           :: auxflux(5)          
      integer(i4b)                       :: ibasin, itime, ivar ! index 
      character(len=strLen)              :: filename
      ! initialize error control
      err=0; message='read_lohmann_sim/'
      open (UNIT=53,file=trim(filelist_name),form='formatted',status='old',iostat=err)
      if (err/=0) then; message=trim(message)//"openError['"//trim(filelist_name)//"']";return;endif
      do ibasin = 1,nbasin
        read (UNIT=53,fmt=*) filename
        filename=trim(sim_dir)//trim(filename)
        open (UNIT=50,file=filename,form='formatted',status='old',iostat=err)
        if (err/=0) then; message=trim(message)//"openError['"//trim(filename)//"']";return;endif
        do itime = 1,sim_len
          read (UNIT=50,fmt=*) (auxflux(ivar), ivar=1,4)
          sim(ibasin,itime) = auxflux(4)
        enddo
        close(UNIT=50)
      enddo
      close(UNIT=53)
      return
    end subroutine

  end subroutine

  !private routine
  subroutine generic_hru_id(hruid, err, message)
    use read_ncdata, only: get_vec_ivar
    implicit none
    ! input variables
    ! output variables
    integer(i4b),         intent(out)  :: hruid(:)     !  
    integer(i4b),         intent(out)  :: err          ! error code
    character(len=strLen)              :: message     ! error message from downward subroutine
    ! local variables 
    character(len=strLen)              :: cmessage     ! error message from downward subroutine

    err=0; message="generic_hru_id/"
    call get_vec_ivar(trim(mpr_input_dir)//trim(fname_smapping), trim("hru_id"), 1, nHru, hruid, err, cmessage) 
    if (err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine

end module model_wrapper 
