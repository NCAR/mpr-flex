module write_param_nc 

  ! collection of netCDF writing routines
  
  use netcdf
  use public_var
  use nrtype
  use data_type                                    ! data strucutre definition
  use var_lookup,   only: nPar 
  
  implicit none
  
  private
  public::defSoilNetCDF
  public::write_vec_ivar
  public::write_vec_dvar
  public::write_array2_ivar
  public::write_array2_dvar
  
  ! define dimension names
  character(len=32),parameter :: sHru_DimName='hruid'    ! output dimension name for HRU
  character(len=32),parameter :: sLyr_DimName='lyr'    ! dimension name for layer

contains

  ! *********************************************************************
  ! public subroutine: define output NetCDF file for soil parameters 
  ! *********************************************************************
  subroutine defSoilNetCDF(fname,           &  ! input: output nc filename
                           nHru,            &  ! input: number of model hrus 
                           nLyr,            &  ! input: number of model soil layer 
                           err, message)      ! output: error control
    use globalData,    only: betaInGamma, parMaster 
    implicit none
    ! input variables
    character(*), intent(in)          :: fname        ! filename
    integer(i4b), intent(in)          :: nHru         ! number of hru polygon 
    integer(i4b), intent(in)          :: nLyr         ! number of soil layer 
    ! output variables
    integer(i4b), intent(out)         :: err         ! error code
    character(*), intent(out)         :: message      ! error message
    ! local variables
    type(par_meta),  allocatable      :: betaMeta(:)  ! meta data for beta parameter estimated via MPR 
    integer(i4b)                      :: ncid         ! NetCDF file ID
    integer(i4b)                      :: hrudimID     ! hru dimension ID
    integer(i4b)                      :: lyrdimID     ! layer dimension ID 
    integer(i4b)                      :: iPar         ! variable index
    integer(i4b)                      :: iBeta        ! variable index
    integer(i4b)                      :: nBeta        ! number of variables
    character(len=strLen)             :: cmessage     ! error message of downwind routine
    
    ! initialize error control
    err=0; message='defSoilNetCDF/'
    ! define file
    err = nf90_create(trim(fname),nf90_classic_model,ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! define dimension (unlimited)
    !err = nf90_def_dim(ncid, trim(sHru_DimName), nf90_unlimited, idimId)
    err = nf90_def_dim(ncid, trim(sHru_DimName), nHru, hrudimId)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    err = nf90_def_dim(ncid, trim(sLyr_DimName), nLyr, lyrdimId)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! define coordinate variable - variable name is the same as dimension 
    call defvar(sHru_DimName,'Model HRU ID','-',  (/sHru_DimName/),nf90_int, err,cmessage)
    call defvar(sLyr_DimName,'Soil Layer ID','-', (/sLyr_DimName/),nf90_int, err,cmessage)
    nBeta=size(betaInGamma)
    allocate(betaMeta(nBeta)) 
    do iBeta=1,nBeta
      do iPar=1,nPar
        if ( parMaster(iPar)%pname==betaInGamma(iBeta) )then
          betaMeta(iBeta)=parMaster(iPar)
        endif
      enddo
      ! define parameter values 
      call defvar(trim(betaMeta(iBeta)%pname),trim(betaMeta(iBeta)%pname),"-",(/sLyr_DimName,sHru_DimName/),nf90_double,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    end do
    ! end definitions
    err = nf90_enddef(ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! close NetCDF file
    err = nf90_close(ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    
    contains

    subroutine defvar(vname,vdesc,vunit,dimNames,ivtype,err,message)
      ! input
      character(*), intent(in)  :: vname                      ! variable name
      character(*), intent(in)  :: vdesc                      ! variable description
      character(*), intent(in)  :: vunit                      ! variable units
      character(*), intent(in)  :: dimNames(:)                ! variable dimension names
      integer(i4b), intent(in)  :: ivtype                     ! variable type
      ! output
      integer(i4b), intent(out) :: err                       ! error code
      character(*), intent(out) :: message                    ! error message
      ! local
      integer(i4b)              :: id                         ! index of dimension loop
      integer(i4b)              :: dimIDs(size(dimNames))     ! vector of dimension IDs
      integer(i4b)              :: iVarId                     ! variable ID
      integer(i4b)              :: nf90_fill_int=-999         ! fill value for integer type variable
      real(dp)                  :: nf90_fill_double=-999.0_dp ! fill value for double type variable
     
      ! define dimension IDs
      do id=1,size(dimNames)
       err=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id))
       if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
      end do
      ! define variable
      err = nf90_def_var(ncid,trim(vname),ivtype,dimIds,iVarId)
      if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
      ! add variable description
      err = nf90_put_att(ncid,iVarId,'long_name',trim(vdesc))
      if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
      ! add variable units
      err = nf90_put_att(ncid,iVarId,'units',trim(vunit))
      if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
      ! add variable fillvalue 
      if (ivtype == nf90_double) then
        err = nf90_put_att(ncid, iVarId, '_FillValue', nf90_fill_double)
        if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
      else if (ivtype == nf90_int) then
        err = nf90_put_att(ncid, iVarId, '_FillValue', nf90_fill_int)
        if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
      endif
      return
    end subroutine defvar
  end subroutine

  ! *********************************************************************
  ! public subroutine: write an integer vector
  ! *********************************************************************
  subroutine write_vec_ivar(fname,           &  ! input: filename
                            vname,           &  ! input: variable name
                            iVec,            &  ! input: variable data
                            iStart,          &  ! input: start index
                            err, message)      ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)        :: fname        ! filename
    character(*), intent(in)        :: vname        ! variable name
    integer(i4b), intent(in)        :: iVec(:)      ! variable data
    integer(i4b), intent(in)        :: iStart       ! start index
    ! output variables
    integer(i4b), intent(out)       :: err         ! error code
    character(*), intent(out)       :: message      ! error message
    ! local variables
    integer(i4b)                    :: ncid         ! NetCDF file ID
    integer(i4b)                    :: iVarId       ! NetCDF variable ID

    ! initialize error control
    err=0; message='write_vec_ivar/'
    ! open NetCDF file
    err = nf90_open(trim(fname),nf90_write,ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! assign variable ID
    err = nf90_inq_varid(ncid,trim(vname),iVarId)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! write data
    err = nf90_put_var(ncid,iVarId,iVec,start=(/iStart/),count=(/size(iVec)/))
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! close output file
    err = nf90_close(ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    return 
  end subroutine

  ! *********************************************************************
  ! public subroutine: write a double precision vector
  ! *********************************************************************
  subroutine write_vec_dvar(fname,           &  ! input: filename
                            vname,           &  ! input: variable name
                            dVec,            &  ! input: variable data
                            iStart,          &  ! input: start index
                            iCount,          &  ! input: length of vector
                            err, message)      ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)        :: fname        ! filename
    character(*), intent(in)        :: vname        ! variable name
    real(dp), intent(in)            :: dVec(:)      ! variable data
    integer(i4b), intent(in)        :: iStart       ! start indices
    integer(i4b), intent(in)        :: iCount       ! length of vector
    ! output variables
    integer(i4b), intent(out)       :: err         ! error code
    character(*), intent(out)       :: message      ! error message
    ! local variables
    integer(i4b)                    :: ncid         ! NetCDF file ID
    integer(i4b)                    :: iVarId       ! NetCDF variable ID
    ! initialize error control
    err=0; message='write_vec_dvar/'
    ! open NetCDF file
    err = nf90_open(trim(fname),nf90_write,ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! get variable ID
    err = nf90_inq_varid(ncid,trim(vname),iVarId)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! write data
    err = nf90_put_var(ncid,iVarId,dVec,start=(/iStart/),count=(/iCount/))
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! close output file
    err = nf90_close(ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    return 
  end subroutine

  ! *********************************************************************
  ! public subroutine: write a double precision 2d array 
  ! *********************************************************************
  subroutine write_array2_dvar(fname,        &  ! input: filename
                               vname,        &  ! input: variable name
                               dvar,         &  ! input: variable data
                               iStart,       &  ! input: start index
                               iCount,       &  ! input: length of vector
                               err, message)   ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)       :: fname        ! filename
    character(*), intent(in)       :: vname        ! variable name
    real(dp), intent(in)           :: dvar(:,:)    ! variable data
    integer(i4b), intent(in)       :: iStart(:)    ! start indices
    integer(i4b), intent(in)       :: iCount(:)    ! length of vector
    ! output variables
    integer(i4b), intent(out)      :: err         ! error code
    character(*), intent(out)      :: message      ! error message
    ! local variables
    integer(i4b)                   :: ncid         ! NetCDF file ID
    integer(i4b)                   :: VarId        ! NetCDF variable ID
    ! initialize error control
    err=0; message='write_array2_dvar/'
    ! open NetCDF file
    err = nf90_open(trim(fname),nf90_write,ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! get variable ID
    err = nf90_inq_varid(ncid,trim(vname),VarId)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! write data
    err = nf90_put_var(ncid,VarId,dvar,start=iStart,count=iCount)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! close output file
    err = nf90_close(ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    return 
  end subroutine

  ! *********************************************************************
  ! public subroutine: write an integer 2d array 
  ! *********************************************************************
  subroutine write_array2_ivar(fname,        &  ! input: filename
                               vname,        &  ! input: variable name
                               ivar,         &  ! input: variable data
                               iStart,       &  ! input: start index
                               iCount,       &  ! input: length of vector
                               err, message)   ! output: error control
    implicit none
    ! input variables
    character(*), intent(in)       :: fname        ! filename
    character(*), intent(in)       :: vname        ! variable name
    integer(i4b), intent(in)       :: ivar(:,:)    ! variable data
    integer(i4b), intent(in)       :: iStart(:)    ! start indices
    integer(i4b), intent(in)       :: iCount(:)    ! length of vector
    ! output variables
    integer(i4b), intent(out)      :: err         ! error code
    character(*), intent(out)      :: message      ! error message
    ! local variables
    integer(i4b)                   :: ncid         ! NetCDF file ID
    integer(i4b)                   :: VarId        ! NetCDF variable ID
    ! initialize error control
    err=0; message='write_array2_dvar/'
    ! open NetCDF file
    err = nf90_open(trim(fname),nf90_write,ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! get variable ID
    err = nf90_inq_varid(ncid,trim(vname),VarId)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! write data
    err = nf90_put_var(ncid,VarId,ivar,start=iStart,count=iCount)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! close output file
    err = nf90_close(ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    return 
  end subroutine

end module write_param_nc 
