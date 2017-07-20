module write_param_nc 

  ! netCDF writing routines
  
  use netcdf
  use public_var
  use nrtype
  use data_type                                    ! data strucutre definition
  
  implicit none
  
  private

  public::write_nc_soil
  public::write_nc_veg
  
contains

  ! *********************************************************************
  ! public subroutine: write soil parameter in netCDF 
  ! *********************************************************************
  subroutine write_nc_soil(fname,           & ! input: output nc filename
                           hruID,           & ! input: hruID array 
                           hModel,          & ! input: model layer thickness 
                           parMxyMz,        & ! input: parameter data structure 
                           err, message)      ! output: error control
    use globalData,    only: nSoilBetaModel, soilBetaCalName
    implicit none
    ! input variables
    character(*),   intent(in)            :: fname        ! input: filename
    integer(i4b),   intent(in)            :: hruID(:)     ! Hru ID
    real(dp),       intent(in)            :: hModel(:,:)  ! input: array of model layer thickness at model layer x model hru 
    type(namedvar2),intent(in)            :: parMxyMz(:)  ! input: data structure for model soil parameter at model layer x model hru 
    ! output variables
    integer(i4b), intent(out)             :: err          ! error code
    character(*), intent(out)             :: message      ! error message
    ! local variables
    type(defDim)                          :: hruDim
    type(defDim)                          :: sLyrDim
    real(dp),                 allocatable :: zModel(:,:)  ! array of model layer depth at model layer x model hru 
    character(len=strLen),    allocatable :: betaNames(:) ! parameter name array 
    integer(i4b)                          :: iPar         ! loop index for parameter  
    integer(i4b)                          :: iLyr         ! loop index for model layer 
    character(len=strLen)                 :: cmessage     ! error message of downwind routine

    ! initialize error control
    err=0; message='write_nc_soil/'
    ! Define dimension attributes  
    hruDim%dimName='hruid'
    hruDim%dimdesc='Model HRU ID'
    hruDim%dimunit='-'
    sLyrDim%dimName='lyr'
    sLyrDim%dimDesc='Soil Layer ID'
    sLyrDim%dimUnit='-'
    ! construct z array
    allocate(zModel,source=hModel)
    do iLyr=2,nLyr
      zModel(iLyr,:)=zModel(iLyr-1,:)+zModel(iLyr,:) 
    enddo
    ! construct soil beta parameter name array 
    allocate(betaNames(nSoilBetaModel+1),stat=err);if(err/=0)then;message=trim(message)//'error allocating betaNames';return;endif
    betaNames(1:nSoilBetaModel)=soilBetaCalName
    betaNames(nSoilBetaModel+1)='z'
    ! Define variables
    call defNetCDF(fname, betaNames, nHru, nLyr, hruDim, sLyrDim, err, cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    ! write 1st Dimension
    call write_vec_ivar(fname,trim(hruDim%dimName),hruID,1,err,cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    ! write 2nd Dimension
    call write_vec_ivar(fname,trim(sLyrDim%dimName),(/(iLyr,iLyr=1,nLyr)/),1,err,cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    ! Write soil parameters
    do iPar=1,nSoilBetaModel
      call write_array2_dvar(fname,trim(soilBetaCalName(iPar)),parMxyMz(iPar)%varData,(/1,1/),(/nLyr,nHru/),err,cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    enddo
    ! Write soil layers depth 
    call write_array2_dvar(fname,'z',zModel,(/1,1/),(/nLyr,nHru/),err,cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    return 
  end subroutine
  
  ! *********************************************************************
  ! public subroutine: write veg parameters 
  ! *********************************************************************
  subroutine write_nc_veg(fname,           & ! input: output nc filename
                          hruID,           & ! input: hruID array 
                          vegParMxy,         & ! input: parameter data structure 
                          err, message)      ! output: error control
    use globalData,    only: nVegBetaModel, vegBetaCalName
    implicit none
    ! input variables
    character(*),   intent(in)            :: fname        ! input: filename
    integer(i4b),   intent(in)            :: hruID(:)     ! Hru ID
    type(namedvar2),intent(in)            :: vegParMxy(:) ! input: data structure for model veg parameter at model layer x model hru 
    ! output variables
    integer(i4b), intent(out)             :: err          ! error code
    character(*), intent(out)             :: message      ! error message
    ! local variables
    type(defDim)                          :: hruDim
    type(defDim)                          :: monDim
    integer(i4b)                          :: iPar         ! loop index for parameter  
    integer(i4b)                          :: iMon         ! loop index for month 
    character(len=strLen)                 :: cmessage     ! error message of downwind routine

    ! initialize error control
    err=0; message='write_nc_veg/'
    ! Define dimension attributes  
    hruDim%dimName='hruid'
    hruDim%dimdesc='Model HRU ID'
    hruDim%dimunit='-'
    monDim%dimName='mon'
    monDim%dimDesc='Month'
    monDim%dimUnit='-'
    ! Define variables
    call defNetCDF(fname, vegBetaCalName, nHru, nMonth, HruDim, monDim, err, cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    ! write 1st Dimension
    call write_vec_ivar(fname,trim(hruDim%dimName),hruID,1,err,cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    ! write 2nd Dimension
    call write_vec_ivar(fname,trim(monDim%dimName),(/(iMon,iMon=1,nMonth)/),1,err,cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    ! Write veg parameters
    do iPar=1,nVegBetaModel
      call write_array2_dvar(fname,trim(vegBetaCalName(iPar)),vegParMxy(iPar)%varData,(/1,1/),(/nMonth,nHru/),err,cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    enddo
    return 
  end subroutine
  
  !private subroutine
  subroutine defNetCDF(fname,           &  ! input: output nc filename
                       betaNames,       &  ! input: character array for beta parameter name
                       nDim1,           &  ! input: number of 1st dimension 
                       nDim2,           &  ! input: number of 2nd dimension 
                       defDim1,         &  ! input: 1st dimension name 
                       defDim2,         &  ! input: 2nd dimension name 
                       err, message)       ! output: error control
    use globalData,    only: betaMeta 
    use var_lookup,    only: nBeta 
    implicit none
    ! input variables
    character(*), intent(in)          :: fname           ! filename
    character(*), intent(in)          :: betaNames(:)    ! parameter name array 
    integer(i4b), intent(in)          :: nDim1           ! number of 1st dimension (hru polygon )
    integer(i4b), intent(in)          :: nDim2           ! number of 2nd Dimension (e.g., soil layer, month) 
    type(defDim), intent(in)          :: defDim1         ! 1st dimension name 
    type(defDim), intent(in)          :: defDim2         ! 2nd dimension name 
    ! output variables
    integer(i4b), intent(out)         :: err             ! error code
    character(*), intent(out)         :: message         ! error message
    ! local variables
    type(par_meta),  allocatable      :: betaMetaTemp(:) ! meta data for beta parameter estimated via MPR
    integer(i4b)                      :: ncid            ! NetCDF file ID
    integer(i4b)                      :: dim1ID          ! 1st dimension ID (hru)
    integer(i4b)                      :: dim2ID          ! 2nd dimension ID (soil layer, month) 
    integer(i4b)                      :: nBetaOut        ! number of beta parameter to be output 
    integer(i4b)                      :: iPar            ! variable index
    integer(i4b)                      :: iBeta           ! variable index
    character(len=strLen)             :: cmessage        ! error message of downwind routine
    
    ! initialize error control
    err=0; message='defNetCDF/'
    ! define file
    err = nf90_create(trim(fname),nf90_classic_model,ncid)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! define dimension (unlimited)
    !err = nf90_def_dim(ncid, trim(sHru_DimName), nf90_unlimited, idimId)
    err = nf90_def_dim(ncid, trim(defDim1%dimName), nDim1, dim1id)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    err = nf90_def_dim(ncid, trim(defDim2%dimName), nDim2, dim2id)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
    ! define coordinate variable - variable name is the same as dimension 
    call defvar(defDim1%dimName,trim(defDim1%dimDesc),'-', (/defDim1%dimName/), nf90_int, err,cmessage)
    call defvar(defDim2%dimName,trim(defDim2%dimDesc),'-', (/defDim2%dimName/), nf90_int, err,cmessage)
    nBetaOut=size(betaNames)
    allocate(betaMetaTemp(nBetaOut)) 
    do iBeta=1,nBetaOut
      do iPar=1,nBeta
        if ( betaMeta(iPar)%pname==betaNames(iBeta) )then;betaMetaTemp(iBeta)=betaMeta(iPar); exit; endif
      enddo
      ! define parameter values 
      call defvar(trim(betaMetaTemp(iBeta)%pname),trim(betaMetaTemp(iBeta)%pname),'-',(/defDim2%dimName,defDim1%dimName/),nf90_double,err,cmessage)
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

  !private subroutine
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

  ! private subroutine: write a double precision vector
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

  ! private subroutine: write a double precision 2d array 
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

  ! private subroutine: write an integer 2d array 
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
