module process_meta 

  use nrtype
  use public_var

  implicit none

  private

  public::read_inParList
  public::get_parm_meta
  public::param_setup 
  public::print_config

contains

! ************************************************************************************************
! Public subroutine: Read 'inParList' input file 
! ************************************************************************************************
subroutine read_inParList(infile, err, message)
  use globalData, only:inParMeta      ! meta for inParList nml input  
  use data_type,  only:input_meta
  use ascii_util, only:file_open
  use var_lookup, only:nBeta
  implicit none
  ! input
  character(*),intent(in)              :: infile            ! input filename
  ! output
  integer(i4b),intent(out)             :: err               ! error code
  character(*),intent(out)             :: message           ! error message
  ! local variables
  type(input_meta),allocatable         :: tempCalParMeta(:) ! temp data structure for inParList input meta
  character(len=strLen)                :: cmessage          ! error message subroutine
  integer(i4b)                         :: ixLocal           ! index for calibrationg parameter list 
  integer(i4b),parameter               :: maxLines=1000     ! maximum lines in the file 
  integer(i4b)                         :: iend              ! check for the end of the file
  integer(i4b)                         :: unt               ! DK: need to either define units globally, or use getSpareUnit
  character(LEN=strLen)                :: temp              ! single lime of information
  character(LEN=strLen)                :: ffmt              ! file format
  character(len=1)                     :: dLim(4)           ! column delimiter
  integer(i4b)                         :: iline             ! loop through lines in the file 

  ! initialize error handling 
  err=0; message="read_inParList/"
  allocate(tempCalParMeta(nBeta))
  call file_open(trim(infile), unt, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! get to the start of the variable descriptions 
  do iline=1,maxLines
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
    if (temp(1:1)/='!') exit                             ! assume first line not comment is format code
  end do ! looping through file to find the format code
  read(temp,*)ffmt ! to get format 
  ixLocal=0_i4b
  line:do iline=1,maxLines
    ! read a line of data and exit iif an error code (character read, so only possible error is end of file)
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
    ! check that the line is not a comment
    if (temp(1:1)=='!')cycle
    ! save data into a temporary structure
    ixLocal = ixLocal+1_i4b
    read(temp,trim(ffmt),iostat=err) tempCalParMeta(ixLocal)%betaname,    dLim(1),&    ! beta parameter name
                                     tempCalParMeta(ixLocal)%calMethod,   dLim(2),&    ! Parameter estimation mehtod - 0, skip calibration, 1. MPR or 2. mulitplier,  
                                     tempCalParMeta(ixLocal)%TF,          dLim(3),&    ! Transfer function type
                                     tempCalParMeta(ixLocal)%isScaleCalH, dLim(4),&    ! calibrating scaling operator for horizontal direction?
                                     tempCalParMeta(ixLocal)%isScaleCalV               ! calibrating scaling operator for vertical direction?
    if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
    ! check that the delimiters are in the correct place
    if(any(dLim /= '|'))then
      message=trim(message)//'delimiter is not in the correct place; line = ['//trim(temp)//']; filename = '//trim(infile)
      err=32; return
    endif
  enddo line
  ! close file unit
  close(unt)
  ! save 'inParMeta'
  allocate(inParMeta(ixLocal))
  inParMeta=tempCalParMeta(1:ixLocal) 
  return
end subroutine

! ************************************************************************************************
! Public subroutine: Prepare calibrating parameter metadata from a meta file 
! ************************************************************************************************
subroutine get_parm_meta( err, message)
  ! Process inParMeta along with betaMeta and gammaMeta (popMeta.f90) 
  ! Saved data:  calParMeta, calGammaMeta, calBetaName 
  use data_type,  only:par_meta, cpar_meta, scale_meta
  use globalData, only:inParMeta,      & ! meta for beta parameter listed in 'inParList' nml input
                       gammaMeta,      & ! meta for all gamma parameters
                       betaMeta,       & ! meta for all beta parameters
                       calParMeta,     & ! meta for only parameter listed in input
                       calGammaMeta,   & ! meta for only gamma parameters listed in input
                       calBetaName,    & ! list of beta parameters calibrated with MPR 
                       calScaleMeta,   & ! meta for beta parameter whose scaling operator(s) is calibrated
                       nCalPar,        & ! sum of beta and gamma parameters to be calibrated 
                       nCalBetaDir,    & ! number of beta parameters to be calibrated "directly"
                       nCalGamma,      & ! number of gamma parameters to be calibrated 
                       nSoilBetaModel, & ! number of soil beta parameters to be calibrated via MPR  
                       nVegBetaModel,  & ! number of veg beta parameters to be calibrated via MPR
                       nCalParSum,     & ! number of total calibrating parameters going to optimization routine
                       soilBetaCalName,& ! list of Soil parameters to be estimated via MPR
                       vegBetaCalName    ! list of vege parameters to be estimated via MPR
  use get_ixname, only:get_ixBeta, get_ixGamma
  use var_lookup, only:nBeta,nGamma
  implicit none
  ! input
  ! output
  integer(i4b),intent(out)             :: err              ! error code
  character(*),intent(out)             :: message          ! error message
  ! local variables
  character(len=strLen)                :: cmessage         ! error message from subroutine
  integer(i4b)                         :: iBeta            ! loop index of lines in inParList nml input
  integer(i4b)                         :: ivar             ! loop index of master parameter  
  integer(i4b)                         :: iPar             ! loop index of master parameter  
 
  err=0; message="get_param_meta/"
  ! update betaMeta(:)%tftype (transfer function type to be used)
  do iBeta=1,size(inParMeta)
    ivar=get_ixBeta(inParMeta(iBeta)%betaname)
    if(ivar<=0)then; err=40; message=trim(message)//"1.variableNotFound[var="//trim(inParMeta(iBeta)%betaname)//"]"; return; endif
    betaMeta(ivar)%tftype=inParMeta(iBeta)%TF
  enddo
  call get_parMetaSubset( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_calBetaName  ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_calScaleMeta ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call total_calParam   ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif

  return

  contains

  ! Private subroutine:   
  subroutine get_parMetaSubset( err, message )
   ! Saved data:  calParMeta, calGammaMeta, nCalPar, nCalBetaDir, nCalGamma
    !input
    !output
    character(len=strLen),intent(out)   :: message           ! error message for current routine
    integer(i4b),         intent(out)   :: err               ! error code
    !local
    integer(i4b)                         :: ivar             ! loop index of master parameter  
    integer(i4b)                         :: iPar            ! loop index of lines in "inParList" nml input
    integer(i4b)                         :: iGamma           ! loop index of all the gamma parameters in master 
    type(cpar_meta),allocatable          :: tempParSubset(:)
    type(cpar_meta),allocatable          :: tempGammaMeta(:)

    err=0; message="get_parMetaSubset/"
    allocate(tempParSubset(nGamma+nBeta))
    allocate(tempGammaMeta(nGamma+nBeta))
    nCalPar=0 ! Count and Save number of sum of beta (directly) and gamma parameters to be calibrated
    nCalGamma=0     ! Count and Save number of gamma parameter to be calibrated
    do iPar=1,size(inParMeta)
      if(inParMeta(iPar)%calMethod==1 .or. inParMeta(iPar)%calMethod==2) then 
        do iGamma=1,nGamma  ! look for gammma parameters associated with this beta parameter
          if (gammaMeta(iGamma)%beta==inParMeta(iPar)%betaname .and. gammaMeta(iGamma)%tftype==inParMeta(iPar)%TF) then
            nCalPar = nCalPar+1
            ivar=get_ixGamma(gammaMeta(iGamma)%pname)
            tempParSubset(nCalPar)%ixMaster = ivar
            tempParSubset(nCalPar)%pname    = gammaMeta(iGamma)%pname
            tempParSubset(nCalPar)%val      = gammaMeta(iGamma)%val
            tempParSubset(nCalPar)%lwr      = gammaMeta(iGamma)%lwr
            tempParSubset(nCalPar)%upr      = gammaMeta(iGamma)%upr
            tempParSubset(nCalPar)%beta     = gammaMeta(iGamma)%beta
            tempParSubset(nCalPar)%tftype   = gammaMeta(iGamma)%tftype
            tempParSubset(nCalPar)%ptype    = gammaMeta(iGamma)%ptype
            if(inParMeta(iPar)%calMethod==1) then !if calMethod is MPR
              tempParSubset(nCalPar)%flag     = .True.
            elseif(inParMeta(iPar)%calMethod==2) then !if calibrating only scale parameters 
              tempParSubset(nCalPar)%flag     = .False.
            endif
            tempParSubset(nCalPar)%hups     = gammaMeta(iGamma)%hups
            tempParSubset(nCalPar)%hpnorm   = gammaMeta(iGamma)%hpnorm
            tempParSubset(nCalPar)%vups     = gammaMeta(iGamma)%vups
            tempParSubset(nCalPar)%vpnorm   = gammaMeta(iGamma)%vpnorm
            tempParSubset(nCalPar)%perLyr   = gammaMeta(iGamma)%perLyr
            nCalGamma = nCalGamma+1
            tempGammaMeta(nCalGamma) = tempParSubset(nCalPar)
          endif
        end do
      elseif (inParMeta(iPar)%calMethod==3) then ! if calMethod is direct calibration
        ivar = get_ixBeta(inParMeta(iPar)%betaname)
        if(ivar<=0)then; err=40; message=trim(message)//"2.variableNotFound[var="//trim(inParMeta(iPar)%betaname)//"]"; return; endif
        nCalPar = nCalPar+1
        tempParSubset(nCalPar)%ixMaster = ivar
        tempParSubset(nCalPar)%pname    = betaMeta(ivar)%pname
        tempParSubset(nCalPar)%val      = betaMeta(ivar)%val
        tempParSubset(nCalPar)%lwr      = betaMeta(ivar)%lwr
        tempParSubset(nCalPar)%upr      = betaMeta(ivar)%upr
        tempParSubset(nCalPar)%beta     = betaMeta(ivar)%beta
        tempParSubset(nCalPar)%tftype   = betaMeta(ivar)%tftype
        tempParSubset(nCalPar)%ptype    = betaMeta(ivar)%ptype
        tempParSubset(nCalPar)%flag     = .True. 
        tempParSubset(nCalPar)%hups     = betaMeta(ivar)%hups
        tempParSubset(nCalPar)%hpnorm   = betaMeta(ivar)%hpnorm
        tempParSubset(nCalPar)%vups     = betaMeta(ivar)%vups
        tempParSubset(nCalPar)%vpnorm   = betaMeta(ivar)%vpnorm
        tempParSubset(nCalPar)%perLyr   = betaMeta(ivar)%perLyr
      endif
    enddo
    nCalBetaDir=nCalPar-nCalGamma      ! Save number of beta parameter to be directly calibrated
    if (nCalPar > 0_i4b) then          ! Save 'calParMeta'
      allocate(calParMeta(nCalPar))
      calParMeta=tempParSubset(1:nCalPar) 
    endif
    if (nCalGamma > 0_i4b) then              ! Save'calGammaMeta'
      allocate(calGammaMeta(nCalGamma))
      calGammaMeta=tempGammaMeta(1:nCalGamma) 
    endif
    return
  end subroutine

  ! Private subroutine: Get list of beta parameters (excluding z and h) computed with MPR.  
  ! populated global data: calBetaName, soilBetaCalName, vegBetaCalName
  subroutine get_calBetaName( err, message )
    !input
    !output
    character(len=strLen),intent(out)    :: message                ! error message for current routine
    integer(i4b),         intent(out)    :: err                    ! error code
    !local
    integer(i4b)                         :: iPar                   ! loop index of parameter 
    integer(i4b)                         :: ivar                   ! index in master parameter
    integer(i4b)                         :: iSoil                  ! counter for soil beta parameters to be computed with MPR 
    integer(i4b)                         :: iVeg                   ! counter for vege beta parameters to be computed with MPR
    character(len=strLen),allocatable    :: tempSoilBetaInGamma(:) ! 
    character(len=strLen),allocatable    :: tempVegBetaInGamma(:)  ! 
    character(len=strLen),allocatable    :: res(:)                 ! 
    integer(i4b)                         :: counter                 ! counter 

    err=0; message="get_calBetaName/"
    allocate(res(size(betaMeta)))
    counter = 0
    do iPar=1,size(inParMeta)
      if(inParMeta(iPar)%calMethod==1 .or. inParMeta(iPar)%calMethod==2) then 
        ivar = get_ixBeta(inParMeta(iPar)%betaname)
        if(ivar<=0)then; err=40; message=trim(message)//"BetaNotFoundInMasterMeta[var="//trim(inParMeta(iPar)%betaname)//"]"; return; endif
        if ( trim(inParMeta(iPar)%betaname) == 'h1' ) cycle 
        if ( trim(inParMeta(iPar)%betaname) == 'h2' ) cycle
        if ( trim(inParMeta(iPar)%betaname) == 'h3' ) cycle
        if ( trim(inParMeta(iPar)%betaname) == 'h4' ) cycle
        if ( trim(inParMeta(iPar)%betaname) == 'h5' ) cycle
        if ( trim(inParMeta(iPar)%betaname) == 'z' )  cycle
        counter=counter+1
        res(counter)=betaMeta(ivar)%pname
      end if
    end do
    if ( counter > 0 ) then
      allocate(calBetaName(counter))
      calBetaName=res(1:counter)
      ! Count number of soil and Vege parameters to be computed with MPR 
      nSoilBetaModel=0
      nVegBetaModel=0
      do iPar=1,size(calBetaName)
        ivar=get_ixBeta(calBetaName(iPar))
        if (betaMeta(ivar)%ptype=='soil') nSoilBetaModel = nSoilBetaModel+1
        if (betaMeta(ivar)%ptype=='veg')  nVegBetaModel  = nVegBetaModel +1
      enddo
      iSoil=0;iVeg=0
      allocate(tempSoilBetaInGamma(nBeta),stat=err)
      allocate(tempVegBetaInGamma(nBeta),stat=err)
      do iPar=1,size(calBetaName)
        ivar = get_ixBeta(calBetaName(iPar))
        if(ivar<=0)then; err=40; message=trim(message)//"2.variableNotFound[var="//trim(calBetaName(iPar))//"]"; return; endif
        if (betaMeta(ivar)%ptype=='soil')then; iSoil=iSoil+1; tempSoilBetaInGamma(iSoil) = betaMeta(ivar)%pname; endif
        if (betaMeta(ivar)%ptype=='veg') then; iVeg=iVeg+1;   tempVegBetaInGamma(iVeg)   = betaMeta(ivar)%pname; endif
      enddo 
      allocate(soilBetaCalName(nSoilBetaModel))
      soilBetaCalName=tempSoilBetaInGamma(1:nSoilBetaModel) 
      allocate(vegBetaCalName(nVegBetaModel))
      vegBetaCalName=tempVegBetaInGamma(1:nVegBetaModel) 
    else
      print*, 'NO gamma parameters included in the list'
    endif
    return
  end subroutine

  ! Private subroutine: to populate scaling parameter meta for beta parameter calibrated with MPR
  ! populated global data: calScaleMeta
  subroutine get_calScaleMeta( err, message ) 
    !input
    !output
    character(len=strLen),intent(out)   :: message               ! error message for current routine
    integer(i4b),         intent(out)   :: err                   ! error code
    !local
    type(scale_meta),     allocatable   :: tempBetaCalScale(:)   ! 
    integer(i4b)                        :: nScaleBeta            ! counter 

    err=0; message="get_scaleInBeta/"
    allocate(tempBetaCalScale(size(inParMeta)),stat=err);if(err/=0)then;message=trim(message)//'error allocating tempBetaCalScale';return;endif
    nScaleBeta=0 
    do iPar=1,size(inParMeta)!if beta parameter is estimated with MPR or even not calibrated, calibrate scaling parameter
      if ( inParMeta(iPar)%calMethod==1 .or. inParMeta(iPar)%calMethod==2 ) then  
        associate( hups      => betaMeta(get_ixBeta(inParMeta(iPar)%betaname))%hups,   &
                   vups      => betaMeta(get_ixBeta(inParMeta(iPar)%betaname))%vups,   &
                   hpower    => betaMeta(get_ixBeta(inParMeta(iPar)%betaname))%hpnorm, &
                   vpower    => betaMeta(get_ixBeta(inParMeta(iPar)%betaname))%vpnorm, &
                   hscaleMask=> inParMeta(iPar)%isScaleCalH, &
                   vscaleMask=> inParMeta(iPar)%isScaleCalV )
        nScaleBeta=nScaleBeta+1
        tempBetaCalScale(nScaleBeta)%betaname    = inParMeta(iPar)%betaname  
        tempBetaCalScale(nScaleBeta)%pdefault(1) = hpower 
        tempBetaCalScale(nScaleBeta)%pdefault(2) = vpower 
        tempBetaCalScale(nScaleBeta)%mask(1)     = hscaleMask 
        tempBetaCalScale(nScaleBeta)%mask(2)     = vscaleMask 
        ! Check beta meta data for scaling and if value is -999 or na, turns off calibration 
        if (hscaleMask .and. (hpower==-999.0_dp .or. hups=='na'))then 
          print*,'Switch horizontal upscale calibration to False for',inParMeta(iPar)%betaname
          tempBetaCalScale(nScaleBeta)%mask(1)=.False.
        endif
        if (vscaleMask .and. (vpower==-999.0_dp .or. vups=='na'))then
          print*,'Switch vertical upscale calibration to False for',inParMeta(iPar)%betaname
          tempBetaCalScale(nScaleBeta)%mask(2)=.False.
        endif
        end associate
      endif
    enddo
    allocate(calScaleMeta(nScaleBeta),stat=err);if(err/=0)then;message=trim(message)//'error allocating calScaleMeta';return;endif
    calScaleMeta=tempBetaCalScale(1:nScaleBeta)
    return
  end subroutine
  
  ! Private subroutine: to count calibrating parameters including:
  ! beta parameters (direct calibration) per layer, gamma parameters, scaling parameters (V and H) 
  ! save nCalParSum in public_var
  subroutine total_calParam( err, message )
    implicit none
    !output variables
    integer(i4b),         intent(out) :: err         ! error code
    character(*),         intent(out) :: message     ! error message
    ! local variables
    integer(i4b)                 :: iPar          ! loop indices
    
    ! initialize error control
    err=0; message='total_calParam/'
    nCalParSum=0_i4b
    do iPar=1,nCalPar        ! cound number of calibration beta (per layer) and gamma
      if (calParMeta(iPar)%perLyr)then
        nCalParSum=nCalParSum+nLyr
      else
        nCalParSum=nCalParSum+1
      endif
    enddo  
    do iPar=1,size(inParMeta) ! Add number of calibration pnorm value - make False if value is -999 or nan
      if (inParMeta(iPar)%calMethod==1 .or. inParMeta(iPar)%calMethod==2 ) nCalParSum=nCalParSum+2
    enddo
    return
  end subroutine

end subroutine

! ************************************************************************************************
! Public subroutine: convert parameter data structure to simple arrays 
! ************************************************************************************************
subroutine param_setup( err, message )
  use globalData,  only:parArray, parMask, calParMeta, nCalPar, calScaleMeta, nCalParSum 
  implicit none
  !output variables
  integer(i4b),     intent(out) :: err                    ! error code
  character(*),     intent(out) :: message                ! error message
  ! local variables
  integer(i4b)                  :: iPar                   ! loop indices
  integer(i4b)                  :: idx                    ! count of calibrating parameter including per layer parameter 
  integer(i4b)                  :: ixHV                   ! count of calibrating parameter including per layer parameter 
  
  ! initialize error control
  err=0; message='param_setput/'
  allocate(parArray(nCalParSum,3),stat=err);if(err/=0)then;message=trim(message)//'error allocating parArray';return;endif
  allocate(parMask(nCalParSum),stat=err);if(err/=0)then;message=trim(message)//'error allocating parMask';return;endif
  idx=0_i4b
  do iPar=1,nCalPar
    if (calParMeta(iPar)%perLyr)then
      idx=idx+nLyr
      parArray(idx-nLyr+1:idx,1) = calParMeta(iPar)%val
      parArray(idx-nLyr+1:idx,2) = calParMeta(iPar)%lwr
      parArray(idx-nLyr+1:idx,3) = calParMeta(iPar)%upr
      parMask (idx-nLyr+1:idx)   = calParMeta(iPar)%flag
    else
      idx=idx+1
      parArray(idx,1) = calParMeta(iPar)%val
      parArray(idx,2) = calParMeta(iPar)%lwr
      parArray(idx,3) = calParMeta(iPar)%upr
      parMask (idx)   = calParMeta(iPar)%flag
    endif
  enddo  
  do iPar=1,size(calScaleMeta) 
    do ixHV=1,2
      idx=idx+1
      parArray(idx,1) = calScaleMeta(iPar)%pdefault(ixHV)        ! default value of pnorm value
      parArray(idx,2) = -100.0_dp                                ! lower bound or pnorm value
      parArray(idx,3) =  100.0_dp                                ! upper bound of pnorm value
      parMask (idx)   = calScaleMeta(iPar)%mask(ixHV)
    enddo
  enddo
  return
end subroutine

!*********************************************************
! Public subroutine: print out calibrating parameter data  
!*********************************************************
subroutine print_config()
  use globaldata,  only: inParMeta,      &
                         betaMeta,       &
                         calParMeta,      &
                         calGammaMeta,    &
                         calScaleMeta,   & 
                         calBetaName,    &
                         calBetaOrderIdx,&
                         parMask,        &
                         parArray,       &
                         nCalPar,  &
                         nBetaNeed
  implicit none

  integer(i4b) :: i,j   ! loop index for writing

  write(*,*) '!-----------------------------------------------------------'
  write(*,*) '!    MPR-flex - configurations of parameter estimations     '
  write(*,*) '!-----------------------------------------------------------'
  write(*,'(A,1X,A)') new_line(' '),'! Beta parameters listed in input'
  write(*,*) '!-----------------------------------------------------------'
  do i=1,size(inParMeta)
    write(*,*) trim(adjustl(inParMeta(i)%betaname))
  end do
  write(*,'(A,1X,A)') new_line(' '),'! Calibrating Beta and Gamma parameters'
  write(*,*) '!-----------------------------------------------------------'
  do i=1,size(calParMeta)
    write(*,*) ( trim(adjustl(calParMeta(i)%pname)) )
  end do
  if (size(calBetaName)/=0)then
    write(*,'(A,1X,A)') new_line(' '),'! Beta parameters to be estimated with MPR excluding z and h'
    write(*,*) '!-----------------------------------------------------------'
    do i=1,size(calBetaName)
      write(*,*) ( trim(adjustl(calBetaName(i))) )
    end do
    write(*,'(A,1X,A)') new_line(' '),'! List of gamma parameters calibrated'
    write(*,*) '!-----------------------------------------------------------'
    do i=1,size(calGammaMeta)
      write(*,*) ( trim(adjustl(calGammaMeta(i)%pname)) )
    end do
    write(*,'(A,1X,A)') new_line(' '),'! All beta parameters computed with MPR including dependent beta parameters'
    write(*,*) '!-----------------------------------------------------------'
    do i=1,nBetaNeed
      write(*,*) ( trim(adjustl(betaMeta(calBetaOrderIdx(i))%pname)) )
    end do
  else
    write(*,'(A,1X,A)') new_line(' '), '! No beta parameters estimated with MPR' 
  endif
  write(*,'(A,1X,A)') new_line(' '),'! Parameter Array input to optimization routine' 
  write(*,*) '!-----------------------------------------------------------'
  write(*,*) 'Parameter Name        (initial)value    cal.flag   Note'
  do i=1,nCalPar
    if (calParMeta(i)%perLyr)then
      do j=1,nLyr
        write(*,100) calParMeta(i)%pname(1:20), parArray(i,1), parMask(i), 'layer=', j 
        100 format(1X,A,1X,ES17.10,1X,L9,1X,'layer=',I2)
      end do
    else
      write(*,200) calParMeta(i)%pname(1:20), parArray(i,1), parMask(i)
      200 format(1X,A,1X,ES17.10,1X,L9)
    endif
  enddo  
  do i=1,size(calScaleMeta) 
     write(*,300) calScaleMeta(i)%betaname(1:20), parArray(nCalPar+2*i-1,1), parMask(nCalPar+2*i-1), 'Horizontal scaling parameter'
     write(*,300) calScaleMeta(i)%betaname(1:20), parArray(nCalPar+2*i  ,1), parMask(nCalPar+2*i),   'Vertical scaling parameter'
     300 format(1X,A,1X,ES17.10,1X,L9,1X,A30)
  end do
  print*,"!-----------------------------------------------------------"
  print*,"!-----------------------------------------------------------"
  return
end subroutine

!**********************************
! Not Used---  Public subroutine: check if h parameters exist in gamma parameter
!**********************************
subroutine check_gammaH( err, message)
  use globalData,   only: calGammaMeta
  implicit none
  !output variables
  integer(i4b),         intent(out) :: err         ! error code
  character(*),         intent(out) :: message     ! error message
  !local variables
  integer(i4b)                      :: id(20) 
  logical(lgc)                      :: mask(20) 
  logical(lgc),allocatable          :: checkH(:) 
  integer(i4b)                      :: i 

  ! initialize error control
  err=0; message='check_gammaH/'
  allocate(checkH(nLyr-1))
  id=-999
  if ( allocated(calGammaMeta) )then
    !check h parameters - now can chcek up to 5 layers
    do i=1,size(calGammaMeta)
      if (calGammaMeta(i)%pname=="h1gamma1")then;id(1)=1;cycle;endif 
      if (calGammaMeta(i)%pname=="h2gamma1")then;id(2)=1;cycle;endif
      if (calGammaMeta(i)%pname=="h3gamma1")then;id(3)=1;cycle;endif
      if (calGammaMeta(i)%pname=="h4gamma1")then;id(4)=1;cycle;endif
    enddo
    mask=(id>0)
    checkH=mask(1:nLyr-1) 
    if ( any(.not. checkH) )then;err=40;message=trim(message)//"Calibrating gamma parameter require (nLyr-1) hgamma parameters"; return;endif
  else
    print*, 'No gamma parameters listed in CalPar->No MPR computation'
  endif
  return
end subroutine

end module process_meta 
