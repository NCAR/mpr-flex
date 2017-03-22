module process_meta 

  use nrtype
  use public_var

  implicit none

  private

  public::read_calPar
  public::get_parm_meta
  public::param_setup 

contains

! ************************************************************************************************
! Public subroutine: Read 'calPar' input file 
! ************************************************************************************************
subroutine read_calPar(infile, err, message)
  use globalData, only:calParMeta      ! meta for calPar input  
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
  type(input_meta),allocatable         :: tempCalParMeta(:) ! temp data structure for calPar input meta
  character(len=strLen)                :: cmessage          ! error message subroutine
  integer(i2b)                         :: ixLocal           ! index for calibrationg parameter list 
  integer(i4b),parameter               :: maxLines=1000     ! maximum lines in the file 
  integer(i4b)                         :: iend              ! check for the end of the file
  integer(i4b)                         :: unt               ! DK: need to either define units globally, or use getSpareUnit
  character(LEN=strLen)                :: temp              ! single lime of information
  character(LEN=strLen)                :: ffmt              ! file format
  character(len=1)                     :: dLim(3)           ! column delimiter
  integer(i4b)                         :: iline             ! loop through lines in the file 

  ! initialize error handling 
  err=0; message="read_calPar/"
  allocate(tempCalParMeta(nBeta))
  call file_open(trim(infile), unt, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! get to the start of the variable descriptions 
  do iline=1,maxLines
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
    if (temp(1:1)/='!') exit                             ! assume first line not comment is format code
  end do ! looping through file to find the format code
  read(temp,*)ffmt ! to get format 
  ixLocal=0_i2b
  line:do iline=1,maxLines
    ! read a line of data and exit iif an error code (character read, so only possible error is end of file)
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
    ! check that the line is not a comment
    if (temp(1:1)=='!')cycle
    ! save data into a temporary structure
    ixLocal = ixLocal+1_i2b
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
  ! save 'calParMeta'
  allocate(calParMeta(ixLocal))
  calParMeta=tempCalParMeta(1:ixLocal) 
  return
end subroutine

! ************************************************************************************************
! Public subroutine: Prepare calibrating parameter metadata from a meta file 
! ************************************************************************************************
subroutine get_parm_meta( err, message)
  ! Process calParMeta along with betaMaster and gammaMaster (popMeta.f90) 
  ! Saved data:  parSubset, gammaSubset, betaInGamma data structure 
  use data_type,  only:par_meta, cpar_meta, scale_meta
  use globalData, only:calParMeta,     & ! meta for beta parameter listed in 'calPar' input
                       gammaMaster,    & ! meta for all gamma parameters
                       betaMaster,     & ! meta for all beta parameters
                       parSubset,      & ! meta for only parameter listed in input
                       gammaSubset,    & ! meta for only gamma parameters listed in input
                       betaInGamma,    & ! list of beta parameters calibrated with MPR 
                       betaCalScale,   & ! meta for beta parameter whose scaling operator(s) is calibrated
                       nBetaGammaCal,  & ! sum of beta and gamma parameters to be calibrated 
                       nBetaDirCal,    & ! number of beta parameters to be calibrated "directly"
                       nGammaCal,      & ! number of gamma parameters to be calibrated 
                       nSoilParModel,  & ! number of soil beta parameters to be calibrated via MPR  
                       nVegParModel,   & ! number of veg beta parameters to be calibrated via MPR
                       nParCalSum,     & ! number of total calibrating parameters going to optimization routine
                       soilBetaInGamma,& ! list of Soil parameters to be estimated via MPR
                       vegBetaInGamma    ! list of vege parameters to be estimated via MPR
  use get_ixname, only:get_ixBeta, get_ixGamma
  use var_lookup, only:nBeta,nGamma
  implicit none
  ! input
  ! output
  integer(i4b),intent(out)             :: err              ! error code
  character(*),intent(out)             :: message          ! error message
  ! local variables
  character(len=strLen)                :: cmessage         ! error message from subroutine
  integer(i4b)                         :: iBeta            ! loop index of lines in calPar input file 
  integer(i4b)                         :: ivar             ! loop index of master parameter  
  integer(i4b)                         :: iPar             ! loop index of master parameter  
 
  err=0; message="get_param_meta/"
  ! update betaMaster(:)%tftype (transfer function type to be used)
  do iBeta=1,size(calParMeta)
    ivar=get_ixBeta(calParMeta(iBeta)%betaname)
    if(ivar<=0)then; err=40; message=trim(message)//"1.variableNotFound[var="//trim(calParMeta(iBeta)%betaname)//"]"; return; endif
    betaMaster(ivar)%tftype=calParMeta(iBeta)%TF
  enddo
  call get_parMetaSubset( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_betaInGamma  ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_betaCalScale ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call total_calParam   ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  ! check that all elements are populated
  if(any(parSubset(:)%pname==''))then
    do iPar=1,size(parSubset)
      print*,iPar,' -> ',trim(parSubset(iPar)%pname)
    end do
    err=40; message=trim(message)//"parSubset is not filled out completely"; return
  endif
  return

  contains

  ! Private subroutine:   
  subroutine get_parMetaSubset( err, message )
   ! Saved data:  parSubset, gammaSubset, nBetaGammaCal, nBetaDirCal, nGammaCal
    !input
    !output
    character(len=strLen),intent(out)   :: message           ! error message for current routine
    integer(i4b),         intent(out)   :: err               ! error code
    !local
    integer(i4b)                         :: ivar             ! loop index of master parameter  
    integer(i4b)                         :: iBeta            ! loop index of lines in calPar input file 
    integer(i4b)                         :: iGamma           ! loop index of all the gamma parameters in master 
    logical(lgc),allocatable             :: mask(:)
    type(cpar_meta),allocatable          :: tempParSubset(:)
    type(cpar_meta),allocatable          :: tempGammaMeta(:)

    err=0; message="get_parMetaSubset/"
    allocate(tempParSubset(nGamma+nBeta))
    allocate(tempGammaMeta(nGamma+nBeta))
    nBetaGammaCal=0 ! Count and Save number of sum of beta (directly) and gamma parameters to be calibrated
    nGammaCal=0     ! Count and Save number of gamma parameter to be calibrated
    do iBeta=1,size(calParMeta)
      if(calParMeta(iBeta)%calMethod==1 .or. calParMeta(iBeta)%calMethod==2) then 
        do iGamma=1,nGamma  ! look for gammma parameters associated with this beta parameter
          if (gammaMaster(iGamma)%beta==calParMeta(iBeta)%betaname .and. gammaMaster(iGamma)%tftype==calParMeta(iBeta)%TF) then
            nBetaGammaCal = nBetaGammaCal+1
            ivar=get_ixGamma(gammaMaster(iGamma)%pname)
            tempParSubset(nBetaGammaCal)%ixMaster = ivar
            tempParSubset(nBetaGammaCal)%pname    = gammaMaster(iGamma)%pname
            tempParSubset(nBetaGammaCal)%val      = gammaMaster(iGamma)%val
            tempParSubset(nBetaGammaCal)%lwr      = gammaMaster(iGamma)%lwr
            tempParSubset(nBetaGammaCal)%upr      = gammaMaster(iGamma)%upr
            tempParSubset(nBetaGammaCal)%beta     = gammaMaster(iGamma)%beta
            tempParSubset(nBetaGammaCal)%tftype   = gammaMaster(iGamma)%tftype
            tempParSubset(nBetaGammaCal)%ptype    = gammaMaster(iGamma)%ptype
            if(calParMeta(iBeta)%calMethod==1) then !if calMethod is MPR
              tempParSubset(nBetaGammaCal)%flag     = .True.
            elseif(calParMeta(iBeta)%calMethod==2) then !if calibrating only scale parameters 
              tempParSubset(nBetaGammaCal)%flag     = .False.
            endif
            tempParSubset(nBetaGammaCal)%hups     = gammaMaster(iGamma)%hups
            tempParSubset(nBetaGammaCal)%hpnorm   = gammaMaster(iGamma)%hpnorm
            tempParSubset(nBetaGammaCal)%vups     = gammaMaster(iGamma)%vups
            tempParSubset(nBetaGammaCal)%vpnorm   = gammaMaster(iGamma)%vpnorm
            tempParSubset(nBetaGammaCal)%perLyr   = gammaMaster(iGamma)%perLyr
            nGammaCal = nGammaCal+1
            tempGammaMeta(nGammaCal) = tempParSubset(nBetaGammaCal)
          endif
        end do
      elseif (calParMeta(iBeta)%calMethod==3) then ! if calMethod is direct calibration
        ivar = get_ixBeta(calParMeta(iBeta)%betaname)
        if(ivar<=0)then; err=40; message=trim(message)//"2.variableNotFound[var="//trim(calParMeta(iBeta)%betaname)//"]"; return; endif
        nBetaGammaCal = nBetaGammaCal+1
        tempParSubset(nBetaGammaCal)%ixMaster = ivar
        tempParSubset(nBetaGammaCal)%pname    = betaMaster(ivar)%pname
        tempParSubset(nBetaGammaCal)%val      = betaMaster(ivar)%val
        tempParSubset(nBetaGammaCal)%lwr      = betaMaster(ivar)%lwr
        tempParSubset(nBetaGammaCal)%upr      = betaMaster(ivar)%upr
        tempParSubset(nBetaGammaCal)%beta     = betaMaster(ivar)%beta
        tempParSubset(nBetaGammaCal)%tftype   = betaMaster(ivar)%tftype
        tempParSubset(nBetaGammaCal)%ptype    = betaMaster(ivar)%ptype
        tempParSubset(nBetaGammaCal)%flag     = .True. 
        tempParSubset(nBetaGammaCal)%hups     = betaMaster(ivar)%hups
        tempParSubset(nBetaGammaCal)%hpnorm   = betaMaster(ivar)%hpnorm
        tempParSubset(nBetaGammaCal)%vups     = betaMaster(ivar)%vups
        tempParSubset(nBetaGammaCal)%vpnorm   = betaMaster(ivar)%vpnorm
        tempParSubset(nBetaGammaCal)%perLyr   = betaMaster(ivar)%perLyr
      endif
    enddo
    nBetaDirCal=nBetaGammaCal-nGammaCal      ! Save number of beta parameter to be directly calibrated
    if (nBetaGammaCal > 0_i2b) then          ! Save 'parSubset'
      allocate(parSubset(nBetaGammaCal))
      parSubset=tempParSubset(1:nBetaGammaCal) 
    endif
    if (nGammaCal > 0_i2b) then              ! Save'gammaSubset'
      allocate(gammaSubset(nGammaCal))
      gammaSubset=tempGammaMeta(1:nGammaCal) 
    endif
    return
  end subroutine

  ! Private subroutine: Get list of beta parameters (excluding z and h) computed with MPR.  
  ! populated global data: betaInGamma, soilBetaInGamma, vegBetaInGamma
  subroutine get_betaInGamma( err, message )
    !input
    !output
    character(len=strLen),intent(out)    :: message                ! error message for current routine
    integer(i4b),         intent(out)    :: err                    ! error code
    !local
    logical(lgc),allocatable             :: mask(:)
    integer(i4b)                         :: iPar                   ! loop index of parameter 
    integer(i4b)                         :: ivar                   ! index in master parameter
    integer(i4b)                         :: iSoil                  ! counter for soil beta parameters to be computed with MPR 
    integer(i4b)                         :: iVeg                   ! counter for vege beta parameters to be computed with MPR
    character(len=strLen),allocatable    :: allbeta(:)             ! 
    character(len=strLen),allocatable    :: tempSoilBetaInGamma(:) ! 
    character(len=strLen),allocatable    :: tempVegBetaInGamma(:)  ! 
    character(len=strLen),allocatable    :: res(:)                 ! 
    integer(i4b)                         :: i,j,k                  ! loop index

    err=0; message="get_betaInGamma/"
    if ( allocated(gammaSubset) ) then ! get beta parameter associated with gamma parameters (excluding soil depth, layer thickness)
      allocate(res(size(gammaSubset)))
      k = 1
      res(1) = gammaSubset(1)%beta
      outer:do i=2,size(gammaSubset)
        inner:do j=1,k
          ! if find a match so start looking again
          if (res(j)==gammaSubset(i)%beta)then; cycle outer; endif
        end do inner
        ! No match found so add it to the output
        k = k + 1
        res(k) = gammaSubset(i)%beta
      end do outer
      allocate(allbeta(k))
      allbeta=res(1:k) 
      ! exclude h and z parameters from allbeta to get betaInGamma
      allocate(mask(k))
      mask=.true.
      do i=1,k
        if ( trim(allbeta(i)) == 'h1' ) mask(i)=.false. 
        if ( trim(allbeta(i)) == 'h2' ) mask(i)=.false. 
        if ( trim(allbeta(i)) == 'h3' ) mask(i)=.false. 
        if ( trim(allbeta(i)) == 'h4' ) mask(i)=.false. 
        if ( trim(allbeta(i)) == 'h5' ) mask(i)=.false. 
        if ( trim(allbeta(i)) == 'z' )  mask(i)=.false. 
      end do
      allocate(betaInGamma(count(mask)))
      betaInGamma=pack(allbeta,mask)
      ! Count number of soil and Vege parameters to be computed with MPR 
      nSoilParModel=0
      nVegParModel=0
      do iPar=1,size(betaInGamma)
        ivar=get_ixBeta(betaInGamma(iPar))
        if (betaMaster(ivar)%ptype=='soil') nSoilParModel = nSoilParModel+1
        if (betaMaster(ivar)%ptype=='veg')  nVegParModel  = nVegParModel +1
      enddo
      iSoil=0;iVeg=0
      allocate(tempSoilBetaInGamma(nBeta),stat=err)
      allocate(tempVegBetaInGamma(nBeta),stat=err)
      do iPar=1,size(betaInGamma)
        ivar = get_ixBeta(betaInGamma(iPar))
        if(ivar<=0)then; err=40; message=trim(message)//"2.variableNotFound[var="//trim(betaInGamma(iPar))//"]"; return; endif
        if (betaMaster(ivar)%ptype=='soil')then; iSoil=iSoil+1; tempSoilBetaInGamma(iSoil) = betaMaster(ivar)%pname; endif
        if (betaMaster(ivar)%ptype=='veg') then; iVeg=iVeg+1;   tempVegBetaInGamma(iVeg)   = betaMaster(ivar)%pname; endif
      enddo 
      allocate(soilBetaInGamma(nSoilParModel))
      soilBetaInGamma=tempSoilBetaInGamma(1:nSoilParModel) 
      allocate(vegBetaInGamma(nVegParModel))
      vegBetaInGamma=tempVegBetaInGamma(1:nVegParModel) 
    else
      print*, 'NO gamma parameters included in the list'
    endif
    return
  end subroutine

  ! Private subroutine: to populate scaling parameter meta for beta parameter calibrated with MPR
  ! populated global data: betaCalScale
  subroutine get_betaCalScale( err, message ) 
    !input
    !output
    character(len=strLen),intent(out)   :: message                ! error message for current routine
    integer(i4b),         intent(out)   :: err                    ! error code
    !local
    type(scale_meta),     allocatable    :: tempBetaCalScale(:)   ! 
    integer(i4b)                         :: nScaleBeta            ! counter 

    err=0; message="get_scaleInBeta/"
    allocate(tempBetaCalScale(size(calParMeta)),stat=err);if(err/=0)then;message=trim(message)//'error allocating tempBetaCalScale';return;endif
    nScaleBeta=0 
    do iPar=1,size(calParMeta)!if beta parameter is estimated with MPR or even not calibrated, calibrate scaling parameter
      if ( calParMeta(iPar)%calMethod==1 .or. calParMeta(iPar)%calMethod==2 ) then  
        associate( hups      => betaMaster(get_ixBeta(calParMeta(iPar)%betaname))%hups,   &
                   vups      => betaMaster(get_ixBeta(calParMeta(iPar)%betaname))%vups,   &
                   hpower    => betaMaster(get_ixBeta(calParMeta(iPar)%betaname))%hpnorm, &
                   vpower    => betaMaster(get_ixBeta(calParMeta(iPar)%betaname))%vpnorm, &
                   hscaleMask=> calParMeta(iPar)%isScaleCalH, &
                   vscaleMask=> calParMeta(iPar)%isScaleCalV )
        nScaleBeta=nScaleBeta+1
        tempBetaCalScale(nScaleBeta)%betaname    = calParMeta(iPar)%betaname  
        tempBetaCalScale(nScaleBeta)%pdefault(1) = hpower 
        tempBetaCalScale(nScaleBeta)%pdefault(2) = vpower 
        tempBetaCalScale(nScaleBeta)%mask(1)     = hscaleMask 
        tempBetaCalScale(nScaleBeta)%mask(2)     = vscaleMask 
        ! Check beta meta data for scaling and if value is -999 or na, turns off calibration 
        if (hscaleMask .and. (hpower==-999.0_dp .or. hups=='na'))then 
          print*,'Switch horizontal upscale calibration to False for',calParMeta(iPar)%betaname
          tempBetaCalScale(nScaleBeta)%mask(1)=.False.
        endif
        if (vscaleMask .and. (vpower==-999.0_dp .or. vups=='na'))then
          print*,'Switch vertical upscale calibration to False for',calParMeta(iPar)%betaname
          tempBetaCalScale(nScaleBeta)%mask(2)=.False.
        endif
        end associate
      endif
    enddo
    allocate(betaCalScale(nScaleBeta),stat=err);if(err/=0)then;message=trim(message)//'error allocating betaCalScale';return;endif
    betaCalScale=tempBetaCalScale(1:nScaleBeta)
    return
  end subroutine
  
  ! Private subroutine: to count calibrating parameters including:
  ! beta parameters (direct calibration) per layer, gamma parameters, scaling parameters (V and H) 
  ! save nParCalSum in public_var
  subroutine total_calParam( err, message )
    implicit none
    !output variables
    integer(i4b),         intent(out) :: err         ! error code
    character(*),         intent(out) :: message     ! error message
    ! local variables
    integer(i4b)                 :: iPar          ! loop indices
    
    ! initialize error control
    err=0; message='total_calParam/'
    nParCalSum=0_i4b
    do iPar=1,nBetaGammaCal        ! cound number of calibration beta (per layer) and gamma
      if (parSubset(iPar)%perLyr)then
        nParCalSum=nParCalSum+nLyr
      else
        nParCalSum=nParCalSum+1
      endif
    enddo  
    do iPar=1,size(calParMeta) ! Add number of calibration pnorm value - make False if value is -999 or nan
      if (calParMeta(iPar)%calMethod==1 .or. calParMeta(iPar)%calMethod==2 ) nParCalSum=nParCalSum+2
    enddo
    print*, nParCalSum
    return
  end subroutine

end subroutine

! ************************************************************************************************
! Public subroutine: convert parameter data structure to simple arrays 
! ************************************************************************************************
subroutine param_setup( err, message )
  use globalData,  only:parArray, parMask, parSubset, nBetaGammaCal, betaCalScale, nParCalSum 
  implicit none
  !output variables
  integer(i4b),     intent(out) :: err                    ! error code
  character(*),     intent(out) :: message                ! error message
  ! local variables
  integer(i2b)                  :: iPar                   ! loop indices
  integer(i2b)                  :: idx                    ! count of calibrating parameter including per layer parameter 
  integer(i2b)                  :: ixHV                   ! count of calibrating parameter including per layer parameter 
  
  ! initialize error control
  err=0; message='param_setput/'
  allocate(parArray(nParCalSum,3),stat=err);if(err/=0)then;message=trim(message)//'error allocating parArray';return;endif
  allocate(parMask(nParCalSum),stat=err);if(err/=0)then;message=trim(message)//'error allocating parMask';return;endif
  idx=0_i2b
  do iPar=1,nBetaGammaCal
    if (parSubset(iPar)%perLyr)then
      idx=idx+nLyr
      parArray(idx-nLyr+1:idx,1) = parSubset(iPar)%val
      parArray(idx-nLyr+1:idx,2) = parSubset(iPar)%lwr
      parArray(idx-nLyr+1:idx,3) = parSubset(iPar)%upr
      parMask (idx-nLyr+1:idx)   = parSubset(iPar)%flag
    else
      idx=idx+1
      parArray(idx,1) = parSubset(iPar)%val
      parArray(idx,2) = parSubset(iPar)%lwr
      parArray(idx,3) = parSubset(iPar)%upr
      parMask (idx)   = parSubset(iPar)%flag
    endif
  enddo  
  do iPar=1,size(betaCalScale) 
    do ixHV=1,2
      idx=idx+1
      parArray(idx,1) = betaCalScale(iPar)%pdefault(ixHV)        ! default value of pnorm value
      parArray(idx,2) = -100.0_dp                                ! lower bound or pnorm value
      parArray(idx,3) =  100.0_dp                                ! upper bound of pnorm value
      parMask (idx)   = betaCalScale(iPar)%mask(ixHV)
    enddo
  enddo
  return
end subroutine

!**********************************
! Not Used---  Public subroutine: check if h parameters exist in gamma parameter
!**********************************
subroutine check_gammaH( err, message)
  use globalData,   only: gammaSubset
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
  if ( allocated(gammaSubset) )then
    !check h parameters - now can chcek up to 5 layers
    do i=1,size(gammaSubset)
      if (gammaSubset(i)%pname=="h1gamma1")then;id(1)=1;cycle;endif 
      if (gammaSubset(i)%pname=="h2gamma1")then;id(2)=1;cycle;endif
      if (gammaSubset(i)%pname=="h3gamma1")then;id(3)=1;cycle;endif
      if (gammaSubset(i)%pname=="h4gamma1")then;id(4)=1;cycle;endif
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
