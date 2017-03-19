module process_meta 

  use nrtype
  use public_var

  implicit none

  private

  public::read_calPar
  public::get_parm_meta
  public::betaCollection
  public::total_calParam
  public::param_setup 
  public::check_gammaZ
  public::check_gammaH

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
  use data_type,  only:par_meta,cpar_meta
  use globalData, only:calParMeta,     & ! meta for beta parameter listed in 'calPar' input
                       gammaMaster,    & ! meta for all gamma parameters
                       betaMaster,     & ! meta for all beta parameters
                       parSubset,      & ! meta for only parameter listed in input
                       gammaSubset,    & ! meta for only gamma parameters listed in input
                       betaInGamma,    & ! list of beta parameters calibrated with MPR 
                       betaCalScale,   & ! meta for beta parameter whose scaling operator(s) is calibrated
                       nBetaGammaCal,  & ! sum of beta and gamma parameters to be calibrated 
                       nBetaCal,       & ! number of beta parameters to be calibrated directly
                       nGammaCal,      & ! number of gamma parameters to be calibrated 
                       nSoilParModel,  & ! number of soil beta parameters to be calibrated 
                       nVegParModel,   & ! number of veg beta parameters to be calibrated 
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
  call get_parMetaSubset      ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_betaInGamma        ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_betaCalScale       ( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif

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
  ! Saved data:  parSubset, gammaSubset, nBetaGammaCal, nBetaCal, nGammaCal
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
      if(calParMeta(iBeta)%calMethod==1) then !if calMethod is MPR
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
            tempParSubset(nBetaGammaCal)%flag     = .True.
            tempParSubset(nBetaGammaCal)%hups     = gammaMaster(iGamma)%hups
            tempParSubset(nBetaGammaCal)%hpnorm   = gammaMaster(iGamma)%hpnorm
            tempParSubset(nBetaGammaCal)%vups     = gammaMaster(iGamma)%vups
            tempParSubset(nBetaGammaCal)%vpnorm   = gammaMaster(iGamma)%vpnorm
            tempParSubset(nBetaGammaCal)%perLyr   = gammaMaster(iGamma)%perLyr
            nGammaCal = nGammaCal+1
            tempGammaMeta(nGammaCal) = tempParSubset(nBetaGammaCal)
          endif
        end do
      elseif (calParMeta(iBeta)%calMethod==2) then ! if calMethod is direct calibration
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
    nBetaCal=nBetaGammaCal-nGammaCal      ! Save number of beta parameter to be directly calibrated
    if (nBetaGammaCal > 0_i2b) then ! Save 'parSubset'
      allocate(parSubset(nBetaGammaCal))
      parSubset=tempParSubset(1:nBetaGammaCal) 
    endif
    if (nGammaCal > 0_i2b) then     ! Save'gammaSubset'
      allocate(gammaSubset(nGammaCal))
      gammaSubset=tempGammaMeta(1:nGammaCal) 
    endif
    return
  end subroutine

  ! Private subroutine:   
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

  ! Private subroutine:   
  subroutine get_betaCalScale( err, message ) 
    use data_type,  only:scale_meta
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
    do iPar=1,size(calParMeta)
      if (calParMeta(iPar)%calMethod==1) then  ! if beta parameter is estimated with MPR, potentially calibrate pnorm value
        associate( ix=> get_ixBeta(calParMeta(iPar)%betaname))
        nScaleBeta=nScaleBeta+1
        tempBetaCalScale(nScaleBeta)%betaname    = calParMeta(iPar)%betaname  
        tempBetaCalScale(nScaleBeta)%pdefault(1) = betaMaster(ix)%hpnorm
        tempBetaCalScale(nScaleBeta)%pdefault(2) = betaMaster(ix)%vpnorm
        tempBetaCalScale(nScaleBeta)%mask(1)     = calParMeta(iPar)%isScaleCalH
        tempBetaCalScale(nScaleBeta)%mask(2)     = calParMeta(iPar)%isScaleCalV
        end associate
      endif
    enddo
    allocate(betaCalScale(nScaleBeta),stat=err);if(err/=0)then;message=trim(message)//'error allocating betaCalScale';return;endif
    betaCalScale=tempBetaCalScale(1:nScaleBeta)
    return
  end subroutine

end subroutine
  
! ********************************************************************************************
! Public subroutine: collect beta parameters needed for parameter estimation
! *********************************************************************************************
! Giving gamma parameters listed in "calPar" to compute associated beta parameters, figure out dependent parameters and make a list of all the gamma parameters  
! cluding gamma parameters in "calPar" and dependent gamma parameters
subroutine betaCollection(err,message)  
  use globalData, only:betaMaster, betaInGamma, betaNeeded
  use get_ixname, only:get_ixBeta
  use var_lookup, only:ixBeta, nBeta
  implicit none
  ! output
  character(len=strLen),intent(out)   :: message                ! error message for current routine
  integer(i4b),         intent(out)   :: err                   ! error code
  ! local variables
  character(len=strLen),allocatable   :: betaTemp(:)            ! temporal holder for name of beta parameters to be estimated and dependent beta parameters
  integer(i4b)                        :: iParm                  ! Loop index of model parameters 
  integer(i4b)                        :: nCheckDone             ! number of checked Beta parameters 
  logical(lgc)                        :: checkDone(nBeta)       ! used to check if the parameter is processed
  
  err=0; message="betaCollection/"
  allocate(betaTemp(nBeta))
  allocate(betaNeeded,source=betaInGamma)
  checkDone=.false.
  nCheckDone=count(checkDone)
  do
    do iParm = 1,size(betaNeeded)
      select case(get_ixBeta(betaNeeded(iParm)))
        case(ixBeta%ks);                   betaTemp(ixBeta%ks)       = betaMaster(ixBeta%ks)%pname;       checkDone(ixBeta%ks)=.true.
        case(ixBeta%bd);                   betaTemp(ixBeta%bd)       = betaMaster(ixBeta%bd)%pname;       checkDone(ixBeta%bd)=.true.
        case(ixBeta%phi)                   
          if(.not.checkDone(ixBeta%bd))    betaTemp(ixBeta%bd)       = betaMaster(ixBeta%bd)%pname;       checkDone(ixBeta%bd)=.true.
                                          betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
        case(ixBeta%b);                    betaTemp(ixBeta%b)        = betaMaster(ixBeta%b)%pname;        checkDone(ixBeta%b)=.true.
        case(ixBeta%psis);                 betaTemp(ixBeta%psis)     = betaMaster(ixBeta%psis)%pname;     checkDone(ixBeta%psis)=.true.
        case(ixBeta%fc)
          if(.not.checkDone(ixBeta%psis))  betaTemp(ixBeta%psis)     = betaMaster(ixBeta%psis)%pname;     checkDone(ixBeta%psis)=.true.
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%b))     betaTemp(ixBeta%b)        = betaMaster(ixBeta%b)%pname;        checkDone(ixBeta%b)=.true.
                                          betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
        case(ixBeta%wp)
          if(.not.checkDone(ixBeta%psis))  betaTemp(ixBeta%psis)     = betaMaster(ixBeta%psis)%pname;     checkDone(ixBeta%psis)=.true.
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%b))     betaTemp(ixBeta%b)        = betaMaster(ixBeta%b)%pname;        checkDone(ixBeta%b)=.true.
                                          betaTemp(ixBeta%wp)       = betaMaster(ixBeta%wp)%pname;       checkDone(ixBeta%wp)=.true.
        case(ixBeta%myu)
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
                                          betaTemp(ixBeta%myu)      = betaMaster(ixBeta%myu)%pname;      checkDone(ixBeta%myu)=.true.
        case(ixBeta%binfilt);              betaTemp(ixBeta%binfilt)  = betaMaster(ixBeta%binfilt)%pname;  checkDone(ixBeta%binfilt)=.true.
        case(ixBeta%D1)
          if(.not.checkDone(ixBeta%ks))    betaTemp(ixBeta%ks)       = betaMaster(ixBeta%ks)%pname;       checkDone(ixBeta%ks)=.true.
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
                                          betaTemp(ixBeta%D1)       = betaMaster(ixBeta%D1)%pname;       checkDone(ixBeta%D1)=.true.
        case(ixBeta%Ds)
          if(.not.checkDone(ixBeta%D1))    betaTemp(ixBeta%D1)       = betaMaster(ixBeta%D1)%pname;       checkDone(ixBeta%D1)=.true.
          if(.not.checkDone(ixBeta%D3))    betaTemp(ixBeta%D3)       = betaMaster(ixBeta%D3)%pname;       checkDone(ixBeta%D3)=.true.
          if(.not.checkDone(ixBeta%Dsmax)) betaTemp(ixBeta%Dsmax)    = betaMaster(ixBeta%Dsmax)%pname;    checkDone(ixBeta%Dsmax)=.true.
                                          betaTemp(ixBeta%Ds)       = betaMaster(ixBeta%Ds)%pname;       checkDone(ixBeta%Ds)=.true.
        case(ixBeta%D4);                   betaTemp(ixBeta%D4)       = betaMaster(ixBeta%D4)%pname;       checkDone(ixBeta%D4)=.true.
        case(ixBeta%c)
          if(.not.checkDone(ixBeta%D4))    betaTemp(ixBeta%D4)       = betaMaster(ixBeta%D4)%pname;       checkDone(ixBeta%D4)=.true.
                                          betaTemp(ixBeta%c)        = betaMaster(ixBeta%c)%pname;        checkDone(ixBeta%c)=.true.
        case(ixBeta%SD);                   betaTemp(ixBeta%SD)       = betaMaster(ixBeta%SD)%pname;       checkDone(ixBeta%SD)=.true.
        case(ixBeta%expt)                 
          if(.not.checkDone(ixBeta%b))     betaTemp(ixBeta%b)        = betaMaster(ixBeta%b)%pname;        checkDone(ixBeta%b)=.true.
                                          betaTemp(ixBeta%expt)     = betaMaster(ixBeta%expt)%pname;     checkDone(ixBeta%expt)=.true.
        case(ixBeta%D2)
          if(.not.checkDone(ixBeta%ks))    betaTemp(ixBeta%ks)       = betaMaster(ixBeta%ks)%pname;       checkDone(ixBeta%ks)=.true.
          if(.not.checkDone(ixBeta%D4))    betaTemp(ixBeta%D4)       = betaMaster(ixBeta%D4)%pname;       checkDone(ixBeta%D4)=.true.
                                          betaTemp(ixBeta%D2)       = betaMaster(ixBeta%D2)%pname;       checkDone(ixBeta%D2)=.true.
        case(ixBeta%Dsmax)
          if(.not.checkDone(ixBeta%D1))    betaTemp(ixBeta%D1)       = betaMaster(ixBeta%D1)%pname;       checkDone(ixBeta%D1)=.true.
          if(.not.checkDone(ixBeta%D2))    betaTemp(ixBeta%D2)       = betaMaster(ixBeta%D2)%pname;       checkDone(ixBeta%D2)=.true.
          if(.not.checkDone(ixBeta%D3))    betaTemp(ixBeta%D3)       = betaMaster(ixBeta%D3)%pname;       checkDone(ixBeta%D3)=.true.
          if(.not.checkDone(ixBeta%c))     betaTemp(ixBeta%c)        = betaMaster(ixBeta%c)%pname;        checkDone(ixBeta%c)=.true.
                                          betaTemp(ixBeta%Dsmax)    = betaMaster(ixBeta%Dsmax)%pname;    checkDone(ixBeta%Dsmax)=.true.
        case(ixBeta%bbl)
          if(.not.checkDone(ixBeta%expt))  betaTemp(ixBeta%expt)     = betaMaster(ixBeta%expt)%pname;     checkDone(ixBeta%expt)=.true.
                                          betaTemp(ixBeta%bbl)      = betaMaster(ixBeta%bbl)%pname;      checkDone(ixBeta%bbl)=.true.
        case(ixBeta%WcrFrac)
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
                                          betaTemp(ixBeta%WcrFrac)  = betaMaster(ixBeta%WcrFrac)%pname;  checkDone(ixBeta%WcrFrac)=.true.
        case(ixBeta%WpwpFrac)
          if(.not.checkDone(ixBeta%wp))    betaTemp(ixBeta%wp)       = betaMaster(ixBeta%wp)%pname;       checkDone(ixBeta%wp)=.true.
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
                                          betaTemp(ixBeta%WpwpFrac) = betaMaster(ixBeta%WpwpFrac)%pname; checkDone(ixBeta%WpwpFrac)=.true.
        case(ixBeta%D3)
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
                                          betaTemp(ixBeta%D3)       = betaMaster(ixBeta%D3)%pname;       checkDone(ixBeta%D3)=.true.
        case(ixBeta%Ws)
          if(.not.checkDone(ixBeta%D3))    betaTemp(ixBeta%D3)       = betaMaster(ixBeta%D3)%pname;       checkDone(ixBeta%D3)=.true.
                                          betaTemp(ixBeta%Ws)       = betaMaster(ixBeta%Ws)%pname;       checkDone(ixBeta%Ws)=.true.
        case(ixBeta%twm)                  
          if(.not.checkDone(ixBeta%wp))    betaTemp(ixBeta%wp)       = betaMaster(ixBeta%wp)%pname;       checkDone(ixBeta%wp)=.true.
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
                                          betaTemp(ixBeta%twm)      = betaMaster(ixBeta%twm)%pname;      checkDone(ixBeta%twm)=.true.
        case(ixBeta%fwm)                  
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
                                          betaTemp(ixBeta%fwm)      = betaMaster(ixBeta%fwm)%pname;      checkDone(ixBeta%fwm)=.true.
        case(ixBeta%fsm)
          if(.not.checkDone(ixBeta%fwm))   betaTemp(ixBeta%fwm)      = betaMaster(ixBeta%fwm)%pname;      checkDone(ixBeta%fwm)=.true.
                                          betaTemp(ixBeta%fsm)      = betaMaster(ixBeta%fsm)%pname;      checkDone(ixBeta%fsm)=.true.
        case(ixBeta%fpm)
          if(.not.checkDone(ixBeta%fwm))   betaTemp(ixBeta%fwm)      = betaMaster(ixBeta%fwm)%pname;      checkDone(ixBeta%fwm)=.true.
          if(.not.checkDone(ixBeta%fsm))   betaTemp(ixBeta%fsm)      = betaMaster(ixBeta%fsm)%pname;      checkDone(ixBeta%fsm)=.true.
                                          betaTemp(ixBeta%fpm)      = betaMaster(ixBeta%fpm)%pname;      checkDone(ixBeta%fpm)=.true.
        case(ixBeta%zk)                   
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
                                          betaTemp(ixBeta%zk)       = betaMaster(ixBeta%zk)%pname;       checkDone(ixBeta%zk)=.true.
        case(ixBeta%zsk)
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%fc))    betaTemp(ixBeta%fc)       = betaMaster(ixBeta%fc)%pname;       checkDone(ixBeta%fc)=.true.
          if(.not.checkDone(ixBeta%wp))    betaTemp(ixBeta%wp)       = betaMaster(ixBeta%wp)%pname;       checkDone(ixBeta%wp)=.true.
                                          betaTemp(ixBeta%zsk)      = betaMaster(ixBeta%zsk)%pname;      checkDone(ixBeta%zsk)=.true.
        case(ixBeta%zpk)
          if(.not.checkDone(ixBeta%ks))    betaTemp(ixBeta%ks)       = betaMaster(ixBeta%ks)%pname;       checkDone(ixBeta%ks)=.true.
          if(.not.checkDone(ixBeta%myu))   betaTemp(ixBeta%myu)      = betaMaster(ixBeta%myu)%pname;      checkDone(ixBeta%myu)=.true.
                                          betaTemp(ixBeta%zpk)      = betaMaster(ixBeta%zpk)%pname;      checkDone(ixBeta%zpk)=.true.
        case(ixBeta%pfree)
          if(.not.checkDone(ixBeta%phi))   betaTemp(ixBeta%phi)      = betaMaster(ixBeta%phi)%pname;      checkDone(ixBeta%phi)=.true.
          if(.not.checkDone(ixBeta%wp))    betaTemp(ixBeta%wp)       = betaMaster(ixBeta%wp)%pname;       checkDone(ixBeta%wp)=.true.
                                          betaTemp(ixBeta%pfree)    = betaMaster(ixBeta%pfree)%pname;    checkDone(ixBeta%pfree)=.true.
        case(ixBeta%zperc)
          if(.not.checkDone(ixBeta%twm))   betaTemp(ixBeta%twm)      = betaMaster(ixBeta%twm)%pname;      checkDone(ixBeta%twm)=.true.
          if(.not.checkDone(ixBeta%fsm))   betaTemp(ixBeta%fsm)      = betaMaster(ixBeta%fsm)%pname;      checkDone(ixBeta%fsm)=.true.
          if(.not.checkDone(ixBeta%zsk))   betaTemp(ixBeta%zsk)      = betaMaster(ixBeta%zsk)%pname;      checkDone(ixBeta%zsk)=.true.
          if(.not.checkDone(ixBeta%fpm))   betaTemp(ixBeta%fpm)      = betaMaster(ixBeta%fpm)%pname;      checkDone(ixBeta%fpm)=.true.
          if(.not.checkDone(ixBeta%zpk))   betaTemp(ixBeta%zpk)      = betaMaster(ixBeta%zpk)%pname;      checkDone(ixBeta%zpk)=.true.
                                          betaTemp(ixBeta%zperc)    = betaMaster(ixBeta%zperc)%pname;    checkDone(ixBeta%zperc)=.true.
        case(ixBeta%rexp);                 
          if(.not.checkDone(ixBeta%psis))  betaTemp(ixBeta%wp)       = betaMaster(ixBeta%wp)%pname;       checkDone(ixBeta%wp)=.true.
                                          betaTemp(ixBeta%rexp)     = betaMaster(ixBeta%rexp)%pname;     checkDone(ixBeta%rexp)=.true.
        case(ixBeta%lai);                  betaTemp(ixBeta%lai)      = betaMaster(ixBeta%lai)%pname;      checkDone(ixBeta%lai)=.true.
      end select ! end of parameter case
    end do ! end of parameter loop
    if (nCheckDone==count(checkDone)) exit
    nCheckDone=count(checkDone)
    deallocate(betaNeeded,stat=err); if(err/=0)then;message=trim(message)//'error deallocating betaNeeded';return;endif
    allocate(betaNeeded(nCheckDone))
    betaNeeded=pack(betaTemp,checkDone)
  end do ! end of infinite loop

end subroutine

!**********************************
! Public subroutine: check if z parameters exist in gamma parameter
!**********************************
subroutine check_gammaZ( err, message)
  use globalData,   only: gammaSubset
  implicit none
  !output variables
  integer(i4b),                      intent(out)   :: err         ! error code
  character(*),                      intent(out)   :: message     ! error message
  !local variables
  logical(lgc)                                     :: checkZ 
  integer(i4b)                                     :: i           ! loop index

  ! initialize error control
  err=0; message='check_gammaZ/'
  checkZ=.false.
  if ( allocated(gammaSubset) )then
    !check z parameter
    do i=1,size(gammaSubset)
      if (gammaSubset(i)%pname=="z1gamma1") checkZ=.true.
    end do
    if ( .not. checkZ )then;err=40;message=trim(message)//"Calibrating gamma parameter require z1gamma1 parameter"; return;endif
  else
    print*, 'No gamma parameters listed in CalPar->No MPR computation'
  endif
  return
end subroutine

!**********************************
! Public subroutine: check if h parameters exist in gamma parameter
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

! ************************************************************************************************
! Public subroutine: Count total number of calibrating parameters 
! ************************************************************************************************
! Calibrating parameters including:
! beta parameters (direct calibration) per layer, gamma parameters, pnorm parameters (V and H) 
subroutine total_calParam( )
  ! save nParCalSum in public_var
  use globalData,  only: parSubset, nBetaGammaCal, calParMeta
  implicit none
  ! local variables
  integer(i4b)                 :: iPar          ! loop indices
  
  nParCalSum=0_i4b
  do iPar=1,nBetaGammaCal        ! cound number of calibration beta (per layer) and gamma
    if (parSubset(iPar)%perLyr)then
      nParCalSum=nParCalSum+nLyr
    else
      nParCalSum=nParCalSum+1
    endif
  enddo  
  do iPar=1,size(calParMeta) ! Add number of calibration pnorm value
    if (calParMeta(iPar)%calMethod==1) nParCalSum=nParCalSum+2
  enddo
  return
end subroutine

! ************************************************************************************************
! Public subroutine: convert parameter data structure to simple arrays 
! ************************************************************************************************
subroutine param_setup( param, mask )
  use globalData,  only:parSubset, nBetaGammaCal, betaCalScale 
  implicit none
  ! output variables
  real(dp),dimension(:,:),   intent(out)     :: param 
  logical,dimension(:),      intent(out)     :: mask
  ! local variables
  integer(i2b)                               :: iPar    ! loop indices
  integer(i2b)                               :: idx     ! count of calibrating parameter including per layer parameter 
  integer(i2b)                               :: ixHV    ! count of calibrating parameter including per layer parameter 
  
  idx=0_i2b
  do iPar=1,nBetaGammaCal
    if (parSubset(iPar)%perLyr)then
      idx=idx+nLyr
      param(idx-nLyr+1:idx,1) = parSubset(iPar)%val
      param(idx-nLyr+1:idx,2) = parSubset(iPar)%lwr
      param(idx-nLyr+1:idx,3) = parSubset(iPar)%upr
      mask (idx-nLyr+1:idx)   = parSubset(iPar)%flag
    else
      idx=idx+1
      param(idx,1) = parSubset(iPar)%val
      param(idx,2) = parSubset(iPar)%lwr
      param(idx,3) = parSubset(iPar)%upr
      mask (idx)   = parSubset(iPar)%flag
    endif
  enddo  

  do iPar=1,size(betaCalScale) 
    do ixHV=1,2
      idx=idx+1
      param(idx,1)=betaCalScale(iPar)%pdefault(ixHV)        ! default value of pnorm value
      param(idx,2)=-100                                     ! lower bound or pnorm value
      param(idx,3)=100                                      ! upper bound of pnorm value
      mask(idx)=betaCalScale(iPar)%mask(ixHV)
    enddo
  enddo

  return
end subroutine

end module process_meta 
