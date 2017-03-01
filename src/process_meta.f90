module process_meta 

  use nrtype
  use public_var

  implicit none

  private

  public::read_parm_master_meta
  public::read_calPar
  public::get_parm_meta
  public::betaCollection
  public::total_calParam
  public::param_setup 
  public::check_gammaZ
  public::check_gammaH

contains

! ************************************************************************************************
! Subroutine: read soil parameter metadata from a meta file 
! ************************************************************************************************
subroutine read_parm_master_meta(infile, err, message)
  ! used to read metadata from an input file and populate the appropriate metadata structure
  use data_type,  only:par_meta            ! metadata structure
  use ascii_util, only:file_open
  use get_ixname, only:get_ixPar
  use globalData, only:parMaster

  implicit none

  ! define input
  character(*),intent(in)              :: infile         ! input filename
  ! define output
  integer(i4b),intent(out)             :: err            ! error code
  character(*),intent(out)             :: message        ! error message
  ! define local variables
  character(len=256)                   :: cmessage       ! error message for downwind routine
  integer(i4b)                         :: unt            ! DK: need to either define units globally, or use getSpareUnit
  integer(i4b)                         :: iline          ! loop through lines in the file 
  integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
  character(LEN=256)                   :: temp           ! single lime of information
  integer(i4b)                         :: iend           ! check for the end of the file
  character(LEN=256)                   :: ffmt           ! file format
  type(par_meta)                       :: parmdTemp      ! temporary metadata structure
  character(len=1)                     :: dLim(12)        ! column delimiter
  integer(i4b)                         :: ivar           ! index of model variable
  ! Start procedure here
  err=0; message="read_param_meta/"
  ! open file
  call file_open(trim(infile),unt,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! get to the start of the variable descriptions 
  do iline=1,maxLines
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
    if (temp(1:1)/='!') exit  ! assume first line not comment is format code
  end do ! looping through file to find the format code
  ! read in format string
  read(temp,*)ffmt
  ! loop through the lines in the file
  do iline=1,maxLines
    ! read a line of data and exit iif an error code (character read, so only possible error is end of file)
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
    ! check that the line is not a comment
    if (temp(1:1)=='!')cycle
    ! save data into a temporary structure
    read(temp,trim(ffmt),iostat=err) parmdTemp%pname,  dLim(1),&
                                     parmdTemp%val,    dLim(2),&
                                     parmdTemp%lwr,    dLim(3),&
                                     parmdTemp%upr,    dLim(4),& 
                                     parmdTemp%beta,   dLim(5),& 
                                     parmdTemp%tftype, dLim(6),& 
                                     parmdTemp%ptype,  dLim(7),&
                                     parmdTemp%flag,   dLim(8),&
                                     parmdTemp%hups,   dLim(9),&
                                     parmdTemp%hpnorm, dLim(10),&
                                     parmdTemp%vups,   dLim(11),&
                                     parmdTemp%vpnorm, dLim(12),&
                                     parmdTemp%perLyr
    if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
    ! check that the delimiters are in the correct place
    if(any(dLim /= '|'))then
     message=trim(message)//'delimiter is not in the correct place; line = ['//trim(temp)//']; filename = '//trim(infile)
     err=32; return
    endif
    ! identify the index of the named variable
    ivar = get_ixPar(parmdTemp%pname)
    if(ivar<=0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(parmdTemp%pname)//"]"; return; endif
    ! check if index is within range
    if(ivar>size(parMaster))then; err=50; message=trim(message)//"variableExceedsVectorSize[par="//trim(parmdTemp%pname)//"]"; return; endif
    ! put data into the metadata vector
    parMaster(ivar) = parmdTemp 
  enddo  ! looping through lines in the file

  ! check that all elements are populated
  if(any(parMaster(:)%pname==''))then
   do iline=1,size(parMaster)
    print*,iline,' -> ',trim(parMaster(iline)%pname)
   end do
   err=40; message=trim(message)//"someVariablesNotPopulated"; return
  endif
  ! close file unit
  close(unt)
 end subroutine 

! ************************************************************************************************
! Public subroutine: Read 'calPar' input file 
! ************************************************************************************************
subroutine read_calPar(infile, err, message)
  use globalData, only:calParMeta      ! meta for calPar input  
  use data_type,  only:input_meta
  use ascii_util, only:file_open
  use var_lookup, only:nPar
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
  allocate(tempCalParMeta(nPar))
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
  ! Process calParMeta along with parMaster (popMeta.f90) 
  ! Saved data:  gammaMaster, betaMaster, parSubset, gammaSubset, betaInGamma data structure 
  use data_type,  only:par_meta,cpar_meta
  use globalData, only:calParMeta,  & ! meta for beta parameter listed in 'calPar' input
                       parMaster,   & ! meta for all gamma and beta parameters
                       gammaMaster, & ! meta for all gamma parameter 
                       betaMaster,  & ! meta for all beta parameter 
                       parSubset,   & ! meta for only parameter listed in input
                       gammaSubset, & ! meta for only gamma parameters listed in input
                       betaInGamma, & ! list of beta parameter associated with gamma parameters in list
                       betaCalScale,& ! meta for beta parameter whose scaling operator(s) is calibrated
                       nBetaGamma,  & ! sum of beta and gamma parameters to be calibrated 
                       nBeta,       & ! number of beta parameters to be calibrated directly
                       nGamma         ! number of gamma parameters to be calibrated 
  use get_ixname, only:get_ixPar
  use var_lookup, only:nPar
  implicit none
  ! input
  ! output
  integer(i4b),intent(out)             :: err              ! error code
  character(*),intent(out)             :: message          ! error message
  ! local variables
  character(len=strLen)                :: cmessage         ! error message from subroutine
  type(cpar_meta),allocatable          :: tempGammaMeta(:)
  type(cpar_meta),allocatable          :: tempParSubset(:)
  logical(lgc),allocatable             :: mask(:)
  integer(i4b)                         :: iBeta            ! loop index of lines in calPar input file 
  integer(i4b)                         :: ivar             ! loop index of master parameter  
  integer(i4b)                         :: iPar             ! loop index of master parameter  
  integer(i4b)                         :: iGamma           ! loop index of all the gamma parameters in master 
  integer(i4b)                         :: iend             ! check for the end of the file
 
  err=0; message="get_param_meta/"
  ! update parMaster(:)%tftype (transfer function type)
  do iBeta=1,size(calParMeta)
    ivar=get_ixPar(calParMeta(iBeta)%betaname)
    if(ivar<=0)then; err=40; message=trim(message)//"1.variableNotFound[var="//trim(calParMeta(iBeta)%betaname)//"]"; return; endif
    parMaster(ivar)%tftype=calParMeta(iBeta)%TF
  enddo
  ! split parMaster into gammaMaster and betaMaster 
  ! gammaMaster
  allocate(mask(nPar),stat=err); if(err/=0)then;message=trim(message)//'1.error allocating mask';return;endif
  mask=parMaster(:)%beta/='beta'
  allocate(gammaMaster(count(mask)))
  gammaMaster=pack(parMaster,mask)
  deallocate(mask,stat=err); if(err/=0)then;message=trim(message)//'1.error deallocating mask';return;endif
  ! betaMaster
  allocate(mask(nPar),stat=err); if(err/=0)then;message=trim(message)//'2.error allocating mask';return;endif
  mask=parMaster(:)%beta=='beta'
  allocate(betaMaster(count(mask)))
  betaMaster=pack(parMaster,mask)
  deallocate(mask,stat=err); if(err/=0)then;message=trim(message)//'2.error deallocating mask';return;endif
  allocate(tempParSubset(nPar))
  allocate(tempGammaMeta(nPar))
  nBetaGamma=0 ! Count and Save number of parameters to be calibrated
  nGamma=0     ! Count and Save number of gamma parameter to be calibrated
  do iBeta=1,size(calParMeta)
    if(calParMeta(iBeta)%calMethod==1) then !if calMethod is MPR
      do iGamma=1,size(gammaMaster) ! look for gammma parameters associated with this beta parameter
        if (gammaMaster(igamma)%beta==calparmeta(ibeta)%betaname .and. gammaMaster(igamma)%tftype==calparmeta(ibeta)%TF) then
          nBetaGamma = nBetaGamma+1
          ivar=get_ixPar(gammaMaster(iGamma)%pname)
          tempParSubset(nBetaGamma)%ixMaster = ivar
          tempParSubset(nBetaGamma)%pname    = gammaMaster(iGamma)%pname
          tempParSubset(nBetaGamma)%val      = gammaMaster(iGamma)%val
          tempParSubset(nBetaGamma)%lwr      = gammaMaster(iGamma)%lwr
          tempParSubset(nBetaGamma)%upr      = gammaMaster(iGamma)%upr
          tempParSubset(nBetaGamma)%beta     = gammaMaster(iGamma)%beta
          tempParSubset(nBetaGamma)%tftype   = gammaMaster(iGamma)%tftype
          tempParSubset(nBetaGamma)%ptype    = gammaMaster(iGamma)%ptype
          tempParSubset(nBetaGamma)%flag     = .True.
          tempParSubset(nBetaGamma)%hups     = gammaMaster(iGamma)%hups
          tempParSubset(nBetaGamma)%hpnorm   = gammaMaster(iGamma)%hpnorm
          tempParSubset(nBetaGamma)%vups     = gammaMaster(iGamma)%vups
          tempParSubset(nBetaGamma)%vpnorm   = gammaMaster(iGamma)%vpnorm
          tempParSubset(nBetaGamma)%perLyr   = gammaMaster(iGamma)%perLyr
          nGamma = nGamma+1
          tempGammaMeta(nGamma) = tempParSubset(nBetaGamma)
        endif
      end do
    elseif (calParMeta(iBeta)%calMethod==2) then ! if calMethod is direct calibration
      ivar = get_ixPar(calParMeta(iBeta)%betaname)
      if(ivar<=0)then; err=40; message=trim(message)//"2.variableNotFound[var="//trim(calParMeta(iBeta)%betaname)//"]"; return; endif
      nBetaGamma = nBetaGamma+1
      tempParSubset(nBetaGamma)%ixMaster = ivar
      tempParSubset(nBetaGamma)%pname    = parMaster(ivar)%pname
      tempParSubset(nBetaGamma)%val      = parMaster(ivar)%val
      tempParSubset(nBetaGamma)%lwr      = parMaster(ivar)%lwr
      tempParSubset(nBetaGamma)%upr      = parMaster(ivar)%upr
      tempParSubset(nBetaGamma)%beta     = parMaster(ivar)%beta
      tempParSubset(nBetaGamma)%tftype   = parMaster(ivar)%tftype
      tempParSubset(nBetaGamma)%ptype    = parMaster(ivar)%ptype
      tempParSubset(nBetaGamma)%flag     = .True. 
      tempParSubset(nBetaGamma)%hups     = parMaster(ivar)%hups
      tempParSubset(nBetaGamma)%hpnorm   = parMaster(ivar)%hpnorm
      tempParSubset(nBetaGamma)%vups     = parMaster(ivar)%vups
      tempParSubset(nBetaGamma)%vpnorm   = parMaster(ivar)%vpnorm
      tempParSubset(nBetaGamma)%perLyr   = parMaster(ivar)%perLyr
    endif
  enddo
  nBeta=nBetaGamma-nGamma      ! Save number of beta parameter to be directly calibrated
  if (nBetaGamma > 0_i2b) then ! Save 'parSubset'
    allocate(parSubset(nBetaGamma))
    parSubset=tempParSubset(1:nBetaGamma) 
  endif
  if (nGamma > 0_i2b) then     ! Save'gammaSubset'
    allocate(gammaSubset(nGamma))
    gammaSubset=tempGammaMeta(1:nGamma) 
  endif
  call get_betaInGamma( err, cmessage ); if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  call get_betaCalScale( err, cmessage );if(err/=0)then;message=trim(message)//trim(cmessage); return; endif
  ! check that all elements are populated
  if(any(parSubset(:)%pname==''))then
    do iPar=1,size(parSubset)
      print*,iPar,' -> ',trim(parSubset(iPar)%pname)
    end do
    err=40; message=trim(message)//"parSubset is not filled out completely"; return
  endif
  print*, betaCalScale
  return

  contains

  ! Private subroutine:   
  subroutine get_betaInGamma( err, message )
    !input
    !output
    character(len=strLen),intent(out)   :: message           ! error message for current routine
    integer(i4b),         intent(out)   :: err               ! error code
    !local
    character(len=strLen),allocatable    :: allbeta(:)       ! 
    character(len=strLen),allocatable    :: res(:)           ! 
    integer(i4b)                         :: i,j,k            ! loop index
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
        associate( ix=> get_ixPar(calParMeta(iPar)%betaname))
        nScaleBeta=nScaleBeta+1
        tempBetaCalScale(nScaleBeta)%betaname    = calParMeta(iPar)%betaname  
        tempBetaCalScale(nScaleBeta)%pdefault(1) = parMaster(ix)%hpnorm
        tempBetaCalScale(nScaleBeta)%pdefault(2) = parMaster(ix)%vpnorm
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
  use globalData, only:parMaster, betaInGamma, betaNeeded
  use get_ixname, only:get_ixPar
  use var_lookup, only:ixPar, nPar
  implicit none
  ! output
  character(len=strLen),intent(out)   :: message                ! error message for current routine
  integer(i4b),         intent(out)   :: err                   ! error code
  ! local variables
  character(len=strLen),allocatable   :: betaTemp(:)            ! temporal holder for name of beta parameters to be estimated and dependent beta parameters
  integer(i4b)                        :: iParm                  ! Loop index of model parameters 
  integer(i4b)                        :: nCheckDone             ! number of checked Beta parameters 
  logical(lgc)                        :: checkDone(nPar)        ! used to check if the parameter is processed
  
  err=0; message="betaCollection/"
  allocate(betaTemp(nPar))
  allocate(betaNeeded,source=betaInGamma)
  checkDone=.false.
  nCheckDone=count(checkDone)
  do
    do iParm = 1,size(betaNeeded)
      select case(get_ixPar(betaNeeded(iParm)))
        case(ixPar%ks);                   betaTemp(ixPar%ks)       = parMaster(ixPar%ks)%pname;       checkDone(ixPar%ks)=.true.
        case(ixPar%bd);                   betaTemp(ixPar%bd)       = parMaster(ixPar%bd)%pname;       checkDone(ixPar%bd)=.true.
        case(ixPar%phi)                   
          if(.not.checkDone(ixPar%bd))    betaTemp(ixPar%bd)       = parMaster(ixPar%bd)%pname;       checkDone(ixPar%bd)=.true.
                                          betaTemp(ixPar%phi)      = parMaster(ixPar%phi)%pname;      checkDone(ixPar%phi)=.true.
        case(ixPar%b);                    betaTemp(ixPar%b)        = parMaster(ixPar%b)%pname;        checkDone(ixPar%b)=.true.
        case(ixPar%psis);                 betaTemp(ixPar%psis)     = parMaster(ixPar%psis)%pname;     checkDone(ixPar%psis)=.true.
        case(ixPar%fc)
          if(.not.checkDone(ixPar%psis))  betaTemp(ixPar%psis)     = parMaster(ixPar%psis)%pname;     checkDone(ixPar%psis)=.true.
          if(.not.checkDone(ixPar%phi))   betaTemp(ixPar%phi)      = parMaster(ixPar%phi)%pname;      checkDone(ixPar%phi)=.true.
          if(.not.checkDone(ixPar%b))     betaTemp(ixPar%b)        = parMaster(ixPar%b)%pname;        checkDone(ixPar%b)=.true.
                                          betaTemp(ixPar%fc)       = parMaster(ixPar%fc)%pname;       checkDone(ixPar%fc)=.true.
        case(ixPar%wp)
          if(.not.checkDone(ixPar%psis))  betaTemp(ixPar%psis)     = parMaster(ixPar%psis)%pname;     checkDone(ixPar%psis)=.true.
          if(.not.checkDone(ixPar%phi))   betaTemp(ixPar%phi)      = parMaster(ixPar%phi)%pname;      checkDone(ixPar%phi)=.true.
          if(.not.checkDone(ixPar%b))     betaTemp(ixPar%b)        = parMaster(ixPar%b)%pname;        checkDone(ixPar%b)=.true.
                                          betaTemp(ixPar%wp)       = parMaster(ixPar%wp)%pname;       checkDone(ixPar%wp)=.true.
        case(ixPar%myu)
          if(.not.checkDone(ixPar%phi))   betaTemp(ixPar%phi)      = parMaster(ixPar%phi)%pname;      checkDone(ixPar%phi)=.true.
          if(.not.checkDone(ixPar%fc))    betaTemp(ixPar%fc)       = parMaster(ixPar%fc)%pname;       checkDone(ixPar%fc)=.true.
                                          betaTemp(ixPar%myu)      = parMaster(ixPar%myu)%pname;      checkDone(ixPar%myu)=.true.
        case(ixPar%binfilt);              betaTemp(ixPar%binfilt)  = parMaster(ixPar%binfilt)%pname;  checkDone(ixPar%binfilt)=.true.
        case(ixPar%D1);
          if(.not.checkDone(ixPar%ks))    betaTemp(ixPar%ks)       = parMaster(ixPar%ks)%pname;       checkDone(ixPar%ks)=.true.
          if(.not.checkDone(ixPar%phi))   betaTemp(ixPar%phi)      = parMaster(ixPar%phi)%pname;      checkDone(ixPar%phi)=.true.
                                          betaTemp(ixPar%D1)       = parMaster(ixPar%D1)%pname;       checkDone(ixPar%D1)=.true.
        case(ixPar%Ds);
          if(.not.checkDone(ixPar%D1))    betaTemp(ixPar%D1)       = parMaster(ixPar%D1)%pname;       checkDone(ixPar%D1)=.true.
          if(.not.checkDone(ixPar%D3))    betaTemp(ixPar%D3)       = parMaster(ixPar%D3)%pname;       checkDone(ixPar%D3)=.true.
          if(.not.checkDone(ixPar%Dsmax)) betaTemp(ixPar%Dsmax)    = parMaster(ixPar%Dsmax)%pname;    checkDone(ixPar%Dsmax)=.true.
                                          betaTemp(ixPar%Ds)       = parMaster(ixPar%Ds)%pname;       checkDone(ixPar%Ds)=.true.
        case(ixPar%D4);                   betaTemp(ixPar%D4)       = parMaster(ixPar%D4)%pname;       checkDone(ixPar%D4)=.true.
        case(ixPar%c)
          if(.not.checkDone(ixPar%D4))    betaTemp(ixPar%D4)       = parMaster(ixPar%D4)%pname;       checkDone(ixPar%D4)=.true.
                                          betaTemp(ixPar%c)        = parMaster(ixPar%c)%pname;        checkDone(ixPar%c)=.true.
        case(ixPar%SD);                   betaTemp(ixPar%SD)       = parMaster(ixPar%SD)%pname;       checkDone(ixPar%SD)=.true.
        case(ixPar%expt);                 betaTemp(ixPar%expt)     = parMaster(ixPar%expt)%pname;     checkDone(ixPar%expt)=.true.
        case(ixPar%D2)
          if(.not.checkDone(ixPar%ks))    betaTemp(ixPar%ks)       = parMaster(ixPar%ks)%pname;       checkDone(ixPar%ks)=.true.
          if(.not.checkDone(ixPar%D4))    betaTemp(ixPar%D4)       = parMaster(ixPar%D4)%pname;       checkDone(ixPar%D4)=.true.
                                          betaTemp(ixPar%D2)       = parMaster(ixPar%D2)%pname;       checkDone(ixPar%D2)=.true.
        case(ixPar%Dsmax);
          if(.not.checkDone(ixPar%D1))    betaTemp(ixPar%D1)       = parMaster(ixPar%D1)%pname;       checkDone(ixPar%D1)=.true.
          if(.not.checkDone(ixPar%D2))    betaTemp(ixPar%D2)       = parMaster(ixPar%D2)%pname;       checkDone(ixPar%D2)=.true.
          if(.not.checkDone(ixPar%D3))    betaTemp(ixPar%D3)       = parMaster(ixPar%D3)%pname;       checkDone(ixPar%D3)=.true.
          if(.not.checkDone(ixPar%c))     betaTemp(ixPar%c)        = parMaster(ixPar%c)%pname;        checkDone(ixPar%c)=.true.
                                          betaTemp(ixPar%Dsmax)    = parMaster(ixPar%Dsmax)%pname;    checkDone(ixPar%Dsmax)=.true.
        case(ixPar%bbl);
          if(.not.checkDone(ixPar%expt))  betaTemp(ixPar%expt)     = parMaster(ixPar%expt)%pname;     checkDone(ixPar%expt)=.true.
                                          betaTemp(ixPar%bbl)      = parMaster(ixPar%bbl)%pname;      checkDone(ixPar%bbl)=.true.
        case(ixPar%WcrFrac)
          if(.not.checkDone(ixPar%D3))    betaTemp(ixPar%D3)       = parMaster(ixPar%D3)%pname;       checkDone(ixPar%D3)=.true.
                                          betaTemp(ixPar%WcrFrac)  = parMaster(ixPar%WcrFrac)%pname;  checkDone(ixPar%WcrFrac)=.true.
        case(ixPar%WpwpFrac)
          if(.not.checkDone(ixPar%D3))    betaTemp(ixPar%D3)       = parMaster(ixPar%D3)%pname;       checkDone(ixPar%D3)=.true.
                                          betaTemp(ixPar%WpwpFrac) = parMaster(ixPar%WpwpFrac)%pname; checkDone(ixPar%WpwpFrac)=.true.
        case(ixPar%D3)
          if(.not.checkDone(ixPar%fc))    betaTemp(ixPar%fc)       = parMaster(ixPar%fc)%pname;       checkDone(ixPar%fc)=.true.
                                          betaTemp(ixPar%D3)       = parMaster(ixPar%D3)%pname;       checkDone(ixPar%D3)=.true.
        case(ixPar%Ws)
          if(.not.checkDone(ixPar%D3))    betaTemp(ixPar%D3)       = parMaster(ixPar%D3)%pname;       checkDone(ixPar%D3)=.true.
                                          betaTemp(ixPar%Ws)       = parMaster(ixPar%Ws)%pname;       checkDone(ixPar%Ws)=.true.
        case(ixPar%twm);                  betaTemp(ixPar%twm)      = parMaster(ixPar%twm)%pname;      checkDone(ixPar%twm)=.true.
        case(ixPar%fwm);                  betaTemp(ixPar%fwm)      = parMaster(ixPar%fwm)%pname;      checkDone(ixPar%fwm)=.true.
        case(ixPar%fsm)
          if(.not.checkDone(ixPar%fwm))   betaTemp(ixPar%fwm)      = parMaster(ixPar%fwm)%pname;      checkDone(ixPar%fwm)=.true.
                                          betaTemp(ixPar%fsm)      = parMaster(ixPar%fsm)%pname;      checkDone(ixPar%fsm)=.true.
        case(ixPar%fpm)
          if(.not.checkDone(ixPar%fwm))   betaTemp(ixPar%fwm)      = parMaster(ixPar%fwm)%pname;      checkDone(ixPar%fwm)=.true.
          if(.not.checkDone(ixPar%fsm))   betaTemp(ixPar%fsm)      = parMaster(ixPar%fsm)%pname;      checkDone(ixPar%fsm)=.true.
                                          betaTemp(ixPar%fpm)      = parMaster(ixPar%fpm)%pname;      checkDone(ixPar%fpm)=.true.
        case(ixPar%zk);                   betaTemp(ixPar%zk)       = parMaster(ixPar%zk)%pname;       checkDone(ixPar%zk)=.true.
        case(ixPar%zsk);                  betaTemp(ixPar%zsk)      = parMaster(ixPar%zsk)%pname;      checkDone(ixPar%zsk)=.true.
        case(ixPar%zpk);                  betaTemp(ixPar%zpk)      = parMaster(ixPar%zpk)%pname;      checkDone(ixPar%zpk)=.true.
        case(ixPar%pfree);                betaTemp(ixPar%pfree)    = parMaster(ixPar%pfree)%pname;    checkDone(ixPar%pfree)=.true.
        case(ixPar%zperc)
          if(.not.checkDone(ixPar%twm))   betaTemp(ixPar%twm)      = parMaster(ixPar%twm)%pname;      checkDone(ixPar%twm)=.true.
          if(.not.checkDone(ixPar%fsm))   betaTemp(ixPar%fsm)      = parMaster(ixPar%fsm)%pname;      checkDone(ixPar%fsm)=.true.
          if(.not.checkDone(ixPar%zsk))   betaTemp(ixPar%zsk)      = parMaster(ixPar%zsk)%pname;      checkDone(ixPar%zsk)=.true.
          if(.not.checkDone(ixPar%fpm))   betaTemp(ixPar%fpm)      = parMaster(ixPar%fpm)%pname;      checkDone(ixPar%fpm)=.true.
          if(.not.checkDone(ixPar%zpk))   betaTemp(ixPar%zpk)      = parMaster(ixPar%zpk)%pname;      checkDone(ixPar%zpk)=.true.
                                          betaTemp(ixPar%zperc)    = parMaster(ixPar%zperc)%pname;    checkDone(ixPar%zperc)=.true.
        case(ixPar%rexp);                 betaTemp(ixPar%rexp)     = parMaster(ixPar%rexp)%pname;     checkDone(ixPar%rexp)=.true.
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
  use globalData,  only: parSubset, nBetaGamma, calParMeta
  implicit none
  ! local variables
  integer(i4b)                 :: iPar          ! loop indices
  
  nParCalSum=0_i4b
  do iPar=1,nBetaGamma        ! cound number of calibration beta (per layer) and gamma
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
  use globalData,  only:parSubset, nBetaGamma, betaCalScale 
  implicit none
  ! output variables
  real(dp),dimension(:,:),   intent(out)     :: param 
  logical,dimension(:),      intent(out)     :: mask
  ! local variables
  integer(i2b)                               :: iPar    ! loop indices
  integer(i2b)                               :: idx     ! count of calibrating parameter including per layer parameter 
  integer(i2b)                               :: ixHV    ! count of calibrating parameter including per layer parameter 
  
  idx=0_i2b
  do iPar=1,nBetaGamma
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
