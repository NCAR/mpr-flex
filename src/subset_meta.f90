module subset_meta 

  use nrtype
  use public_var

  implicit none

  private
  public::get_parm_meta
  public::param_setup 
  public::check_gammaZ
  public::check_gammaH

contains

! ************************************************************************************************
! Subroutine: read calibrating parameter metadata from a meta file 
! ************************************************************************************************
subroutine get_parm_meta(infile, err, message)
  ! used to read parameter list ascii and obtain meta data from master parameter list (paramMaster.f90) 
  use data_type,  only:cpar_meta
  use globalData, only:parMaster, parSubset, gammaSubset, betaInGamma  ! data type definition
  use ascii_util, only:file_open
  use get_ixname, only:get_ixPar
 
  implicit none
 
  ! input
  character(*),intent(in)              :: infile         ! input filename
  ! output
  integer(i4b),intent(out)             :: err            ! error code
  character(*),intent(out)             :: message        ! error message
  ! local variables
  character(len=256)                   :: cmessage       ! error message for downwind routine
  type(cpar_meta),allocatable          :: tempMeta(:)
  character(len=strLen),allocatable    :: res(:)         ! 
  integer(i4b)                         :: unt            ! DK: need to either define units globally, or use getSpareUnit
  integer(i4b)                         :: iline          ! loop through lines in the file 
  integer(i4b)                         :: ixLocal        ! index for calibrationg parameter list 
  integer(i4b)                         :: ixGamma        ! index for calibrationg gamma parameter list 
  integer(i4b)                         :: ivar           ! index of master parameter  
  integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
  character(LEN=256)                   :: temp           ! single lime of information
  integer(i4b)                         :: iend           ! check for the end of the file
  character(LEN=256)                   :: ffmt           ! file format
  character(len=1)                     :: dLim(1)        ! column delimiter
  integer(i4b)                         :: i,j,k          ! loop index
 
  ! Start procedure here
  err=0; message="get_param_meta/"
  allocate(parSubset(nParCal))
  allocate(tempMeta(nParCal))
  ! open file
  call file_open(trim(infile),unt,err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! get to the start of the variable descriptions 
  do iline=1,maxLines
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
    if (temp(1:1)/='!') exit                             ! assume first line not comment is format code
  end do ! looping through file to find the format code
  ! read in format string
  read(temp,*)ffmt
  ! loop through the lines in the file
  ixLocal=1
  ixGamma=1
  do iline=1,maxLines
    ! read a line of data and exit iif an error code (character read, so only possible error is end of file)
    read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
    ! check that the line is not a comment
    if (temp(1:1)=='!')cycle
    ! save data into a temporary structure
    read(temp,trim(ffmt),iostat=err) parSubset(ixLocal)%pname ,dLim(1),&   ! parameter name
                                     parSubset(ixLocal)%flag               ! flag
    if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
    ! check that the delimiters are in the correct place
    if(any(dLim /= '|'))then
      message=trim(message)//'delimiter is not in the correct place; line = ['//trim(temp)//']; filename = '//trim(infile)
      err=32; return
    endif
    ! identify the index of the named variable from master parameter 
    ivar = get_ixPar(parSubset(ixLocal)%pname)
    if(ivar<=0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(parSubset(ixLocal)%pname)//"]"; return; endif
    ! copy the index of the named variable from master parameter and their info (init.val, lwr, upr)
    parSubset(ixLocal)%ixMaster = ivar
    parSubset(ixLocal)%val      = parMaster(ivar)%val
    parSubset(ixLocal)%lwr      = parMaster(ivar)%lwr
    parSubset(ixLocal)%upr      = parMaster(ivar)%upr
    parSubset(ixLocal)%beta     = parMaster(ivar)%beta
    parSubset(ixLocal)%ptype    = parMaster(ivar)%ptype
    parSubset(ixLocal)%hups     = parMaster(ivar)%hups
    parSubset(ixLocal)%vups     = parMaster(ivar)%vups
    ixLocal = ixLocal+1
    ! extract only gamma parameter list (consider putting this structure in private here)
    associate( parentp => parSubset(ixLocal-1)%beta )
    if ( trim(parentp)/="beta") then 
      tempMeta(ixGamma)%ixMaster = ivar
      tempMeta(ixGamma)%pname    = parMaster(ivar)%pname
      tempMeta(ixGamma)%val      = parMaster(ivar)%val
      tempMeta(ixGamma)%lwr      = parMaster(ivar)%lwr
      tempMeta(ixGamma)%upr      = parMaster(ivar)%upr
      tempMeta(ixGamma)%beta     = trim(parentp) 
      tempMeta(ixGamma)%ptype    = parMaster(ivar)%ptype
      tempMeta(ixGamma)%hups     = parMaster(ivar)%hups
      tempMeta(ixGamma)%vups     = parMaster(ivar)%vups
      ixGamma = ixGamma+1
    endif
    end associate
  enddo  ! looping through lines in the file
  if (ixGamma > 1) then
    allocate(gammaSubset(ixGamma-1))
    gammaSubset=tempMeta(1:ixGamma) 
  endif
  ! get beta parameter associated with gamma parameter
  if ( allocated(gammaSubset) ) then
    allocate(res(size(gammaSubset)))
    k = 1
    res(1) = gammaSubset(1)%beta
    outer: do i=2,size(gammaSubset)
      do j=1,k
        ! if find a match so start looking again
        if (res(j)==gammaSubset(i)%beta)then; cycle outer; endif
      end do
      ! No match found so add it to the output
      k = k + 1
      res(k) = gammaSubset(i)%beta
    end do outer
    allocate(betaInGamma(k))
    betaInGamma=res(1:k) 
  else
    print*, 'NO gamma parameters to be calibrated'
  endif
  ! check that all elements are populated
  if(any(parSubset(:)%pname==''))then
    do iline=1,size(parSubset)
      print*,iline,' -> ',trim(parSubset(iline)%pname)
    end do
    err=40; message=trim(message)//"someVariablesNotPopulated"; return
  endif
  ! close file unit
  close(unt)
  return
end subroutine get_parm_meta

!**********************************
! Subroutine: check if z parameters exist in gamma parameter
!**********************************
subroutine check_gammaZ( err, message)
  use globalData,   only: gammaSubset
  implicit none
  !output variables
  integer(i4b),                      intent(out)   :: err         ! error code
  character(*),                      intent(out)   :: message     ! error message
  !local variables
  logical(lgc)                                     :: checkZ 
  character(len=strLen)                            :: cmessage    ! error message from downward subroutine
  integer(i4b)                                     :: unt         ! DK: need to either define units globally, or use getSpareUnit
  integer(i4b)                                     :: i           ! loop index

  ! initialize error control
  err=0; message='check_gammaZ/'
  checkZ=.false.
  if ( allocated(gammaSubset) )then
    !check z parameter
    do i=1,size(gammaSubset)
      if (gammaSubset(i)%pname=="z1gamma1") checkZ=.true.
    end do
    if ( .not. checkZ ) stop 'Calibrating Gamma parameter require z1gamma1'
  else
    print*, 'NO gamma parameters'
  endif
  return
end subroutine check_gammaZ

!**********************************
! Subroutine: check if h parameters exist in gamma parameter
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
  character(len=strLen)             :: cmessage    ! error message from downward subroutine
  integer(i4b)                      :: unt         ! DK: need to either define units globally, or use getSpareUnit
  integer(i4b)                      :: i,j 

  ! initialize error control
  err=0; message='check_gammaH/'
  allocate(checkH(nLyr-1))
  id=-999
  if ( allocated(gammaSubset) )then
    !check h parameters - now can chcek up to 5 layers
    do i=1,size(gammaSubset)
      if (gammaSubset(i)%pname=="h1gamma1")then;id(1)=1;cycle;endif 
      if (gammaSubset(i)%pname=="h1gamma2")then;id(2)=1;cycle;endif
      if (gammaSubset(i)%pname=="h1gamma3")then;id(3)=1;cycle;endif
      if (gammaSubset(i)%pname=="h1gamma4")then;id(4)=1;cycle;endif
    enddo
    mask=(id>0)
    checkH=mask(1:nLyr-1) 
    if ( any(.not. checkH) ) stop 'Calibrating Gamma parameter require (nLyr-1) hgamma parameters '
  else
    print*, 'NO gamma parameters'
  endif
  return
end subroutine check_gammaH

! ************************************************************************************************
! Subroutine: convert parameter data structure to simple arrays 
! ************************************************************************************************
subroutine param_setup( param, mask)
  use globalData,  only:parSubset
  implicit none
  ! output variables
  real(dp),dimension(:,:),   intent(out)     :: param 
  logical,dimension(:),      intent(out)     :: mask
  ! local variables
  integer(i4b)                               :: iPar  ! loop indices
  
  do iPar=1,nParCal
    param(iPar,1) = parSubset(iPar)%val
    param(iPar,2) = parSubset(iPar)%lwr
    param(iPar,3) = parSubset(iPar)%upr
    mask(iPar)    = parSubset(iPar)%flag
  enddo  

end subroutine param_setup

end module subset_meta 
