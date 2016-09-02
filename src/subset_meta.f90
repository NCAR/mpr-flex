module subset_meta 

  use nrtype
  use public_var

  implicit none

  private
  public::get_parm_meta
  public::param_setup 

contains

! ************************************************************************************************
! Subroutine: read calibrating parameter metadata from a meta file 
! ************************************************************************************************
subroutine get_parm_meta(infile, err, message)

  ! used to read metadata from an input file and populate the appropriate metadata structure
  use data_type,  only:cpar_meta              ! metadata structure
  use globalData, only:parMaster, parSubset
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
  integer(i4b)                         :: unt            ! DK: need to either define units globally, or use getSpareUnit
  integer(i4b)                         :: iline          ! loop through lines in the file 
  integer(i4b)                         :: ixLocal        ! index for calibrationg parameter list 
  integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
  character(LEN=256)                   :: temp           ! single lime of information
  integer(i4b)                         :: iend           ! check for the end of the file
  character(LEN=256)                   :: ffmt           ! file format
  character(len=1)                     :: dLim(1)        ! column delimiter
  integer(i4b)                         :: ivar           ! index of master parameter  
 
  ! Start procedure here
  err=0; message="get_param_meta/"
  
  allocate(parSubset(nParCal))

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
    parSubset(ixLocal)%ptype    = parMaster(ivar)%ptype
    ixLocal = ixLocal+1
  enddo  ! looping through lines in the file
 
  ! check that all elements are populated
  if(any(parSubset(:)%pname==''))then
    do iline=1,size(parSubset)
      print*,iline,' -> ',trim(parSubset(iline)%pname)
    end do
    err=40; message=trim(message)//"someVariablesNotPopulated"; return
  endif
 
  ! close file unit
  close(unt)

end subroutine get_parm_meta

! ************************************************************************************************
! Subroutine: convert parameter data structure to simple arrays 
! ************************************************************************************************
subroutine param_setup( param,ptype, mask)
  use globalData,  only:parSubset
  implicit none
  ! output variables
  real(dp),dimension(:,:),   intent(out)     :: param 
  integer(i4b),dimension(:), intent(out)     :: ptype 
  logical,dimension(:),      intent(out)     :: mask
  ! local variables
  integer(i4b)                               :: iPar  ! loop indices
  
  do iPar=1,nParCal
    param(iPar,1) = parSubset(iPar)%val
    param(iPar,2) = parSubset(iPar)%lwr
    param(iPar,3) = parSubset(iPar)%upr
    ptype(iPar)   = parSubset(iPar)%ptype
    mask(iPar)    = parSubset(iPar)%flag
  enddo  

end subroutine param_setup

end module subset_meta 
