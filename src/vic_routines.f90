module vic_routines
! Routines specific to VIC model
  use nrtype 
  use public_var
  use strings
  use data_type 

  implicit none

  private

  public :: adj_soil_param_vic 
  public :: adj_vege_param_vic 
  public :: replace_soil_param_vic 
  public :: vic_hru_id
  public :: vic_soil_layer
  public :: read_vic_sim
  public :: read_soil_param_vic
  public :: write_soil_param_vic

contains

!***************************
! write VIC soil parameters 
!***************************
subroutine write_soil_param_vic(hruid, param, ierr, message)
  implicit none
  !input variables
  integer(i4b),intent(in)            :: hruid(:)     ! hru ID
  real(dp),    intent(in)            :: param(:,:)   ! 
  ! output
  integer(i4b),intent(out)           :: ierr         ! error code
  character(*),intent(out)           :: message      ! error message
  ! local variables
  integer(i4b)                       :: iHru         ! loop index

  ! initialize error control
  ierr=0; message='write_soil_param_vic/'
  if(size(param,2)/=TotNpar)then;ierr=10;message=trim(message)//"params 2nd dimension size different than TotNpar";return;endif
  if(size(param,1)/=nHru)then;ierr=11;message=trim(message)//'params 1st dimension size different than nHru';return;endif
  if(size(hruid)/=nHru)then;ierr=12;message=trim(message)//'hruid size different than nHru';return;endif
  open(UNIT=51,file=trim(calibparam_name),action='write',status='unknown' )
  hru:do iHru = 1,nHru
    write(51,'(I,2X)',advance='no') 1
    write(51,'(I8,2X)',advance='no') hruid(iHru)
    write(51,'(f9.4,X)',advance='no') param(iHru,3:4)
    write(51,'(f9.5,X)',advance='no') param(iHru,5)
    write(51,'(f9.4,X)',advance='no') param(iHru,6:9)
    write(51,'(f9.4,X)',advance='no') param(iHru,10:12)
    write(51,'(f10.4,X)',advance='no')param(iHru,13:15)
    write(51,'(f7.1,X)',advance='no') param(iHru,16:18)
    write(51,'(f10.4,X)',advance='no')param(iHru,19:52)
    write(51,'(I2,X)',advance='no') int(param(iHru,53))
    write(51,'(f9.4)') param(iHru,54)
  enddo hru 
  close(UNIT=51)
  return
end subroutine

!***************************
! Read VIC hru IDs 
!***************************
subroutine vic_hru_id(hruid, err, message)
  implicit none
  ! input 
  ! output
  integer(i4b),intent(out)                   :: hruid(:)     ! list of hru where calibration is performed 
  integer(i4b),intent(out)                   :: err          ! error code
  character(*),intent(out)                   :: message      ! error message
  ! local variables
  real(dp),dimension(TotNpar)                :: realline
  integer(i4b)                               :: ipar,iHru    ! loop index
  integer(i4b)                               :: stat

  ! initialize error control
  err=0; message='vic_hru_id/'
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
 ! Read original soil parameter file
  do iHru = 1,nHru
    read(unit=50,fmt=*) (realline(ipar), ipar=1,TotNpar)
    hruid(iHru)=realline(2)
  end do
  close(UNIT=50)
  return
end subroutine

!***************************
! Read VIC soil layer parameters 
!***************************
subroutine vic_soil_layer(hlyr, err, message)
  implicit none
  ! input 
  ! output
  real(dp),    intent(out)           :: hlyr(:,:) ! calibrating parameter list 
  integer(i4b),intent(out)           :: err       ! error code
  character(*),intent(out)           :: message   ! error message
  ! local variables
  real(dp),dimension(TotNpar)        :: realline
  integer(i4b)                       :: ipar,iHru  ! loop index
  integer(i4b)                       :: stat

  ! initialize error control
  err=0; message='vic_soil_layer/'
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
 ! Read original soil parameter file
  do iHru = 1,nHru
    read(unit=50,fmt=*) (realline(ipar), ipar=1,TotNpar)
    hlyr(iHru,:)=realline(4*nLyr+11:5*nLyr+10)
  end do
  close(UNIT=50)
  return
end subroutine

!***************************
! Read VIC soil parameters 
!***************************
subroutine read_soil_param_vic(param, err, message)
  implicit none
  ! input 
  ! output
  real(dp),    intent(out)                   :: param(:,:)   ! calibrating parameter list 
  integer(i4b),intent(out)                   :: err          ! error code
  character(*),intent(out)                   :: message      ! error message
  ! local variables
  integer(i4b)                               :: ipar,iHru    ! loop index
  integer(i4b)                               :: stat

  ! initialize error control
  err=0; message='vic_soil_param/'
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
 ! Read original soil parameter file
  do iHru = 1,nHru
    read(unit=50,fmt=*) (param(iHru,ipar), ipar=1,TotNpar)
  end do
  close(UNIT=50)
  return
end subroutine

!***************************
! replace VIC soil parameters 
!***************************
subroutine replace_soil_param_vic(param, hModel, parMxyMz, adjParam, ierr, message)
  use globalData, only: parMaster, betaInGamma 
  use get_ixname, only: get_ixPar
  implicit none
  !input variables
  real(dp),         intent(in)   :: param(:,:)    ! original soil parameters 
  real(dp),         intent(in)   :: hModel(:,:)   ! Model layer thickness at model hrus 
  type(namedvar2),  intent(in)   :: parMxyMz(:)   ! soil model parameter at model layer x model hrus
  ! output
  real(dp),         intent(out)  :: adjParam(:,:) ! adjusted soil parameter
  integer(i4b),     intent(out)  :: ierr          ! error code
  character(*),     intent(out)  :: message       ! error message
  ! local variables
  integer(i4b)                   :: ipar,iHru     ! loop index
  integer(i4b)                   :: nSoilParModel ! number of parameters 

  ! initialize error control
  ierr=0; message='replace_soil_param_vic/'
  nSoilParModel=size(parMxyMz)
  adjParam=param
  hru: do iHru = 1,nHru
    ! replace parameter values
    do iPar=1,nSoilParModel
      associate( ix=>get_ixPar(trim(betaInGamma(iPar))) )
      select case( parMaster(ix)%pname )
        case('binfilt');  adjParam(iHru,5)     = parMxyMz(iPar)%varData(1, iHru) 
        case('D1');       adjParam(iHru,6)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('D2');       adjParam(iHru,7)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('D3');       adjParam(iHru,8)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('D4');       adjParam(iHru,9)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('expt');     adjParam(iHru,10:12) = parMxyMz(iPar)%varData(:,iHru)
        case('ks');       adjParam(iHru,13:15) = parMxyMz(iPar)%varData(:,iHru)
        case('h1');       adjParam(iHru,23)    = hModel(1,iHru)
        case('h2');       adjParam(iHru,24)    = hModel(2,iHru)
        case('h3');       adjParam(iHru,25)    = hModel(3,iHru)
        case('bbl');      adjParam(iHru,28:30) = parMxyMz(iPar)%varData(:,iHru)
        case('BD');       adjParam(iHru,34:36) = parMxyMz(iPar)%varData(:,iHru)
        case('SD');       adjParam(iHru,37:39) = parMxyMz(iPar)%varData(:,iHru)
        case('WcrFrac');  adjParam(iHru,41:43) = parMxyMz(iPar)%varData(:,iHru)
        case('WpwpFrac'); adjParam(iHru,44:46) = parMxyMz(iPar)%varData(:,iHru)
      end select
      end associate
    end do
  enddo hru
  return
end subroutine

!***************************
! Adjust VIC soil parameters with multipliers 
!***************************
subroutine adj_soil_param_vic(param, multiplier, adjParam,  err, message)
!! This routine takes the adjustable parameter set "param" from namelist, reads into "origparam_name",
!! computes the new parameters, writes them into "calibparam_name" 
  use globalData, only: parSubset
  implicit none
  !input variables
  real(dp),    intent(in)    :: param(:,:)    ! original soil parameters 
  real(dp),    intent(in)    :: multiplier(:) ! mulitpliers for calibrating soil parameter 
  ! output
  real(dp),    intent(out)   :: adjParam(:,:) ! adjusted soil parameter
  integer(i4b),intent(out)   :: err           ! error code
  character(*),intent(out)   :: message       ! error message
  ! local variables
  integer(i4b)               :: iPar,iHru     ! loop index

  ! initialize error control
  err=0; message='adj_soil_param_vic/'
  adjParam=param
  do iHru = 1,nHru
    ! Modify parameter values
    do iPar=1,nParCal
      select case( parSubset(iPar)%pname )
        case('binfilt');  adjParam(iHru,5)     = multiplier( iPar )*Param(iHru,5)
        case('D1');       adjParam(iHru,6)     = multiplier( iPar )*Param(iHru,6)
        case('D2');       adjParam(iHru,7)     = multiplier( iPar )*Param(iHru,7)
        case('D3');       adjParam(iHru,8)     = multiplier( iPar )*Param(iHru,8)
        case('D4');       adjParam(iHru,9)     = multiplier( iPar )*Param(iHru,9)
        case('expt');     adjParam(iHru,10:12) = multiplier( iPar )*Param(iHru,10:12)
        case('ks');       adjParam(iHru,13:15) = multiplier( iPar )*Param(iHru,13:15)
        case('h1');       adjParam(iHru,23)    = multiplier( iPar )*Param(iHru,23)
        case('h2');       adjParam(iHru,24)    = multiplier( iPar )*Param(iHru,24)
        case('h3');       adjParam(iHru,25)    = multiplier( iPar )*Param(iHru,25)
        case('bbl');      adjParam(iHru,28:30) = multiplier( iPar )*Param(iHru,28:30)
        case('BD');       adjParam(iHru,34:36) = multiplier( iPar )*Param(iHru,34:36)
        case('SD');       adjParam(iHru,37:39) = multiplier( iPar )*Param(iHru,37:39)
        case('WcrFrac');  adjParam(iHru,41:43) = multiplier( iPar )*Param(iHru,41:43)
        case('WpwpFrac'); adjParam(iHru,44:46) = multiplier( iPar )*Param(iHru,44:46)
       end select
    end do
    ! Limit parameters to correct possible values without physical meaning: this applies for all configurations
    !binfilt
    if(adjParam(iHru,5) .lt. 0.001) then
      adjParam(iHru,5) = 0.001
    elseif(adjParam(iHru,5) .gt. 0.5) then
      adjParam(iHru,5) = 0.5
    endif
    !Ds
    if(adjParam(iHru,6) .lt. 0.0001) then
      adjParam(iHru,6) = 0.0001
    elseif(adjParam(iHru,6) .gt. 1.0) then
      adjParam(iHru,6) = 1.0
    endif
    !Dsmax
    if(adjParam(iHru,7) .lt. 0.0001) then
      adjParam(iHru,7) = 0.00001
    elseif(adjParam(iHru,7) .gt. 1.0) then
      adjParam(iHru,7) = 1.0
    endif
    !Ws
    if(adjParam(iHru,8) .lt. 0.0001) then
      adjParam(iHru,8) = 0.0001
    elseif(adjParam(iHru,8) .gt. 1000) then
      adjParam(iHru,8) = 1000.0 
    endif
    !bulk density for each layer
    do iPar = 34,36
      if(adjParam(iHru,iPar) .lt. 805.) then
        adjParam(iHru,iPar) = 805.
      elseif(adjParam(iHru,iPar) .gt. 1880.) then
        adjParam(iHru,iPar) = 1880.
      endif
    enddo
  enddo 
  return
end subroutine

!***************************
! Adjust VIC vege parameters with multiplier
!***************************
subroutine adj_vege_param_vic(multiplier, err, message)
  use globalData, only: parSubset
  implicit none

  ! input variables
  real(dp),             intent(in) :: multiplier(:)             ! list of calibratin parameters 
  ! output
  integer(i4b),intent(out)         :: err                       ! error code
  character(*),intent(out)         :: message                   ! error message
  ! local variables
  integer(i4b)                     :: vegClass                  ! vegetation class 
  real(dp)                         :: vegFrac                   ! fraction of vage class
  real(dp),dimension(nLyr)         :: rootDepth                 ! root zone depth
  real(dp),dimension(nLyr)         :: rootFrac                  ! root zone fraction
  real(dp),dimension(12)           :: laiMonth                  ! monthly LAI
  integer(i4b)                     :: hruID                     ! hru ID
  integer(i4b)                     :: nTile                     ! number of vege tile 
  integer(i4b)                     :: iPar,iHru,iTile,iMon,iLyr ! loop index
  character(50)                    :: rowfmt                    ! string specifying write format for real value
  integer(i4b)                     :: stat

  ! initialize error control
  err=0; message='adj_vege_param_vic/'
  !Open original and modified vege parameter files
  open (UNIT=50,file=origvege_name,form='formatted',status='old',IOSTAT=stat)
  open (UNIT=51,file=calivege_name,action='write',status='replace' )
  write(rowfmt,'(A,I2,A)') '(',nLyr,'(1X,F4.2))'
 ! Read original vege parameter file
  hru:do iHru = 1,nHru
    read(unit=50,fmt=*) hruID,nTile
    write(51,'(I10,1X,I2)') (hruID,nTile)
    tile:do iTile = 1,nTile
      read(unit=50,fmt=*) vegClass,vegFrac,(rootDepth(iLyr), iLyr=1,nLyr),(rootFrac(iLyr), iLyr=1,nLyr)
      read(unit=50,fmt=*) (laiMonth(iMon), iMon=1,12)
      ! Modify parameter values
      par:do iPar=1,nParCal
        select case( parSubset(iPar)%pname )
          case('lai');    laiMonth = multiplier( iPar )*laiMonth
        end select
      enddo par
      ! Write the modified parameter file for the entire basin/region for traditional upscaling
      write(51,'(3X,I2,1X,F8.6)',advance='no') (vegClass,vegFrac)
      write(51,rowfmt,advance='no')            (rootDepth(iLyr), iLyr=1,nLyr)
      write(51,rowfmt)                         (rootFrac(iLyr), iLyr=1,nLyr)
      write(51,'(5X,12(1X,F6.3))')             (laiMonth(iMon), iMon=1,12)
    enddo tile 
  enddo hru
  ! Close original and modified basin parameter files
  close(UNIT=50)
  close(UNIT=51)
  return
end subroutine

!***************************
! Read VIC output file
!***************************
subroutine read_vic_sim(sim, err, message)
  implicit none
  !output variables
  real(dp),              intent(out) :: sim(:,:)
  integer(i4b),          intent(out) :: err            ! error code
  character(*),          intent(out) :: message        ! error message
  !local variables
  character(len=strLen)              :: filename
  real(dp)                           :: cellfraction,basin_area
  real(dp)                           :: auxflux(5)                 ! This is only in case of water balance mode
  integer(i4b)                       :: ibasin, itime, ivar, icell ! index 
  integer(i4b)                       :: ncell
  integer(i4b)                       :: dum,c_cell

  ! initialize error control
  err=0; message='read_vic_sim/'
  !set output variable to zero
  sim = 0.0_dp
  !cell counter
  c_cell = 1
  !open a few files
  open (UNIT=53,file=trim(filelist_name),form='formatted',status='old')
  open (UNIT=54,file=trim(cellfrac_name),form='formatted',status='old')
  open (UNIT=51,file=trim(region_info),form='formatted',status='old')
  do ibasin = 1,nbasin
    read (UNIT=51,fmt=*) dum,dum,basin_area,ncell
    do icell = 1,ncell
      read (UNIT=53,fmt=*) filename
      read (UNIT=54,fmt=*) cellfraction
      filename=trim(sim_dir)//trim(filename)
      open (UNIT=55,file= filename,form='formatted',status='old')
      do itime = 1,sim_len
        read (UNIT=55,fmt=*) (auxflux(ivar), ivar=1,5)
        sim(c_cell,itime) = (auxflux(4) + auxflux(5))*cellfraction
      enddo
      close(UNIT=55)
      c_cell = c_cell + 1
    enddo
  enddo
  close(UNIT=51)
  close(UNIT=53)
  close(UNIT=54)
  return
end subroutine

end module vic_routines
