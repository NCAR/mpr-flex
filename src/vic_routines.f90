module vic_routines
! Routines specific to VIC model
  use nrtype 
  use public_var
  use strings
  use data_type 

  implicit none
  public :: adj_soil_param_vic 
  public :: adj_vege_param_vic 
  public :: replace_soil_param_vic 
  public :: vic_hru_id
  public :: vic_soil_layer
  public :: read_vic_sim
  public :: read_soil_param_vic
  private

contains

!***************************
! write VIC soil parameters 
!***************************
subroutine write_soil_param_vic(fname, hruid, param, ierr, message)
  implicit none
  !input variables
  character(*),intent(in)            :: fname        ! file name to be written 
  integer(i4b),intent(in)            :: hruid(:)     ! hru ID
  real(dp),    intent(in)            :: param(:,:)   ! 
  ! output
  integer(i4b),intent(out)           :: ierr         ! error code
  character(*),intent(out)           :: message      ! error message
  ! local variables
  integer(i4b)                       :: iHru         ! loop index
  integer(i4b)                       :: nHru         ! number of hrus 
  integer(i4b)                       :: nPar         ! number of parameters 
  integer(i4b)                       :: stat

  ! initialize error control
  ierr=0; message='write_soil_param_vic/'
  nHru=size(param,1)
  nPar=size(param,2)
  if ( nHru/=size(hruid) )then;ierr=10;message=trim(message)//'hruid size and param size not matched';return;endif
  open(UNIT=51,file=trim(fname),action='write',status='unknown' )
  do iHru = 1,nHru
    ! Write the modified parameter file for the entire basin/region for traditional upscaling
      write(51,'(I,2X)',advance='no') 1
      write(51,'(I8,2X)',advance='no') hruid(2)
      write(51,'(f9.4,X)',advance='no') param(iHru,3:4)
      write(51,'(f9.5,X)',advance='no') param(iHru,5)
      write(51,'(f9.4,X)',advance='no') param(iHru,6:8)
      write(51,'(I3,2X)',advance='no') int(param(iHru,9))
      write(51,'(f9.4,X)',advance='no') param(iHru,10:12)
      write(51,'(f10.4,X)',advance='no')param(iHru,13:15)
      write(51,'(f7.1,X)',advance='no') param(iHru,16:18)
      write(51,'(f10.4,X)',advance='no')param(iHru,19:52)
      write(51,'(I2,X)',advance='no') int(param(iHru,53))
      write(51,'(f9.4)') param(iHru,54)
  enddo  !end cell loop
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
    read(unit=50,*) (realline(ipar), ipar=1,TotNpar)
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
  real(dp),    intent(out)                   :: hlyr(:,:) ! calibrating parameter list 
  integer(i4b),intent(out)                   :: err             ! error code
  character(*),intent(out)                   :: message         ! error message
  ! local variables
  real(dp),dimension(TotNpar)             :: realline
  integer(i4b)                               :: ipar,iHru,iLyr  ! loop index
  integer(i4b)                               :: stat

  ! initialize error control
  err=0; message='vic_soil_layer/'
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
 ! Read original soil parameter file
  do iHru = 1,nHru
    read(unit=50,*) (realline(ipar), ipar=1,TotNpar)
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
    read(unit=50,*) (param(iHru,ipar), ipar=1,TotNpar)
  end do
  close(UNIT=50)
  return
end subroutine

!***************************
! replace VIC soil parameters 
!***************************
subroutine replace_soil_param_vic(hModel, parMxyMz, ierr, message)
  use globalData, only: parMaster, betaInGamma 
  use get_ixname, only: get_ixPar
  implicit none
  !input variables
  real(dp),         intent(in)   :: hModel(:,:)  ! Model layer thickness for one hru
  type(namedvar2),  intent(in)   :: parMxyMz(:)  ! 
  ! output
  integer(i4b),     intent(out)  :: ierr         ! error code
  character(*),     intent(out)  :: message      ! error message
  ! local variables
  integer(i4b)                   :: ipar,iHru     ! loop index
  integer(i4b)                   :: nSoilParModel ! number of parameters 
  integer(i4b)                   :: stat
  real(dp),dimension(TotNpar)    :: realline

  ! initialize error control
  ierr=0; message='replace_soil_param_vic/'
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
  open (UNIT=51,file=calibparam_name,action='write',status='unknown' )
  nSoilParModel=size(parMxyMz)
 ! Read original soil parameter file
  do iHru = 1,nHru
    read(unit=50,*) (realline(ipar), ipar=1,TotNpar)
    ! replace parameter values
    do iPar=1,nSoilParModel
      associate( ix=>get_ixPar(trim(betaInGamma(iPar))) )
      select case( parMaster(ix)%pname )
        case('binfilt');  realline(5)     = parMxyMz(iPar)%varData(1, iHru) 
        case('D1');       realline(6)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('D2');       realline(7)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('D3');       realline(8)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('D4');       realline(9)     = parMxyMz(iPar)%varData(nLyr,iHru)
        case('expt');     realline(10:12) = parMxyMz(iPar)%varData(:,iHru)
        case('ks');       realline(13:15) = parMxyMz(iPar)%varData(:,iHru)
        case('h1');       realline(23)    = hModel(1,iHru)
        case('h2');       realline(24)    = hModel(2,iHru)
        case('h3');       realline(25)    = hModel(3,iHru)
        case('bbl');      realline(28:30) = parMxyMz(iPar)%varData(:,iHru)
        case('BD');       realline(34:36) = parMxyMz(iPar)%varData(:,iHru)
        case('SD');       realline(37:39) = parMxyMz(iPar)%varData(:,iHru)
        case('WcrFrac');  realline(41:43) = parMxyMz(iPar)%varData(:,iHru)
        case('WpwpFrac'); realline(44:46) = parMxyMz(iPar)%varData(:,iHru)
      end select
      end associate
    end do
    ! Write the modified parameter file for the entire basin/region for traditional upscaling
    write(51,'(I,2X)',advance='no') 1
    write(51,'(I8,2X)',advance='no') int(realline(2))
    write(51,'(f9.4,X)',advance='no') realline(3:4)
    write(51,'(f9.5,X)',advance='no') realline(5)
    write(51,'(f9.4,X)',advance='no') realline(6:8)
    write(51,'(I3,2X)',advance='no') int(realline(9))
    write(51,'(f9.4,X)',advance='no') realline(10:12)
    write(51,'(f10.4,X)',advance='no') realline(13:15)
    write(51,'(f7.1,X)',advance='no') realline(16:18)
    write(51,'(f10.4,X)',advance='no') realline(19:52)
    write(51,'(I2,X)',advance='no') int(realline(53))
    write(51,'(f9.4)') realline(54)
  enddo  !end cell loop
  ! Close original and modified basin parameter files
  close(UNIT=50)
  close(UNIT=51)
  return
end subroutine

!***************************
! Adjust VIC soil parameters with multipliers 
!***************************
subroutine adj_soil_param_vic(param, err, message)
!! This routine takes the adjustable parameter set "param" from namelist, reads into "origparam_name",
!! computes the new parameters, writes them into "calibparam_name" 
  use globalData, only: parSubset
  implicit none

  !input variables
  real(dp),dimension(:),intent(in)   :: param        ! calibrating parameter list 
  ! output
  integer(i4b),intent(out)           :: err          ! error code
  character(*),intent(out)           :: message      ! error message
  ! local variables
  integer(i4b)                       :: ipar,iHru    ! loop index
  integer(i4b)                       :: stat
  real(dp),dimension(TotNpar)        :: realline

  ! initialize error control
  err=0; message='adj_soil_param_vic/'
 !Open original and modified basin parameter files
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
  open (UNIT=51,file=calibparam_name,action='write',status='unknown' )
 ! Read original soil parameter file
  do iHru = 1,nHru
    read(unit=50,*) (realline(ipar), ipar=1,TotNpar)
    ! Modify parameter values
    do iPar=1,nParCal
      select case( parSubset(iPar)%pname )
        case('binfilt');  realline(5)     = param( iPar )*realline(5)
        case('D1');       realline(6)     = param( iPar )*realline(6)
        case('D2');       realline(7)     = param( iPar )*realline(7)
        case('D3');       realline(8)     = param( iPar )*realline(8)
        case('D4');       realline(9)     = param( iPar )*realline(9)
        case('expt');     realline(10:12) = param( iPar )*realline(10:12)
        case('ks');       realline(13:15) = param( iPar )*realline(13:15)
        case('h1');       realline(23)    = param( iPar )*realline(23)
        case('h2');       realline(24)    = param( iPar )*realline(24)
        case('h3');       realline(25)    = param( iPar )*realline(25)
        case('bbl');      realline(28:30) = param( iPar )*realline(28:30)
        case('BD');       realline(34:36) = param( iPar )*realline(34:36)
        case('SD');       realline(37:39) = param( iPar )*realline(37:39)
        case('WcrFrac');  realline(41:43) = param( iPar )*realline(41:43)
        case('WpwpFrac'); realline(44:46) = param( iPar )*realline(44:46)
       end select
    end do
    ! Limit parameters to correct possible values without physical meaning: this applies for all configurations
    !binfilt
    if(realline(5) .lt. 0.001) then
      realline(5) = 0.001
    elseif(realline(5) .gt. 0.5) then
      realline(5) = 0.5
    endif
    !Ds
    if(realline(6) .lt. 0.0001) then
      realline(6) = 0.0001
    elseif(realline(6) .gt. 1.0) then
      realline(6) = 1.0
    endif
    !Dsmax
    if(realline(7) .lt. 0.0001) then
      realline(7) = 0.00001
    elseif(realline(7) .gt. 1.0) then
      realline(7) = 1.0
    endif
    !Ws
    if(realline(8) .lt. 0.0001) then
      realline(8) = 0.0001
    elseif(realline(8) .gt. 1000) then
      realline(8) = 1000.0 
    endif
    !bulk density for each layer
    do iPar = 34,36
      if(realline(iPar) .lt. 805.) then
        realline(iPar) = 805.
      elseif(realline(iPar) .gt. 1880.) then
        realline(iPar) = 1880.
      endif
    enddo
    ! Write the modified parameter file for the entire basin/region for traditional upscaling
      write(51,'(I,2X)',advance='no') 1
      write(51,'(I8,2X)',advance='no') int(realline(2))
      write(51,'(f9.4,X)',advance='no') realline(3:4)
      write(51,'(f9.5,X)',advance='no') realline(5)
      write(51,'(f9.4,X)',advance='no') realline(6:8)
      write(51,'(I3,2X)',advance='no') int(realline(9))
      write(51,'(f9.4,X)',advance='no') realline(10:12)
      write(51,'(f10.4,X)',advance='no') realline(13:15)
      write(51,'(f7.1,X)',advance='no') realline(16:18)
      write(51,'(f10.4,X)',advance='no') realline(19:52)
      write(51,'(I2,X)',advance='no') int(realline(53))
      write(51,'(f9.4)') realline(54)
  enddo  !end cell loop
  ! Close original and modified basin parameter files
  close(UNIT=50)
  close(UNIT=51)
  return
end subroutine

!***************************
! Adjust VIC vege parameters with multiplier
!***************************
subroutine adj_vege_param_vic(param, err, message)
  use globalData, only: parSubset
  implicit none

  ! input variables
  real(dp),dimension(:),intent(in) :: param                     ! list of calibratin parameters 
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
    read(unit=50,*) hruID,nTile
    write(51,'(I10,1X,I2)') (hruID,nTile)
    tile:do iTile = 1,nTile
      read(unit=50,*) vegClass,vegFrac,(rootDepth(iLyr), iLyr=1,nLyr),(rootFrac(iLyr), iLyr=1,nLyr)
      read(unit=50,*) (laiMonth(iMon), iMon=1,12)
      ! Modify parameter values
      par:do iPar=1,nParCal
        select case( parSubset(iPar)%pname )
          case('lai');    laiMonth = param( iPar )*laiMonth
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
    read (UNIT=51,*) dum,dum,basin_area,ncell
    do icell = 1,ncell
      read (UNIT=53,*) filename
      read (UNIT=54,*) cellfraction
      filename=trim(sim_dir)//trim(filename)
      open (UNIT=55,file= filename,form='formatted',status='old')
      do itime = 1,sim_len
        read (UNIT=55,*) (auxflux(ivar), ivar=1,5)
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
