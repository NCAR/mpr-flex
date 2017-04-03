module sac_routines
! Routines specific to sac model
  use nrtype 
  use public_var
  use data_type 

  implicit none

  private

  public :: adj_soil_param_sac 
  public :: adj_snow_param_sac 
  public :: replace_soil_param_sac 
  public :: sac_hru_id
  public :: sac_soil_layer
  public :: read_sac_sim
  public :: read_soil_param_sac
  public :: write_soil_param_sac

contains

!***************************
! write sac soil parameters 
!***************************
subroutine write_soil_param_sac(hruid, param, err, message)
  implicit none
  !input variables
  integer(i4b),intent(in)            :: hruid(:)     ! hru ID
  real(dp),    intent(in)            :: param(:,:)   ! parameter matrix (nHru x nParamInModel) 
  ! output
  integer(i4b),intent(out)           :: err         ! error code
  character(*),intent(out)           :: message      ! error message
  ! local variables
  character(len=strLen)              :: rowfmt
  integer(i4b)                       :: iHru         ! loop index

  ! initialize error control
  err=0; message='write_soil_param_sac/'
  if(size(param,2)/=TotNpar)then;err=10;message=trim(message)//"params 2nd dimension size different than TotNpar";return;endif
  if(size(param,1)/=nHru)then;err=11;message=trim(message)//'params 1st dimension size different than nHru';return;endif
  if(size(hruid)/=nHru)then;err=12;message=trim(message)//'hruid size different than nHru';return;endif
  open(UNIT=51,file=trim(calibparam_name),action='write',status='unknown' )
  write(rowfmt,'(A,I3,A)') '(A,',nHru,'(1X,F10.5))' 
  write(51,fmt=trim(rowfmt)) ('uztwm', param(iHru,1), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('uzfwm', param(iHru,2), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('lztwm', param(iHru,3), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('lzfpm', param(iHru,4), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('lzfsm', param(iHru,5), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adimp', param(iHru,6), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('uzk',   param(iHru,7), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('lzpk',  param(iHru,8), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('lzsk',  param(iHru,9), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('zperc', param(iHru,10),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('rexp',  param(iHru,11),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('pctim', param(iHru,12),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('pfree', param(iHru,13),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('riva',  param(iHru,14),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('side',  param(iHru,15),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('rserv', param(iHru,16),iHru=1,nHru )
  close(UNIT=51)
  return
end subroutine

!***************************
! Read sac hru IDs 
!***************************
subroutine sac_hru_id(hruid, err, message)
  implicit none
  ! input 
  ! output
  integer(i4b),intent(out)                   :: hruid(:)     ! list of hru where calibration is performed 
  integer(i4b),intent(out)                   :: err          ! error code
  character(*),intent(out)                   :: message      ! error message
  ! local variables
  character(len=strLen)                      :: dummy 
  integer(i4b)                               :: iHru         ! loop index
  integer(i4b)                               :: stat

  ! initialize error control
  err=0; message='sac_hru_id/'
  open (UNIT=50,file=trim(region_info),form='formatted',status='old',IOSTAT=stat)
  do iHru = 1,nHru
    read(unit=50,fmt=*) dummy, hruid(iHru), dummy, dummy 
  end do
  return
end subroutine

!***************************
! Read sac soil layer parameters 
!***************************
subroutine sac_soil_layer(hlyr, err, message)
  implicit none
  ! input 
  ! output
  real(dp),    intent(out)           :: hlyr(:,:)       ! soil layer thickness (bucket size) matrix (nHru x nLyr)
  integer(i4b),intent(out)           :: err            ! error code
  character(*),intent(out)           :: message         ! error message
  ! local variables
  real(dp)                           :: paramTemp(nHru) ! temporal parameter vector (nHru) 
  real(dp)                           :: hlyrTemp(nHru,5)! calibrating parameter list 
  character(len=strLen)              :: parName         ! parameter name
  integer(i4b)                       :: iPar,iHru       ! loop index
  integer(i4b)                       :: stat

  ! initialize error control
  err=0; message='sac_soil_layer/'
  if(size(hlyr,2)/=nLyr)then;err=10;message=trim(message)//"hlyr 2nd dimension size different than nLyr";return;endif
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
  ! Read original soil parameter file
  do iPar=1,TotNpar
    read(unit=50,fmt=*) parName, (paramTemp(iHru), iHru=1,nHru)
    select case(trim(parName))
      case('uztwm'); hlyrTemp(:,1)=paramTemp
      case('uzfwm'); hlyrTemp(:,2)=paramTemp
      case('lztwm'); hlyrTemp(:,3)=paramTemp
      case('lzfpm'); hlyrTemp(:,4)=paramTemp
      case('lzfsm'); hlyrTemp(:,5)=paramTemp
    end select
  end do
  hlyr(:,1)=sum(hlyrTemp(:,1:2),dim=2)
  hlyr(:,2)=sum(hlyrTemp(:,3:5),dim=2)
  return
end subroutine

!***************************
! Read sac soil parameters 
!***************************
subroutine read_soil_param_sac(param, err, message)
  implicit none
  ! input 
  ! output
  real(dp),    intent(out)                   :: param(:,:)   ! parameter matrix (nHru x nParamInModel) 
  integer(i4b),intent(out)                   :: err          ! error code
  character(*),intent(out)                   :: message      ! error message
  ! local variables
  real(dp)                                   :: paramTemp(nHru)  ! temporal parameter vector (nHru) 
  character(len=strLen)                      :: parName      ! parameter name
  integer(i4b)                               :: ipar,iHru    ! loop index
  integer(i4b)                               :: stat

  ! initialize error control
  err=0; message='sac_soil_param/'
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
  ! Read original soil parameter file
  do iPar=1,size(param,2)
    read(unit=50,fmt=*) parName, (paramTemp(iHru), iHru=1,nHru)
    select case(trim(parName))
      case('uztwm'); param(:,1)=paramTemp
      case('uzfwm'); param(:,2)=paramTemp
      case('lztwm'); param(:,3)=paramTemp
      case('lzfpm'); param(:,4)=paramTemp
      case('lzfsm'); param(:,5)=paramTemp
      case('adimp'); param(:,6)=paramTemp
      case('uzk');   param(:,7)=paramTemp
      case('lzpk');  param(:,8)=paramTemp
      case('lzsk');  param(:,9)=paramTemp
      case('zperc'); param(:,10)=paramTemp
      case('rexp');  param(:,11)=paramTemp
      case('pctim'); param(:,12)=paramTemp
      case('pfree'); param(:,13)=paramTemp
      case('riva');  param(:,14)=paramTemp
      case('side');  param(:,15)=paramTemp
      case('rserv'); param(:,16)=paramTemp
    end select
  end do
  close(UNIT=50)
  return
end subroutine

!***************************
! replace sac soil parameters 
!***************************
subroutine replace_soil_param_sac(param, parMxyMz, adjParam, err, message)
  use globalData, only: betaMaster, betaInGamma, nSoilParModel 
  use get_ixname, only: get_ixBeta
  implicit none
  !input variables
  real(dp),         intent(in)   :: param(:,:)    ! original soil parameters matrix (nHru x nParamInModel) 
  type(namedvar2),  intent(in)   :: parMxyMz(:)   ! soil model parameter at model layer x model hrus
  ! output
  real(dp),         intent(out)  :: adjParam(:,:) ! adjusted soil parameter
  integer(i4b),     intent(out)  :: err          ! error code
  character(*),     intent(out)  :: message       ! error message
  ! local variables
  integer(i4b)                   :: ipar,iHru     ! loop index

  ! initialize error control
  err=0; message='replace_soil_param_sac/'
  adjParam=param
  hru: do iHru = 1,nHru
    ! replace parameter values
    do iPar=1,nSoilParModel
      associate( ix=>get_ixBeta(trim(betaInGamma(iPar))) )
      select case( betaMaster(ix)%pname )
        case('twm')
          adjParam(iHru,1)  = parMxyMz(iPar)%varData(1,iHru)
          adjParam(iHru,3)  = parMxyMz(iPar)%varData(nLyr,iHru)
        case('fwm')
          adjParam(iHru,2)  = parMxyMz(iPar)%varData(1,iHru)
        case('fpm')
          adjParam(iHru,4)  = parMxyMz(iPar)%varData(nLyr,iHru)
        case('fsm')
          adjParam(iHru,5)  = parMxyMz(iPar)%varData(nLyr,iHru)
        case('zk')
          adjParam(iHru,7)  = parMxyMz(iPar)%varData(1,iHru)
        case('zpk')
          adjParam(iHru,8) = parMxyMz(iPar)%varData(nLyr,iHru)
        case('zsk')
          adjParam(iHru,9) = parMxyMz(iPar)%varData(nLyr,iHru)
        case('zperc')
          adjParam(iHru,10) = parMxyMz(iPar)%varData(1,iHru)
        case('rexp')
          adjParam(iHru,11) = parMxyMz(iPar)%varData(1,iHru)
        case('pfree')
          adjParam(iHru,13) = parMxyMz(iPar)%varData(1,iHru)
      end select
      end associate
    end do
  enddo hru
  return
end subroutine

subroutine adj_snow_param_sac(multiplier, err, message)
  use globalData, only: parSubset, nBetaGammaCal
  implicit none
  !input
  type(var_d), intent(in)    :: multiplier(:)   ! mulitpliers for calibrating soil parameter 
  ! output
  integer(i4b),intent(out)   :: err             ! error code
  character(*),intent(out)   :: message         ! error message
  ! Local
  character(len=strLen)      :: parName         ! parameter name
  character(len=strLen)      :: rowfmt
  integer(i4b)               :: iPar,iHru       ! loop index
  real(dp)                   :: param(nHru,22)      ! original snow parameters matrix (nHru x nParamInModel) 
  real(dp)                   :: paramTemp(nHru) ! temporal parameter vector (nHru) 
  integer(i4b)               :: stat,io

  ! initialize error control
  err=0; message='adj_snow_param_sac/'
  !Open original and modified vege parameter files
  open (UNIT=50,file=origvege_name,form='formatted',status='old',IOSTAT=stat)
  open (UNIT=51,file=calivege_name,action='write',status='replace' )
  ! Read original soil parameter file
  do 
    read(unit=50,fmt=*,iostat=io) parName, (paramTemp(iHru), iHru=1,nHru)
    if (io>0) then
      stop 'something wrong in input'
    elseif (io<0) then
      exit
    else
      select case(trim(parName))
        case('scf');   param(:,1)=paramTemp
        case('mfmax'); param(:,2)=paramTemp
        case('mfmin'); param(:,3)=paramTemp
        case('uadj');  param(:,4)=paramTemp
        case('si');    param(:,5)=paramTemp
        case('pxtemp');param(:,6)=paramTemp
        case('nmf');   param(:,7)=paramTemp
        case('tipm');  param(:,8)=paramTemp
        case('mbase'); param(:,9)=paramTemp
        case('plwhc'); param(:,10)=paramTemp
        case('daygm'); param(:,11)=paramTemp
        case('adc1');  param(:,12)=paramTemp
        case('adc2');  param(:,13)=paramTemp
        case('adc3');  param(:,14)=paramTemp
        case('adc4');  param(:,15)=paramTemp
        case('adc5');  param(:,16)=paramTemp
        case('adc6');  param(:,17)=paramTemp
        case('adc7');  param(:,18)=paramTemp
        case('adc8');  param(:,19)=paramTemp
        case('adc9');  param(:,20)=paramTemp
        case('adc10'); param(:,21)=paramTemp
        case('adc11'); param(:,22)=paramTemp
      end select
    endif
  enddo
  do iPar=1,nBetaGammaCal
    select case( parSubset(iPar)%pname )
      case('scf');    param(:,1)=multiplier(iPar)%var(1)*param(:,1) 
      case('mfmax');  param(:,2)=multiplier(iPar)%var(1)*param(:,1) 
      case('mfmin');  param(:,3)=multiplier(iPar)%var(1)*param(:,1) 
      case('uadj');   param(:,4)=multiplier(iPar)%var(1)*param(:,1) 
      case('si');     param(:,5)=multiplier(iPar)%var(1)*param(:,1) 
      case('pxtemp'); param(:,6)=multiplier(iPar)%var(1)*param(:,1) 
      case('nmf');    param(:,7)=multiplier(iPar)%var(1)*param(:,1) 
      case('tipm');   param(:,8)=multiplier(iPar)%var(1)*param(:,1) 
      case('plwhc');  param(:,9)=multiplier(iPar)%var(1)*param(:,1) 
      case('daygm');  param(:,10)=multiplier(iPar)%var(1)*param(:,1) 
    end select
  enddo
  if(    param(iHru, 1) .lt. 0.70_dp) then; param(iHru,1) = 0.70_dp !scf
  elseif(param(iHru, 1) .gt. 1.40_dp) then; param(iHru,1) = 1.40_dp
  elseif(param(iHru, 2) .lt. 0.50_dp) then; param(iHru,2) = 0.50_dp !mfmax 
  elseif(param(iHru, 2) .gt. 2.00_dp) then; param(iHru,2) = 2.00_dp
  elseif(param(iHru, 3) .lt. 0.05_dp) then; param(iHru,3) = 0.05_dp !mfmin
  elseif(param(iHru, 3) .gt. 0.49_dp) then; param(iHru,3) = 0.49_dp
  elseif(param(iHru, 4) .lt. 0.03_dp) then; param(iHru,4) = 0.03_dp !uadj
  elseif(param(iHru, 4) .gt. 0.19_dp) then; param(iHru,4) = 0.19_dp
  elseif(param(iHru, 5) .lt. 0.00_dp) then; param(iHru,5) = 0.00_dp !si
  elseif(param(iHru, 5) .gt. 2000_dp) then; param(iHru,5) = 2000_dp
  elseif(param(iHru, 6) .lt. -2.0_dp) then; param(iHru,6) = -2.0_dp !pxtemp
  elseif(param(iHru, 6) .gt. 2.00_dp) then; param(iHru,6) = 2.00_dp
  elseif(param(iHru, 7) .lt. 0.05_dp) then; param(iHru,7) = 0.05_dp !nmf
  elseif(param(iHru, 7) .gt. 0.50_dp) then; param(iHru,7) = 0.50_dp
  elseif(param(iHru, 8) .lt. 0.10_dp) then; param(iHru,8) = 0.10_dp !tipm
  elseif(param(iHru, 8) .gt. 1.00_dp) then; param(iHru,8) = 1.00_dp
  elseif(param(iHru, 9) .lt. 0.02_dp) then; param(iHru,9) = 0.02_dp !plwhc
  elseif(param(iHru, 9) .gt. 0.30_dp) then; param(iHru,9) = 0.30_dp
  elseif(param(iHru,10) .lt. 0.00_dp) then;param(iHru,10) = 0.00_dp !daygm
  elseif(param(iHru,10) .gt. 0.30_dp) then;param(iHru,10) = 0.30_dp
  endif

  write(rowfmt,'(A,I3,A)') '(A,',nHru,'(1X,F10.5))' 
  write(51,fmt=trim(rowfmt)) ('scf',    param(iHru,1), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('mfmax',  param(iHru,2), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('mfmin',  param(iHru,3), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('uadj',   param(iHru,4), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('si',     param(iHru,5), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('pxtemp', param(iHru,6), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('nmf',    param(iHru,7), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('tipm',   param(iHru,8), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('plwhc',  param(iHru,9), iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('daygm',  param(iHru,10),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc1',   param(iHru,11),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc2',   param(iHru,12),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc3',   param(iHru,13),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc4',   param(iHru,14),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc5',   param(iHru,15),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc6',   param(iHru,16),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc7',   param(iHru,17),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc8',   param(iHru,18),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc9',   param(iHru,19),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc10',  param(iHru,20),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc11',  param(iHru,21),iHru=1,nHru )
  close(UNIT=50)
  close(UNIT=51)
  return
end subroutine

!***************************
! Adjust sac soil parameters with multipliers 
!***************************
subroutine adj_soil_param_sac(param, multiplier, adjParam,  err, message)
!! This routine takes the adjustable parameter set "param" from namelist, reads into "origparam_name",
!! computes the new parameters, writes them into "calibparam_name" 
  use globalData, only: parSubset, nBetaGammaCal
  implicit none
  !input variables
  real(dp),    intent(in)    :: param(:,:)    ! original soil parameters matrix (nHru x nParamInModel) 
  type(var_d), intent(in)    :: multiplier(:) ! mulitpliers for calibrating soil parameter 
  ! output
  real(dp),    intent(out)   :: adjParam(:,:) ! adjusted soil parameter
  integer(i4b),intent(out)   :: err           ! error code
  character(*),intent(out)   :: message       ! error message
  ! local variables
  integer(i4b)               :: iPar,iHru     ! loop index

  ! initialize error control
  err=0; message='adj_soil_param_sac/'
  adjParam=param
  do iHru = 1,nHru
    ! Modify parameter values
    do iPar=1,nBetaGammaCal
      select case( parSubset(iPar)%pname )
        case('twm')
          adjParam(iHru,1) = multiplier( iPar )%var(1)*Param(iHru,1)
          adjParam(iHru,3) = multiplier( iPar )%var(2)*Param(iHru,3)
        case('fwm')
          adjParam(iHru,2) = multiplier( iPar )%var(1)*Param(iHru,2)
        case('fpm')
          adjParam(iHru,4) = multiplier( iPar )%var(1)*Param(iHru,4)
        case('fsm')
          adjParam(iHru,5) = multiplier( iPar )%var(1)*Param(iHru,5)
        case('zk')
          adjParam(iHru,7) = multiplier( iPar )%var(1)*Param(iHru,7)
        case('zpk')
          adjParam(iHru,8) = multiplier( iPar )%var(1)*Param(iHru,8)
        case('zsk')
          adjParam(iHru,9) = multiplier( iPar )%var(1)*Param(iHru,9)
        case('zperc')
          adjParam(iHru,10)= multiplier( iPar )%var(1)*Param(iHru,10)
        case('rexp')
          adjParam(iHru,11)= multiplier( iPar )%var(1)*Param(iHru,11)
        case('pfree')
          adjParam(iHru,13)= multiplier( iPar )%var(1)*Param(iHru,12)
       end select
    end do
    ! Limit parameters to correct possible values without physical meaning: this applies for all configurations
    if(adjParam(iHru,1) .lt. 10.0_dp) then         ! UZTWM
      adjParam(iHru,1) = 10.0_dp
    elseif(adjParam(iHru,1) .gt. 300.0_dp) then
      adjParam(iHru,1) = 300.0_dp
    elseif(adjParam(iHru,2) .lt. 5.0_dp) then      ! UZFWM
      adjParam(iHru,2) = 5.0_dp
    elseif(adjParam(iHru,2) .gt. 150.0_dp) then
      adjParam(iHru,2) = 150.0_dp
    elseif(adjParam(iHru,3) .lt. 10.0_dp) then     ! LZTMM
      adjParam(iHru,3) = 10.0_dp
    elseif(adjParam(iHru,3) .gt. 500.0_dp) then
      adjParam(iHru,3) = 500.0_dp
    elseif(adjParam(iHru,4) .lt. 10.0_dp) then     ! LZFPM
      adjParam(iHru,4) = 10.0_dp
    elseif(adjParam(iHru,4) .gt. 1000.0_dp) then
      adjParam(iHru,4) = 1000.0_dp
    elseif(adjParam(iHru,5) .lt. 5.0_dp) then      ! LZFSM
      adjParam(iHru,5) = 5.0_dp
    elseif(adjParam(iHru,5) .gt. 400.0_dp) then
      adjParam(iHru,5) = 400.0_dp
    elseif(adjParam(iHru,7) .lt. 0.1_dp) then      ! UZK
      adjParam(iHru,7) = 0.1_dp
    elseif(adjParam(iHru,7) .gt. 0.75_dp) then
      adjParam(iHru,7) = 0.75_dp
    elseif(adjParam(iHru,8) .lt. 0.001_dp) then    ! LZPK
      adjParam(iHru,8) = 0.001_dp
    elseif(adjParam(iHru,8) .gt. 0.05_dp) then
      adjParam(iHru,8) = 0.05_dp
    elseif(adjParam(iHru,9) .lt. 0.01_dp) then     ! LZSK
      adjParam(iHru,9) = 0.01_dp
    elseif(adjParam(iHru,9) .gt. 0.35_dp) then
      adjParam(iHru,9) = 0.35_dp
    elseif(adjParam(iHru,10) .lt. 5.0_dp) then     ! ZPERC 
      adjParam(iHru,10) = 5.0_dp
    elseif(adjParam(iHru,10) .gt. 350.0_dp) then
      adjParam(iHru,10) = 350.0_dp
    elseif(adjParam(iHru,11) .lt. 1.0_dp) then     ! rexp 
      adjParam(iHru,11) = 1.0_dp
    elseif(adjParam(iHru,11) .gt. 5.0_dp) then
      adjParam(iHru,11) = 5.0_dp
    elseif(adjParam(iHru,13) .lt. 0.0_dp) then     ! pfree 
      adjParam(iHru,13) = 0.0_dp
    elseif(adjParam(iHru,13) .gt. 0.8_dp) then
      adjParam(iHru,13) = 0.8_dp
    endif
  enddo 
  return
end subroutine

!***************************
! Read sac output file
!***************************
subroutine read_sac_sim(sim, err, message)
  implicit none
  !output variables
  real(dp),              intent(out) :: sim(:,:)
  integer(i4b),          intent(out) :: err            ! error code
  character(*),          intent(out) :: message        ! error message
  !local variables
  character(len=strLen)              :: filename
  real(dp)                           :: cellfraction,basin_area
  real(dp)                           :: auxflux(18)                ! This is only in case of water balance mode
  integer(i4b)                       :: ibasin, itime, ivar, icell ! index 
  integer(i4b)                       :: ncell
  integer(i4b)                       :: dum,c_cell
  character(len=strLen)              :: strDum

  ! initialize error control
  err=0; message='read_sac_sim/'
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
      read (UNIT=55,fmt=*) (strDum, ivar=1,18)
      do itime = 1,sim_len
        read (UNIT=55,fmt=*) (auxflux(ivar), ivar=1,18)
        sim(c_cell,itime) = (auxflux(18))*cellfraction
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

end module sac_routines
