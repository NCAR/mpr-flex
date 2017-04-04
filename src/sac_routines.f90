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
  character(len=strLen)              :: rowfmt_int
  integer(i4b)                       :: iHru         ! loop index

  ! initialize error control
  err=0; message='write_soil_param_sac/'
  if(size(param,2)/=TotNpar)then;err=10;message=trim(message)//"params 2nd dimension size different than TotNpar";return;endif
  if(size(param,1)/=nHru)then;err=11;message=trim(message)//'params 1st dimension size different than nHru';return;endif
  if(size(hruid)/=nHru)then;err=12;message=trim(message)//'hruid size different than nHru';return;endif
  open(UNIT=51,file=trim(calibparam_name),action='write',status='unknown' )
  write(rowfmt,'(A,I3,A)') '(A,',nHru,'(1X,F10.5))' 
  write(rowfmt_int,'(A,I3,A)') '(A,',nHru,'(1X,F8.0))' 
  write(51,fmt=trim(rowfmt_int)) ('hru_id',  param(iHru, 1),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('hru_area',param(iHru, 2),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('uztwm',   param(iHru, 3),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('uzfwm',   param(iHru, 4),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('lztwm',   param(iHru, 5),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('lzfpm',   param(iHru, 6),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('lzfsm',   param(iHru, 7),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('adimp',   param(iHru, 8),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('uzk',     param(iHru, 9),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('lzpk',    param(iHru,10),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('lzsk',    param(iHru,11),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('zperc',   param(iHru,12),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('rexp',    param(iHru,13),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('pctim',   param(iHru,14),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('pfree',   param(iHru,15),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('riva',    param(iHru,16),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('side',    param(iHru,17),iHru=1,nHru )
  write(51,fmt=trim(rowfmt))     ('rserv',   param(iHru,18),iHru=1,nHru )
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
      case('hru_id');   param(:,1)=paramTemp
      case('hru_area'); param(:,2)=paramTemp
      case('uztwm');    param(:,3)=paramTemp
      case('uzfwm');    param(:,4)=paramTemp
      case('lztwm');    param(:,5)=paramTemp
      case('lzfpm');    param(:,6)=paramTemp
      case('lzfsm');    param(:,7)=paramTemp
      case('adimp');    param(:,8)=paramTemp
      case('uzk');      param(:,9)=paramTemp
      case('lzpk');     param(:,10)=paramTemp
      case('lzsk');     param(:,11)=paramTemp
      case('zperc');    param(:,12)=paramTemp
      case('rexp');     param(:,13)=paramTemp
      case('pctim');    param(:,14)=paramTemp
      case('pfree');    param(:,15)=paramTemp
      case('riva');     param(:,16)=paramTemp
      case('side');     param(:,17)=paramTemp
      case('rserv');    param(:,18)=paramTemp
    end select
  end do
  close(UNIT=50)
  return
end subroutine

!***************************
! replace sac soil parameters 
!***************************
subroutine replace_soil_param_sac(param, parMxyMz, adjParam, err, message)
  use globalData, only: betaMeta, calBetaName, nSoilBetaModel 
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
    do iPar=1,nSoilBetaModel
      associate( ix=>get_ixBeta(trim(calBetaName(iPar))) )
      select case( betaMeta(ix)%pname )
        case('twm')
          adjParam(iHru,3)  = parMxyMz(iPar)%varData(1,iHru)
          adjParam(iHru,5)  = parMxyMz(iPar)%varData(nLyr,iHru)
        case('fwm')
          adjParam(iHru,4)  = parMxyMz(iPar)%varData(1,iHru)
        case('fpm')
          adjParam(iHru,6)  = parMxyMz(iPar)%varData(nLyr,iHru)
        case('fsm')
          adjParam(iHru,7)  = parMxyMz(iPar)%varData(nLyr,iHru)
        case('zk')
          adjParam(iHru,9)  = parMxyMz(iPar)%varData(1,iHru)
        case('zpk')
          adjParam(iHru,10) = parMxyMz(iPar)%varData(nLyr,iHru)
        case('zsk')
          adjParam(iHru,11) = parMxyMz(iPar)%varData(nLyr,iHru)
        case('zperc')
          adjParam(iHru,12) = parMxyMz(iPar)%varData(1,iHru)
        case('rexp')
          adjParam(iHru,13) = parMxyMz(iPar)%varData(1,iHru)
        case('pfree')
          adjParam(iHru,15) = parMxyMz(iPar)%varData(1,iHru)
      end select
      end associate
    end do
  enddo hru
  return
end subroutine

subroutine adj_snow_param_sac(multiplier, err, message)
  use globalData, only: calParMeta, nCalPar
  implicit none
  !input
  type(var_d), intent(in)    :: multiplier(:)   ! mulitpliers for calibrating soil parameter 
  ! output
  integer(i4b),intent(out)   :: err             ! error code
  character(*),intent(out)   :: message         ! error message
  ! Local
  character(len=strLen)      :: parName         ! parameter name
  character(len=strLen)      :: rowfmt
  character(len=strLen)      :: rowfmt_int
  integer(i4b)               :: iPar,iHru       ! loop index
  real(dp)                   :: param(nHru,25)  ! original snow parameters matrix (nHru x nParamInModel) 
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
        case('hru_id');  param(:,1)=paramTemp
        case('latitude');param(:,2)=paramTemp
        case('elev');    param(:,3)=paramTemp
        case('scf');     param(:,4)=paramTemp
        case('mfmax');   param(:,5)=paramTemp
        case('mfmin');   param(:,6)=paramTemp
        case('uadj');    param(:,7)=paramTemp
        case('si');      param(:,8)=paramTemp
        case('pxtemp');  param(:,9)=paramTemp
        case('nmf');     param(:,10)=paramTemp
        case('tipm');    param(:,11)=paramTemp
        case('mbase');   param(:,12)=paramTemp
        case('plwhc');   param(:,13)=paramTemp
        case('daygm');   param(:,14)=paramTemp
        case('adc1');    param(:,15)=paramTemp
        case('adc2');    param(:,16)=paramTemp
        case('adc3');    param(:,17)=paramTemp
        case('adc4');    param(:,18)=paramTemp
        case('adc5');    param(:,19)=paramTemp
        case('adc6');    param(:,20)=paramTemp
        case('adc7');    param(:,21)=paramTemp
        case('adc8');    param(:,22)=paramTemp
        case('adc9');    param(:,23)=paramTemp
        case('adc10');   param(:,24)=paramTemp
        case('adc11');   param(:,25)=paramTemp
      end select
    endif
  enddo
  do iPar=1,nCalPar
    select case( calParMeta(iPar)%pname )
      case('scf');    param(:,4)=multiplier(iPar)%var(1)*param(:,4) 
      case('mfmax');  param(:,5)=multiplier(iPar)%var(1)*param(:,5) 
      case('mfmin');  param(:,6)=multiplier(iPar)%var(1)*param(:,6) 
      case('uadj');   param(:,7)=multiplier(iPar)%var(1)*param(:,7) 
      case('si');     param(:,8)=multiplier(iPar)%var(1)*param(:,8) 
      case('pxtemp'); param(:,9)=multiplier(iPar)%var(1)*param(:,9) 
      case('nmf');    param(:,10)=multiplier(iPar)%var(1)*param(:,11) 
      case('tipm');   param(:,11)=multiplier(iPar)%var(1)*param(:,12) 
      case('plwhc');  param(:,13)=multiplier(iPar)%var(1)*param(:,13) 
      case('daygm');  param(:,14)=multiplier(iPar)%var(1)*param(:,14) 
    end select
  enddo
  where (param(:, 4) .lt. 0.70_dp) param(:,4) = 0.70_dp !scf
  where (param(:, 4) .gt. 1.40_dp) param(:,4) = 1.40_dp
  where (param(:, 5) .lt. 0.50_dp) param(:,5) = 0.50_dp !mfmax 
  where (param(:, 5) .gt. 2.00_dp) param(:,5) = 2.00_dp
  where (param(:, 6) .lt. 0.05_dp) param(:,6) = 0.05_dp !mfmin
  where (param(:, 6) .gt. 0.49_dp) param(:,6) = 0.49_dp
  where (param(:, 7) .lt. 0.03_dp) param(:,7) = 0.03_dp !uadj
  where (param(:, 7) .gt. 0.19_dp) param(:,7) = 0.19_dp
  where (param(:, 8) .lt. 0.00_dp) param(:,8) = 0.00_dp !si
  where (param(:, 8) .gt. 2000_dp) param(:,8) = 2000_dp
  where (param(:, 9) .lt. -2.0_dp) param(:,9) = -2.0_dp !pxtemp
  where (param(:, 9) .gt. 2.00_dp) param(:,9) = 2.00_dp
  where (param(:,10) .lt. 0.05_dp) param(:,10) = 0.05_dp !nmf
  where (param(:,10) .gt. 0.50_dp) param(:,10) = 0.50_dp
  where (param(:,11) .lt. 0.10_dp) param(:,11) = 0.10_dp !tipm
  where (param(:,11) .gt. 1.00_dp) param(:,11) = 1.00_dp
  where (param(:,13) .lt. 0.02_dp) param(:,13) = 0.02_dp !plwhc
  where (param(:,13) .gt. 0.30_dp) param(:,13) = 0.30_dp
  where (param(:,14) .lt. 0.00_dp) param(:,14) = 0.00_dp !daygm
  where (param(:,14) .gt. 0.30_dp) param(:,14) = 0.30_dp

  write(rowfmt,'(A,I3,A)') '(A,',nHru,'(1X,F10.5))' 
  write(rowfmt_int,'(A,I3,A)') '(A,',nHru,'(1X,F8.0))' 
  write(51,fmt=trim(rowfmt_int)) ('hru_id',  param(iHru, 1),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('latitude',param(iHru, 2),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('elev',    param(iHru, 3),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('scf',     param(iHru, 4),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('mfmax',   param(iHru, 5),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('mfmin',   param(iHru, 6),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('uadj',    param(iHru, 7),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('si',      param(iHru, 8),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('pxtemp',  param(iHru, 9),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('nmf',     param(iHru,10),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('tipm',    param(iHru,11),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('mbase',   param(iHru,12),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('plwhc',   param(iHru,13),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('daygm',   param(iHru,14),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc1',    param(iHru,15),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc2',    param(iHru,16),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc3',    param(iHru,17),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc4',    param(iHru,18),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc5',    param(iHru,19),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc6',    param(iHru,20),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc7',    param(iHru,21),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc8',    param(iHru,22),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc9',    param(iHru,23),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc10',   param(iHru,24),iHru=1,nHru )
  write(51,fmt=trim(rowfmt)) ('adc11',   param(iHru,25),iHru=1,nHru )
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
  use globalData, only: calParMeta, nCalPar
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
    do iPar=1,nCalPar
      select case( calParMeta(iPar)%pname )
        case('twm')
          adjParam(iHru,3) = multiplier( iPar )%var(1)*Param(iHru,3)
          adjParam(iHru,5) = multiplier( iPar )%var(2)*Param(iHru,5)
        case('fwm')
          adjParam(iHru,4) = multiplier( iPar )%var(1)*Param(iHru,4)
        case('fpm')
          adjParam(iHru,6) = multiplier( iPar )%var(1)*Param(iHru,6)
        case('fsm')
          adjParam(iHru,7) = multiplier( iPar )%var(1)*Param(iHru,7)
        case('zk')
          adjParam(iHru,9) = multiplier( iPar )%var(1)*Param(iHru,9)
        case('zpk')
          adjParam(iHru,10) = multiplier( iPar )%var(1)*Param(iHru,10)
        case('zsk')
          adjParam(iHru,11) = multiplier( iPar )%var(1)*Param(iHru,11)
        case('zperc')
          adjParam(iHru,12)= multiplier( iPar )%var(1)*Param(iHru,12)
        case('rexp')
          adjParam(iHru,13)= multiplier( iPar )%var(1)*Param(iHru,13)
        case('pfree')
          adjParam(iHru,15)= multiplier( iPar )%var(1)*Param(iHru,15)
       end select
    end do
    ! Limit parameters to correct possible values without physical meaning: this applies for all configurations
    if(adjParam(iHru,3) .lt. 10.0_dp) then         ! UZTWM
      adjParam(iHru,3) = 10.0_dp
    elseif(adjParam(iHru,3) .gt. 300.0_dp) then
      adjParam(iHru,3) = 300.0_dp
    elseif(adjParam(iHru,4) .lt. 5.0_dp) then      ! UZFWM
      adjParam(iHru,4) = 5.0_dp
    elseif(adjParam(iHru,4) .gt. 150.0_dp) then
      adjParam(iHru,4) = 150.0_dp
    elseif(adjParam(iHru,5) .lt. 10.0_dp) then     ! LZTMM
      adjParam(iHru,5) = 10.0_dp
    elseif(adjParam(iHru,5) .gt. 500.0_dp) then
      adjParam(iHru,5) = 500.0_dp
    elseif(adjParam(iHru,6) .lt. 10.0_dp) then     ! LZFPM
      adjParam(iHru,6) = 10.0_dp
    elseif(adjParam(iHru,6) .gt. 1000.0_dp) then
      adjParam(iHru,6) = 1000.0_dp
    elseif(adjParam(iHru,7) .lt. 5.0_dp) then      ! LZFSM
      adjParam(iHru,7) = 5.0_dp
    elseif(adjParam(iHru,7) .gt. 400.0_dp) then
      adjParam(iHru,7) = 400.0_dp
    elseif(adjParam(iHru,9) .lt. 0.1_dp) then      ! UZK
      adjParam(iHru,9) = 0.1_dp
    elseif(adjParam(iHru,7) .gt. 0.75_dp) then
      adjParam(iHru,9) = 0.75_dp
    elseif(adjParam(iHru,10) .lt. 0.001_dp) then    ! LZPK
      adjParam(iHru,10) = 0.001_dp
    elseif(adjParam(iHru,10) .gt. 0.05_dp) then
      adjParam(iHru,10) = 0.05_dp
    elseif(adjParam(iHru,11) .lt. 0.01_dp) then     ! LZSK
      adjParam(iHru,11) = 0.01_dp
    elseif(adjParam(iHru,11) .gt. 0.35_dp) then
      adjParam(iHru,11) = 0.35_dp
    elseif(adjParam(iHru,12) .lt. 5.0_dp) then     ! ZPERC 
      adjParam(iHru,12) = 5.0_dp
    elseif(adjParam(iHru,12) .gt. 350.0_dp) then
      adjParam(iHru,12) = 350.0_dp
    elseif(adjParam(iHru,13) .lt. 1.0_dp) then     ! rexp 
      adjParam(iHru,13) = 1.0_dp
    elseif(adjParam(iHru,13) .gt. 5.0_dp) then
      adjParam(iHru,13) = 5.0_dp
    elseif(adjParam(iHru,15) .lt. 0.0_dp) then     ! pfree 
      adjParam(iHru,15) = 0.0_dp
    elseif(adjParam(iHru,15) .gt. 0.8_dp) then
      adjParam(iHru,15) = 0.8_dp
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
