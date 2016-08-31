module vic_routines

  use nrtype 
  use public_var
  use strings
  use data_type 

  implicit none
  public :: vic_soil_param 

contains

subroutine vic_soil_param(param, err, message)
!! This routine takes the adjustable parameter set "param" from namelist, reads into "origparam_name",
!! computes the new parameters, writes them into "calibparam_name" 

  use globalData, only: parSubset

  implicit none

!input variables
  real(dp),dimension(:),intent(in)    :: param            ! parameter in namelist, not necessarily all parameters are calibrated
! output
  integer(i4b),intent(out)             :: err            ! error code
  character(*),intent(out)             :: message        ! error message
!local variables
  integer(i4b)                        :: ipar,iHru    ! loop index
  integer(i4b)                        :: stat
  real(dp),dimension(TotNparVic)      :: realline

  ! initialize error control
  err=0; message='vic_soil_param/'

 !Open original and modified basin parameter files
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
  open (UNIT=51,file=calibparam_name,action='write',status='unknown' )

 ! Read original soil parameter file
  do iHru = 1,Ncells
    read(unit=50,*) (realline(ipar), ipar=1,TotNparVic)
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
        case('bbl');      realline(28:30) = param( iPar )*realline(28:30)
        case('BD');       realline(34:36) = param( iPar )*realline(34:36)
        case('SD');       realline(37:39) = param( iPar )*realline(37:39)
        case('WcrFrac');  realline(41:43) = param( iPar )*realline(41:43)
        case('WpwpFrac'); realline(44:46) = param( iPar )*realline(44:46)
        case default; err=10; message=trim(message)//'parameter name not found';return 
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
      print*, realline(7)
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

end subroutine vic_soil_param

end module vic_routines
