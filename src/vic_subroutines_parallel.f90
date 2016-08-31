module vic_subroutines_parallel

  use nrtype 
  use public_var
  use strings
  use data_type 
  use uh,            only: duamel

  implicit none
  public :: eval_objfn

contains

! -----------------------------------------
subroutine eval_objfn(param,objfnval)
! -----------------------------------------

  implicit none

!! This routine takes the adjustable parameter set "param" from namelist, reads into "origparam_name",
!! computes the new parameters, writes them into "calibparam_name" and finally
!! runs the model

!input variables
  real(dp),dimension(:),intent(in)    :: param            ! parameter in namelist, not necessarily all parameters are calibrated
!output variables
  real(dp),             intent(out)   :: objfnval 
!local variables
  character(len=strLen)               :: rdline,filename
  integer(i4b)                        :: ipar,jcell    ! loop index
  integer(i4b)                        :: i,ntau
  integer(i4b)                        :: stat
  real(dp),dimension(TotNparVic)      :: realline
  real(dp),dimension(:)  ,allocatable :: route_qmod,qmod,qobs
  real(dp),dimension(:,:),allocatable :: route_qmod_region,qmod_region
  real(dp),dimension(:,:),allocatable :: route_qmod_region_cell,qmod_region_cell
  real(dp)                            :: dtuh
  real(dp)                            :: total_depth  !total depth of soil column
  real(dp)                            :: ushape,uscale
  !integer(i8b)                         :: t1, t2, count_rate,j

 !Open original and modified basin parameter files
  open (UNIT=50,file=origparam_name,form='formatted',status='old',IOSTAT=stat)
  open (UNIT=51,file=calibparam_name,form='formatted',status='old')

  !allocate qobs and qsim 
  if(region_flag .eq. 0) then
    ALLOCATE (qobs(sim_len))
    ALLOCATE (qmod(sim_len))
    !! Initialize variables
    qmod = 0.0
    qobs = 0.0
  elseif(region_flag .eq. 1) then
    ALLOCATE(qmod_region_cell(ntot,sim_len))
    ALLOCATE(qmod_region(nbasin,sim_len))
    ALLOCATE(qobs(nbasin*sim_len))
    !! Initialize variables
    qobs = 0.0
    qmod_region = 0.0
    qmod_region_cell = 0.0
  endif

  !allocate routed flow array
  if(region_flag .eq. 0) then
    allocate(route_qmod(size(qobs)))
    route_qmod = 0.0
  elseif(region_flag .eq. 1) then
    allocate(route_qmod_region_cell(ntot,sim_len))
    allocate(route_qmod_region(nbasin,sim_len))
    route_qmod_region = 0.0
    route_qmod_region_cell = 0.0
  endif

 ! Read original soil parameter file
  do jcell = 1,Ncells
    read(unit=50,*) (realline(ipar), ipar=1,TotNparVic)
    !apply multipliers directly to VIC parameters if running traditional upscaling 
    if(upscale_flag .eq. 0) then
      !B_infilt
      realline(5) = param(3)*realline(5)
      !Ds
      realline(6) = param(4)*realline(6)
      !Ds_max
      realline(7) = param(5)*realline(7)
      !Ws
      realline(8) = param(6)*realline(8)
      !Ksat
      realline(13:15) = param(10)*realline(13:15)
      !total depth of soil column
      total_depth = param(9) * (realline(23)+realline(24)+realline(25))
      !depth of 1st layer
      realline(23) = param(7)*total_depth
      !depth of 2nd layer
      realline(24) = param(8)*total_depth
      !depth of 3rd layer
      realline(25) = total_depth - (realline(23)+realline(24))
      !bulk density
      realline(34:36) = param(11)*realline(34:36)
    end if

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
    do i = 34,36
      if(realline(i) .lt. 805.) then
        realline(i) = 805.
      elseif(realline(i) .gt. 1880.) then
        realline(i) = 1880.
      endif
    enddo
    ! Write the modified parameter file for the entire basin/region for traditional upscaling
    if(upscale_flag .eq. 0) then
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
    end if
  enddo  !end cell loop

  ! Close original and modified basin parameter files
  close(UNIT=50)
  close(UNIT=51)

! if upscale flag is selected, run MPR to rebuild VIC paramters from statsgo and transfer functions
  !don't need to run this for each cell, only one time since these apply globally
  if(upscale_flag .ne. 0 .or. upscale_flag .eq. 0) then
    !want to write out hydro params to update the hydro_params namelist with
    open(unit=52,file='./hydro_pedotransfer.txt',form='formatted')
    !Binfilt 
    write(unit=52,*) 'a_infilt',param(3)
    !ds (D1) 
    if(param(4) .ne. 0) then
      write(unit=52,*) 'a_Ds',param(4)
    else
      write(unit=52,*) 'a_Ds',1.0
    endif
    ! Dsmx (D2)
    if(param(5) .ne. 0) then
      write(unit=52,*) 'a_Dsmax',param(5)
    else
      write(unit=52,*) 'a_Dsmax',1.0
    endif
    !Ws (D3) 
    if(param(6) .ne. 0) then
      write(unit=52,*) 'a_Ws',param(6)
    else
      write(unit=52,*) 'a_Ws',1.0
    endif
    !depth1
    if(param(7) .ne. 0) then
      write(unit=52,*) 'a_d1',param(7)
    else
      write(unit=52,*) 'a_d1',0.05
    endif
    !depth2
    if(param(8) .ne. 0) then
      write(unit=52,*) 'a_d2',param(8)
    else
      write(unit=52,*) 'a_d2',0.3
    endif
    !depth tot
    if(param(9) .ne. 0) then
      write(unit=52,*) 'a_h1',param(9)
      write(unit=52,*) 'a_h2',param(9)
      write(unit=52,*) 'a_h3',param(9)
      write(unit=52,*) 'a_h4',param(9)
      write(unit=52,*) 'a_h5',param(9)
      write(unit=52,*) 'a_h6',param(9)
      write(unit=52,*) 'a_h7',param(9)
      write(unit=52,*) 'a_h8',param(9)
      write(unit=52,*) 'a_h9',param(9)
      write(unit=52,*) 'a_h10',param(9)
      write(unit=52,*) 'a_h11',param(9)
    else
      write(unit=52,*) 'a_h1',1.0
      write(unit=52,*) 'a_h2',1.0
      write(unit=52,*) 'a_h3',1.0
      write(unit=52,*) 'a_h4',1.0
      write(unit=52,*) 'a_h5',1.0
      write(unit=52,*) 'a_h6',1.0
      write(unit=52,*) 'a_h7',1.0
      write(unit=52,*) 'a_h8',1.0
      write(unit=52,*) 'a_h9',1.0
      write(unit=52,*) 'a_h10',1.0
      write(unit=52,*) 'a_h11',1.0
    endif
    !Ksat
    if(param(10) .ne. 0) then
      write(unit=52,*) 'a_Ks',param(10)
    else
      write(unit=52,*) 'a_Ks',1.0
    endif
    !Bulk density
    if(param(11) .ne. 0) then
      write(unit=52,*) 'a_bulkDensity',param(11)
    else
      write(unit=52,*) 'a_bulkDensity',1.0
    endif
    !min. stomatal
    if(param(12) .ne. 0) then
      write(unit=52,*) 'a_rsmin',param(12)
    else
      write(unit=52,*) 'a_rsmin',1.0
    endif
    ! Lai
    if(param(13) .ne. 0) then
      write(unit=52,'(A,1x,F8.6)') 'a_lai',param(13)
    else
      write(unit=52,'(A,1x,F8.6)') 'a_lai',1.0
    endif
    ! new snow albedo
    if(param(14) .ne. 0) then
      write(unit=52,'(A,1x,F8.6)') 'a_new_snow',param(14)
    else
      write(unit=52,'(A,1x,F8.6)') 'a_new_snow',1.0
    endif
    ! albedo decay coef. during accumulation
    if(param(15) .ne. 0) then
      write(unit=52,'(A,1x,F8.6)') 'a_accum_a',param(15)
    else
      write(unit=52,'(A,1x,F8.6)') 'a_accum_a',1.0
    endif
    ! albedo decay coef. during melting  
    if(param(16) .ne. 0) then
      write(unit=52,'(A,1x,F8.6)') 'a_melt_a',param(16)
    else
      write(unit=52,'(A,1x,F8.6)') 'a_melt_a',1.0
    endif
    ! rain/snow threshold (mid-point)
    if(param(17) .ne. 0) then
      write(unit=52,'(A,1x,F8.6)') 'a_rain_snow',param(17)
    else
      write(unit=52,'(A,1x,F8.6)') 'a_rain_snow',1.0
    endif

    close(unit=52)

  end if

  call system(executable) !< run VIC

 if(region_flag .eq. 0) then
    !read observations
    call read_obsflow(outputobs_name,sim_len,qobs)
    !read simulations 
    call read_modflow(outputmod_name,cellfrac_name,filelist_name,sim_len,qmod)
  ! route flow now
    dtuh = real(86400./86400.) 
    ntau = 0
    if (param(1) .le. 0.0 .and. param(2) .le. 0.0) THEN
      route_qmod = qmod
    else
      call duamel(qmod, param(1), param(2),dtuh,size(qmod,1)-1,route_qmod,ntau)
    end if
  ! calculate objective function
    objfnval = 0.0
    call calc_rmse(route_qmod,qobs,objfnval)
  elseif(region_flag .eq. 1) then
    !read observations
    call read_obsflow_region(outputobs_name,sim_len,nbasin,qobs)
    !read simulations 
    call read_modflow_region(outputmod_name,cellfrac_name,filelist_name,region_info,sim_len,qmod_region)
    !read model output  place into array of ncells 
    call read_vic_cells(outputmod_name,cellfrac_name,filelist_name,region_info,sim_len,qmod_region_cell)
    !call function to route flow for each basin
    ushape = param(1)
    uscale = param(2)
    !call function to perform UH on every grid cell
    call route_grid_cell(qmod_region_cell,route_qmod_region_cell,ushape,uscale)
    !aggregate UH routed grid cell runoff to basin total runoff
    call agg_cells_to_basin(region_info,route_qmod_region_cell,route_qmod_region)
    !call rmse calculation
    objfnval = 0.0
    call calc_rmse_region(route_qmod_region,qobs,objfnval)
  endif

  if(region_flag .eq. 0) then
    deallocate(qobs)
    deallocate(qmod)
    deallocate(route_qmod)
  elseif(region_flag .eq. 1) then
    deallocate(qobs)
    deallocate(qmod_region)
    deallocate(route_qmod_region)
    deallocate(route_qmod_region_cell)
    deallocate(qmod_region_cell)
  endif

  return

end subroutine eval_objfn

! Calculate Objective functions
! -------------------------------------------
subroutine calc_rmse(model,streamflow,rmse)
! -------------------------------------------

  implicit none

!input variables (model: simulations, streamflow: observations)
  real(dp),dimension(:),intent(in)    :: model
  real(dp),dimension(:),intent(in)    :: streamflow
!output variables
  real(dp), intent(out)               :: rmse
!local variables
  integer                             :: itime
  real(dp)                            :: sum_sqr
  integer                             :: nargs
  character(len=strLen)               :: out_name
  character(len=strLen),dimension(10) :: tokens
  character(len=strLen)               :: last_token
  character(len=1)                    :: delims
  real(dp),allocatable,dimension(:)   :: log_model
  real(dp),allocatable,dimension(:)   :: log_streamflow
  integer                             :: nmonths  !for monthly rmse calculation
  integer                             :: start_ind, end_ind
  real                                :: month_model, month_streamflow

  sum_sqr = 0.0

  allocate(log_model(sim_len))
  allocate(log_streamflow(sim_len))

!want to output basin total streamflow to a different file of opt .ne. 1
  delims='/'
  if(opt .ne. 1) then
    call parse(outputobs_name,delims,tokens,nargs)
    last_token = tokens(nargs)

    out_name = trim(outputmod_name)//last_token(1:9)//"_flow.txt"
    print *,trim(out_name)
    open(unit=88,file=out_name)
  endif

  !calculate rmse using log(flow)
!  log_streamflow = streamflow + 0.0001
!  log_model      = model + 0.0001
  log_streamflow = streamflow
  log_model      = model
!  where(log_model .eq. 0.0) log_model = 0.001
!  where(log_streamflow .eq. 0.0) log_streamflow = 0.001

!  log_model = log(log_model)
!  log_streamflow = log(log_streamflow)

  select case (trim(eval_length))
    case ("daily","Daily","DAILY")
      sum_sqr = 0.0
      do itime = start_cal,end_cal-1
        sum_sqr = sum_sqr + (log_model(itime)-log_streamflow(itime))**2
        !print basin total flow to file
        if(opt .ne. 1) then
          write(unit=88,*) model(itime),streamflow(itime)
        endif
      enddo
      rmse = sqrt(sum_sqr/real((end_cal-1) - start_cal + 1))
    case ("monthly","Monthly","MONTHLY")
      sum_sqr = 0.0
      nmonths = floor(((end_cal-start_cal)+1)/30.0)  !use 30 day months uniformly to make it easier
      start_ind = start_cal
      end_ind = start_cal + 29
      do itime = 1,nmonths
        month_model = sum(log_model(start_ind:end_ind))
        month_streamflow = sum(log_streamflow(start_ind:end_ind))
        sum_sqr = sum_sqr + (month_model-month_streamflow)**2
        start_ind = end_ind+1
        end_ind = end_ind + 29
      enddo
      !grab remainder portion and weight it by number of days
      month_model = sum(log_model(end_ind+1:end_cal-1))
      month_streamflow = sum(log_streamflow(end_ind+1:end_cal-1))
      sum_sqr = sum_sqr + ((month_model-month_streamflow)**2) * real((end_cal-1-(end_ind+1))/30.0)
      rmse = sqrt(sum_sqr/real(nmonths+1))
    case default
      sum_sqr = 0.0
      do itime = start_cal,end_cal-1
        sum_sqr = sum_sqr + (log_model(itime)-log_streamflow(itime))**2
        ! print basin total flow to file
        if(opt .ne. 1) then
          write(unit=88,*) model(itime),streamflow(itime)
        endif
      enddo
      rmse = sqrt(sum_sqr/real((end_cal-1) - start_cal + 1))
  end select

  if(opt .ne. 1) then
    close(unit=88)
  endif

  return
end subroutine calc_rmse

! -----------------------------------
subroutine calc_nse(qsim,qobs,nse)

  implicit none

!input variables (Qsim, qobs, length of evaluation time period at daily)
  real(dp), dimension(:), intent(in) :: qsim
  real(dp), dimension(:), intent(in) :: qobs
!output variables 
  real(dp), intent(out)              :: nse 
!local variables
  integer                            :: itime      ! index of time looping
  real(dp)                           :: sumSqrErr  ! sum of squre of Qsim minus Qobs for entire time period
  real(dp)                           :: sumSqrDev  ! sum of squre of deviation of obs from mean for entire time period
  real(dp)                           :: sumQ       ! sum of Qobs for entire time period
  real(dp)                           :: meanQ      ! mean of Qobs for entire time period
  integer                            :: nargs
  character(len=2000)                :: out_name
  character(len=1000),dimension(10)  :: tokens
  character(len=1000)                :: last_token
  character(len=1)                   :: delims
  integer(i4b)                        :: nmonths    ! number of months within evaluation time period for monthly nse calculation
  integer(i4b)                        :: start_ind  ! inex of first day in a month ( monthly loop) 
  integer(i4b)                        :: end_ind    ! inex of last day in a month ( monthly loop) 
  real(dp)                           :: month_qsim, month_qobs
  
  ! initinalize variabes
  sumSqrErr = 0.0
  sumSqrDev = 0.0
  sumQ      = 0.0

  ! output basin total qsim and qobs to a different file if opt .ne. 1
  delims='/'
  if(opt .ne. 1) then
    call parse(outputobs_name,delims,tokens,nargs)
    last_token = tokens(nargs)

    out_name = trim(outputmod_name)//last_token(1:9)//"_flow.txt"
    print *,trim(out_name)
    open(unit=88,file=out_name)
    do itime = start_cal,end_cal-1
      if(opt .ne. 1) then
        write(unit=88,*) qsim(itime),qobs(itime)
      endif
    enddo
  endif
  
  select case (trim(eval_length))
    case ("daily","Daily","DAILY") 
      ! Compute mean of obs
      sumQ = sum(qobs(start_cal:end_cal-1))
      meanQ=sumQ/real((end_cal-1)-start_cal+1)
      ! Compute sum of square of error and std of obs 
      do itime = start_cal,end_cal-1
        sumSqrErr = sumSqrErr + (qsim(itime)-qobs(itime))**2
        sumSqrDev = sumSqrDev + (qobs(itime)-meanQ)**2
      enddo
      nse = sumSqrErr/sumSqrDev

    case ("monthly","Monthly","MONTHLY")
      nmonths = floor(((end_cal-start_cal)+1)/30.0)  !use 30 day months uniformly to make it easier
      ! Compute montly observed Q
      ! Indices of start and end for first month
      start_ind = start_cal
      end_ind   = start_cal + 29
      do itime = 1,nmonths
        month_qobs = sum(qobs(start_ind:end_ind))
        sumQ       = sumQ + month_qobs
        ! Update starting and ending indice for next month step
        start_ind = end_ind+1
        end_ind   = end_ind + 29
      enddo
      !grab remainder portion and weight it by number of days
      sumQ  = sumQ + (month_qobs) * real((end_cal-1-(end_ind+1))/30.0)
      meanQ = sumQ /real(nmonths+1)
       
      ! Compute  sum of square of monthly error and std of monthly obs
      ! Indices of start and end for first month
      start_ind = start_cal
      end_ind   = start_cal + 29
      do itime = 1,nmonths
        month_qsim = sum(qsim(start_ind:end_ind))
        month_qobs = sum(qobs(start_ind:end_ind))
        sumSqrErr  = sumSqrErr + (month_qsim-month_qobs)**2
        sumSqrDev  = sumSqrDev + (month_qobs-meanQ)**2
        ! Update starting and ending indice for next month step
        start_ind = end_ind+1
        end_ind   = end_ind + 29
      enddo
      !grab remainder portion and weight it by number of days
      month_qsim = sum(qsim(end_ind+1:end_cal-1))
      month_qobs = sum(qobs(end_ind+1:end_cal-1))
      sumSqrErr = sumSqrErr + ((month_qsim-month_qobs)**2) * real((end_cal-1-(end_ind+1))/30.0)
      sumSqrDev = sumSqrDev + ((month_qobs-meanQ)**2) * real((end_cal-1-(end_ind+1))/30.0)

      nse = sumSqrErr/sumSqrDev

    case default  ! default is daily
      ! Compute mean of obs
      do itime = start_cal,end_cal-1
        sumQ = sumQ + qobs(itime)
      enddo
      meanQ=sumQ/real((end_cal-1)-start_cal+1)
      ! Compute sum of square of error and std of obs 
      do itime = start_cal,end_cal-1
        sumSqrErr = sumSqrErr + (qsim(itime)-qobs(itime))**2
        sumSqrDev = sumSqrDev + (qobs(itime)-meanQ)**2
      enddo
      nse = sumSqrErr/sumSqrDev
  end select

  if(opt .ne. 1) then
    close(unit=88)
  endif

  return
end subroutine calc_nse

!--------------------------------------------------
subroutine calc_rmse_region(model,streamflow,rmse)

  implicit none

!input variables (model: simulations, streamflow: observations)
  real(dp), dimension(:,:), intent(in)  :: model
  real(dp), dimension(:),   intent(in)  :: streamflow
!output variables
  real(dp),                 intent(out) :: rmse
!local variables
  integer                               :: itime,ibasin,total_len,offset
  real(dp)                              :: sum_sqr
  integer(i4b)                           :: nargs,nstream,nb,nday
  character(len=strLen)                 :: out_name
  character(len=strlen),dimension(10)   :: tokens
  character(len=strlen)                 :: last_token
  character(len=1)                      :: delims
  real(dp),allocatable,dimension(:,:)   :: log_model
  real(dp),allocatable,dimension(:)     :: log_streamflow
  integer,allocatable,dimension(:)      :: basin_id
  real(dp),allocatable,dimension(:)     :: obj_fun_weight
  real(dp),allocatable,dimension(:)     :: basin_rmse
  integer(i4b)                           :: nmonths  !for monthly rmse calculation
  integer(i4b)                           :: start_ind, end_ind
  integer(i4b)                           :: start_obs, end_obs
  integer(i4b)                           :: rmse_period
  real(dp)                              :: month_model, month_streamflow

  total_len = (end_cal-start_cal)*nbasin
  sum_sqr = 0.0

  allocate(log_model(nbasin,sim_len))
  allocate(log_streamflow(nbasin*sim_len))
  allocate(obj_fun_weight(nbasin))
  allocate(basin_rmse(nbasin))
  allocate(basin_id(nbasin))

!want to output streamflow to a different file if opt .ne. 1
  delims='/'
  if(opt .ne. 1) then
    call parse(outputobs_name,delims,tokens,nargs)
    last_token = tokens(nargs)

    out_name = trim(outputmod_name)//last_token(1:9)//"_flow.txt"
    print *,trim(out_name)
    open(unit=88,file=out_name)
  endif

!read in basin weight file
!this file determines how much each basin contributes to the total rmse
!weights need to sum to 1 in the file
  open (UNIT=58,file=trim(basin_objfun_weight_file),form='formatted',status='old')

  do ibasin = 1,nbasin
    read (UNIT=58,*) basin_id(ibasin),obj_fun_weight(ibasin)
  enddo
  close(UNIT=58)

  log_streamflow = streamflow
  log_model      = model

  !need to make sure i'm using the appropriate parts of the catenated region observed streamflow timeseries
  do ibasin = 0,nbasin-1
    !offset places me at the start of each basin
    offset = ibasin*sim_len

    !then run from the starting calibration point to the ending calibration point
    select case (trim(eval_length))
      case ("daily","Daily","DAILY")
        sum_sqr = 0.0
        do itime = start_cal,end_cal-1
          sum_sqr = sum_sqr + ((log_model(ibasin+1,itime)-log_streamflow(itime+offset))**2)
          if(opt .ne. 1) then
            write(unit=88,*) model(ibasin+1,itime),streamflow(itime+offset)
          endif
        enddo
        basin_rmse(ibasin+1) = sqrt(sum_sqr/real((end_cal-start_cal)+1))
      case ("monthly","Monthly","MONTHLY")
        sum_sqr = 0.0
        nmonths = floor(((end_cal-start_cal)+1)/30.0)  !use 30 day months uniformly to make it easier
        start_ind = start_cal
        end_ind = start_cal + 29
        do itime = 1,nmonths
          month_model = sum(log_model(ibasin+1,start_ind:end_ind))
          month_streamflow = sum(log_streamflow(offset+start_ind:offset+end_ind))
          sum_sqr = sum_sqr + ((month_model-month_streamflow)**2)
          start_ind = end_ind+1
          end_ind = end_ind + 29
        enddo
        basin_rmse(ibasin+1) = sqrt(sum_sqr/real((nmonths)))

      case ("weekly","Weekly","WEEKLY")
        sum_sqr = 0.0
        nmonths = floor(((end_cal-start_cal)+1)/7.0) !7 days in a week
        start_ind = start_cal
        end_ind = start_cal + 6
        do itime = 1,nmonths
          month_model = sum(log_model(ibasin+1,start_ind:end_ind))
          month_streamflow = sum(log_streamflow(offset+start_ind:offset+end_ind))
          sum_sqr = sum_sqr + ((month_model-month_streamflow)**2)
          start_ind = end_ind+1
          end_ind = end_ind + 6
        enddo
        basin_rmse(ibasin+1) = sqrt(sum_sqr/real((nmonths)))

      case ("pentad","Pentad","PENTAD")
        sum_sqr = 0.0
        nmonths = floor(((end_cal-start_cal)+1)/5.0) !5 days in a pentad
        start_ind = start_cal
        end_ind = start_cal + 4
        do itime = 1,nmonths
          month_model = sum(log_model(ibasin+1,start_ind:end_ind))
          month_streamflow = sum(log_streamflow(offset+start_ind:offset+end_ind))
          sum_sqr = sum_sqr + ((month_model-month_streamflow)**2)
          start_ind = end_ind+1
          end_ind = end_ind + 4
        enddo
        basin_rmse(ibasin+1) = sqrt(sum_sqr/real((nmonths)))

      case ("custom")
        rmse_period = 5
        sum_sqr = 0.0
        nmonths = floor(((end_cal-start_cal)+1)/real(rmse_period,kind(sp))) !
        start_ind = start_cal
        end_ind = start_cal + (rmse_period -1)
        do itime = 1,nmonths
          month_model = sum(log_model(ibasin+1,start_ind:end_ind))
          month_streamflow = sum(log_streamflow(offset+start_ind:offset+end_ind))
          sum_sqr = sum_sqr + ((month_model-month_streamflow)**2)
          start_ind = end_ind+1
          end_ind = end_ind + (rmse_period -1)
        enddo
        basin_rmse(ibasin+1) = sqrt(sum_sqr/real((nmonths)))

      case DEFAULT
        sum_sqr = 0.0
        do itime = start_cal,end_cal-1
          sum_sqr = sum_sqr + ((log_model(ibasin+1,itime)-log_streamflow(itime+offset))**2.0)
          if(opt .ne. 1) then
            write(unit=88,*) model(ibasin+1,itime),streamflow(itime+offset)
          endif
        enddo
    end select
  enddo

  if(opt .ne. 1) then
    close(unit=88)
  endif

  !calculate rmse
  rmse = sum(basin_rmse*obj_fun_weight)

  return

end subroutine calc_rmse_region

!----------------------------------------------------
subroutine calc_nse_region(qsim,qobs,nse)

  implicit none

!input variables (model: simulations, qobs: observations)
  real(dp), dimension(:,:), intent(in)  :: qsim 
  real(dp), dimension(:),   intent(in)  :: qobs
!output variables
  real(dp),                 intent(out) :: nse 
!local variables
  integer(i4b)                           :: itime,ibasin,total_len,offset
  real(dp)                              :: sumSqrErr
  real(dp)                              :: sumSqrDev
  real(dp)                              :: sumQ
  real(dp)                              :: meanQ
  integer(i4b)                           :: nargs,nstream,nb,nday
  character(len=strLen)                 :: out_name
  character(len=strLen),dimension(10)   :: tokens
  character(len=strLen)                 :: last_token
  character(len=1)                      :: delims
  integer(i4b),allocatable,dimension(:)  :: basin_id
  real(dp),allocatable,dimension(:)     :: obj_fun_weight
  real(dp),allocatable,dimension(:)     :: basin_nse          ! nse for individual basin
  integer(i4b)                           :: nmonths            ! for monthly rmse calculation
  integer(i4b)                           :: start_ind, end_ind
  integer(i4b)                           :: start_obs, end_obs
  real(dp)                              :: month_qsim, month_qobs

  total_len = (end_cal-start_cal)*nbasin
  
  ! variable allocation
  allocate(obj_fun_weight(nbasin))
  allocate(basin_nse(nbasin))
  allocate(basin_id(nbasin))

  ! Output qobs to a different file if opt .ne. 1
  delims='/'
  if(opt .ne. 1) then
    call parse(outputobs_name,delims,tokens,nargs)
    last_token = tokens(nargs)
    out_name = trim(outputmod_name)//last_token(1:9)//"_flow.txt"
    print *,trim(out_name)
    open(unit=88,file=out_name)

    do ibasin = 0,nbasin-1
      offset = ibasin*sim_len
      do itime = start_cal,end_cal-1
        write(unit=88,*) qsim(ibasin+1,itime),qobs(itime+offset)
      enddo
    enddo
    close(unit=88)
  endif

  ! Read basin weight file
  ! this file determines how much each basin contributes to the total objective function 
  ! weights need to sum to 1 in the file
  open (UNIT=58,file=trim(basin_objfun_weight_file),form='formatted',status='old')
  do ibasin = 1,nbasin
    read (UNIT=58,*) basin_id(ibasin),obj_fun_weight(ibasin)
  enddo
  close(UNIT=58)

  !need to make sure i'm using the appropriate parts of the catenated region observed qobs timeseries
  do ibasin = 0,nbasin-1
    !offset places me at the start of each basin
    offset = ibasin*sim_len

    sumQ = 0.0
    sumSqrDev = 0.0
    sumSqrErr = 0.0

    !run from the starting calibration point to the ending calibration point
    select case (eval_length)
      case ("daily","Daily","DAILY")
        ! Compute Qob mean
        sumQ = sum(qobs(start_cal+offset:end_cal-1+offset))
        meanQ = sumQ/real((end_cal-start_cal))
        ! Compute sum of squre of error and deviation from menan (for qobs) 
        do itime = start_cal,end_cal-1
          sumSqrDev = sumSqrDev + (qobs(itime+offset)-meanQ)**2
          sumSqrErr = sumSqrErr + (qsim(ibasin+1,itime)-qobs(itime+offset))**2
        enddo
        ! Compute nse for current basin 
        basin_nse(ibasin+1) = sumSqrErr/sumSqrDev
     
      case ("monthly","Monthly","MONTHLY")
        nmonths = floor(((end_cal-start_cal)+1)/30.0)  !use 30 day months uniformly to make it easier
        ! Compute montly observed Q
        ! Indices of start and end for first month
        start_ind = start_cal
        end_ind = start_cal + 29
        do itime = 1,nmonths
          month_qobs = sum(qobs(offset+start_ind:offset+end_ind))
          sumQ       = sumQ+month_qobs
        enddo
        meanQ = sumQ/real(nmonths)
        ! Compute sum of squre of error and deviation from menan (for qobs) 
        start_ind = start_cal
        end_ind = start_cal + 29
        do itime = 1,nmonths
          month_qsim = sum(qsim(ibasin+1,start_ind:end_ind))
          month_qobs = sum(qobs(offset+start_ind:offset+end_ind))
          sumSqrErr = sumSqrErr + ((month_qsim-month_qobs)**2)
          sumSqrDev = sumSqrDev + ((month_qobs-meanQ)**2)
          !update starting and ending indice for next month step
          start_ind = end_ind+1
          end_ind = end_ind + 29
        enddo
        !grab remainder portion and weight it by number of days
        month_qsim = sum(qsim(ibasin+1,end_ind+1:end_cal-1))
        month_qobs = sum(qobs(offset+end_ind+1:end_cal-1))
        sumSqrErr = sumSqrErr + ((month_qsim-month_qobs)**2) * real((end_cal-1-(offset+end_ind+1))/30.0)
        sumSqrDev = sumSqrDev + ((month_qobs-meanQ)**2) * real((end_cal-1-(offset+end_ind+1))/30.0)
        ! Compute nse for current basin 
        basin_nse(ibasin+1) = sumSqrErr/sumSqrDev

      case default
        ! Compute Qob mean
        sumQ = sum(qobs(start_cal+offset:end_cal-1+offset))
        meanQ = sumQ/real((end_cal-start_cal)+1)
        ! Compute sum of squre of error and deviation from menan (for qobs) 
        do itime = start_cal,end_cal-1
          sumSqrDev = sumSqrDev + (qobs(itime+offset)-meanQ)**2
          sumSqrErr = sumSqrErr + (qsim(ibasin+1,itime)-qobs(itime+offset))**2
        enddo
        ! Compute nse for current basin 
        basin_nse(ibasin+1) = sumSqrErr/sumSqrDev
    end select

  enddo

  !calculate rmse
  nse = sum(basin_nse*obj_fun_weight)

  return

end subroutine calc_nse_region

!***************************
subroutine calc_kge(model,streamflow,kge)
!***************************
! This subroutine calculates the negative of the Kling-Gupta Efficiency
  implicit none
!input variables (model: simulations, streamflow: observations)
  real(dp), dimension(:),     intent(in)  :: model
  real(dp), dimension(:),     intent(in)  :: streamflow
!output variables
  real(dp),                   intent(out) :: kge
!local variables
  integer(i4b)                             :: itime
  real(dp)                                :: cc,alpha,betha,mu_s,mu_o,sigma_s,sigma_o
  
  mu_s = 0.0
  mu_o = 0.0
  sigma_s = 0.0
  sigma_o = 0.0  

  ! first compute the mean
  mu_s = sum(model(start_cal:end_cal))/real(end_cal-start_cal+1)
  mu_o = sum(streamflow(start_cal:end_cal))/real(end_cal-start_cal+1)
  betha = mu_s/mu_o

  ! Now compute the standard deviation
  do itime = start_cal,end_cal-1
    sigma_s = sigma_s + (model(itime)-mu_s)**2
    sigma_o = sigma_o + (streamflow(itime)-mu_o)**2
  enddo
  sigma_s = sqrt(mu_s/real(end_cal-start_cal+1))
  sigma_o = sqrt(mu_s/real(end_cal-start_cal+1))
  alpha = sigma_s/sigma_o
  
  !Compute linear correlation coefficient
  call pearsn(model(start_cal:end_cal-1),streamflow(start_cal:end_cal-1),cc)

  kge = -( 1.0 - sqrt((cc-1.0)**2 + (alpha-1.0)**2 + (betha-1.0)**2) )

  return
end subroutine

!***********************************************************************
! This subroutine calculates the negative of the Kling-Gupta Efficiency
!***********************************************************************
subroutine calc_kge_region(model,streamflow,kge)
  implicit none

!  use strings
!input variables (model: simulations, streamflow: observations)
  real(dp), dimension(:,:), intent(in)  :: model
  real(dp), dimension(:),   intent(in)  :: streamflow
!output variables
  real(dp),                 intent(out) :: kge
!local variables
  integer                               :: itime
  real(dp)                              :: cc,alpha,betha,mu_s,mu_o,sigma_s,sigma_o
  integer                               :: ibasin,total_len,offset,cnt
  double precision                      :: sum_sqr
  integer                               :: nargs,nstream,nb,nday
  character(len=2000)                   :: out_name
  character(len=1000),dimension(10)     :: tokens
  character(len=1000)                   :: last_token
  character(len=1)                      :: delims
  real(dp),dimension(:),allocatable     :: model_local
  real(dp),dimension(:),allocatable     :: obs_local

  total_len = (end_cal-start_cal)*nbasin

  allocate(model_local(total_len))
  allocate(obs_local(total_len))

!want to output streamflow to a different file of opt .ne. 1
  delims='/'

  if(opt .ne. 1) then
    call parse(outputobs_name,delims,tokens,nargs)
    last_token = tokens(nargs)

    out_name = trim(outputmod_name)//last_token(1:9)//"_flow.txt"
    print *,trim(out_name)
    open(unit=88,file=out_name)
  endif

  cnt = 1
  do ibasin = 0,nbasin-1
    !offset places me at the start of each basin
    offset = ibasin*sim_len
    do itime = start_cal,end_cal-1
      model_local(cnt) = model(ibasin+1,itime)
      obs_local(cnt) = streamflow(offset+itime)
      cnt = cnt + 1
    enddo
  enddo

  !!set variables to zero
  mu_s = 0.0
  mu_o = 0.0
  sigma_s = 0.0
  sigma_o = 0.0 

  !offset places me at the start of each basin
  offset = ibasin*sim_len
  mu_s = sum(model_local)/real(total_len)
  mu_o = sum(obs_local)/real(total_len)
  betha = mu_s/mu_o

  !Now we compute the standard deviation
  do itime = 1,total_len
    sigma_s = sigma_s + (model_local(itime)-mu_s)**2
    sigma_o = sigma_o + (obs_local(itime)-mu_o)**2
    if(opt .ne. 1) then
      write(unit=88,*) model_local(itime),obs_local(itime)
    endif
  enddo   !end itime loop

  sigma_s = sqrt(mu_s/real(total_len))
  sigma_o = sqrt(mu_s/real(total_len))
  alpha = sigma_s/sigma_o
    
  !Compute linear correlation coefficient
  call pearsn(model_local,obs_local,cc)

  kge = -( 1.0 - sqrt((cc-1.0)**2 + (alpha-1.0)**2 + (betha-1.0)**2) )

  if(opt .ne. 1) then
    close(unit=88)
  endif
  
  return
end subroutine calc_kge_region

! -----------------------------
subroutine pearsn(x,y,r)

  implicit none

!input variables
  real(dp), dimension(:), intent(in)  :: x
  real(dp), dimension(:), intent(in)  :: y
!output variables
  real(dp),               intent(out) :: r
!local variables
  real(dp)                            :: tiny = 1.0e-20
  real(dp), dimension(size(x))        :: xt,yt
  real(dp)                            :: ax,ay,sxx,sxy,syy
  integer(i4b)                        :: n

  n=size(x)
  ax=sum(x)/n
  ay=sum(y)/n
  xt(:)=x(:)-ax
  yt(:)=y(:)-ay
  sxx=dot_product(xt,xt)
  syy=dot_product(yt,yt)
  sxy=dot_product(xt,yt)
  r=sxy/(sqrt(sxx*syy)+tiny)

end subroutine

!******************************
subroutine read_obsflow(fname, nTime, obs)
!******************************
  implicit none
! This subroutine reads observed streamflow

!input variables
  character(len=strLen), intent(in)    :: fname
  integer(i4b),intent(in)              :: nTime 
!output variables
  real(dp), dimension(:), intent(out)  :: obs
!local variables
  integer(i4b)                         :: itime

!read observed streamflow
  open (UNIT=52,file=trim(fname),form='formatted',status='old')
  do itime = 1,nTime
    read (UNIT=52,*) obs(itime)
  enddo
  close(UNIT=52)

  return
end subroutine

!******************************
subroutine read_obsflow_region(fname,nTime,nbsn,obs)
!******************************
  implicit none
! This subroutine reads observed streamflow

!input variables
  character(len=strLen), intent(in)    :: fname
  integer(i4b),          intent(in)    :: nTime
  integer(i4b),          intent(in)    :: nbsn
!output variables
  real(dp), dimension(:),intent(out)   :: obs
!local variables
  integer                              :: itime   ! loop index

!read observed streamflow
  open (UNIT=52,file=trim(fname),form='formatted',status='old')
  do itime = 1,nTime*nbsn
    read (UNIT=52,*) obs(itime)
  enddo
  close(UNIT=52)

  return
end subroutine

!******************************
subroutine read_vic_cells(outputdir,cellfile,listfile,infofile,nTime,vic_runoff)
!******************************
  implicit none

!input variables
  character(len=strLen),intent(in)    :: outputdir
  character(len=strLen),intent(in)    :: cellfile
  character(len=strLen),intent(in)    :: listfile
  character(len=strLen),intent(in)    :: infofile
  integer(i4b),          intent(in)    :: nTime
!output variables
  real(dp),dimension(:,:),intent(out) :: vic_runoff
!local variables
  character(len = 200)                :: filename
  real(dp)                            :: cellfraction,basin_area
  real(dp)                            :: auxflux(5) ! This is only in case of water balance mode
  integer                             :: ibasin, itime, ivar,icell,ncell
  integer                             :: dum,j,c_cell

!set output variable to zero
  vic_runoff = 0.0
!cell counter
  c_cell = 1
!open a few files
  open (UNIT=53,file=listfile,form='formatted',status='old')
  open (UNIT=54,file=cellfile,form='formatted',status='old')
  open (UNIT=51,file=infofile,form='formatted',status='old')

  do ibasin = 1,nbasin
    read (UNIT=51,*) dum,dum,basin_area,ncell
    do icell = 1,ncell
      read (UNIT=53,*) filename
      read (UNIT=54,*) cellfraction
      filename = TRIM(outputdir) // TRIM(filename)
      open (UNIT=55,file= filename,form='formatted',status='old')
      do itime = 1,nTime
        read (UNIT=55,*) (auxflux(ivar), ivar=1,5)
        vic_runoff(c_cell,itime) = (auxflux(4) + auxflux(5))*cellfraction
      enddo  !end time loop
      close(UNIT=55)
      c_cell = c_cell + 1
    enddo  !end cell loop
  enddo  !end basin loop

  close(UNIT=51)
  close(UNIT=53)
  close(UNIT=54)

  return
end subroutine

!******************************
subroutine route_grid_cell(qmod_in,route_mod,ushape,uscale)
!******************************
  implicit none

  !input variables
  real(dp), dimension(:,:), intent(in)  :: qmod_in
  real(dp), intent(in)                  :: ushape,uscale
  !output variables
  real(dp), dimension(:,:), intent(out) :: route_mod
  !local variables
  integer(i4b)                           :: i,k,m,ntau,icell,j,ncell
  real(dp)                              :: dtuh

!number of grid cells coming in
  ncell = size(qmod_in,1)
! route flow for each grid cell
  dtuh = real(86400./86400.) 
  ntau = 0
  if (ushape .le. 0.0 .and. uscale .le. 0.0) THEN
    do icell=1,ncell
      route_mod(icell,:) = qmod_in(icell,:)
    enddo
  else
    do icell=1,ncell
      call duamel(qmod_in(icell,1:sim_len-1),ushape,uscale,dtuh,sim_len-1,route_mod(icell,1:sim_len-1),ntau)
    enddo
  end if

  return
end subroutine

!******************************
subroutine agg_cells_to_basin(infofile,qmod_in,vic_runoff)
!******************************
  implicit none

!input variables
  character(len=strLen),intent(in)    :: infofile
  real(dp),dimension(:,:),intent(in)  :: qmod_in
!output variables
  real(dp),dimension(:,:),intent(out) :: vic_runoff
!local variables
  real(dp)                            :: cellfraction,basin_area
  real(dp)                            :: auxflux(5) ! This is only in case of water balance mode
  integer                             :: ibasin, itime, ivar,icell,ncell
  integer                             :: dum,j,c_cell

!set output variable to zero
  vic_runoff = 0.0
!cell counter
  c_cell = 1
!open a few files
  open (UNIT=51,file=infofile,form='formatted',status='old')

  do ibasin = 1,nbasin
    read (UNIT=51,*) dum,dum,basin_area,ncell
    do icell = 1,ncell
      vic_runoff(ibasin,:) = vic_runoff(ibasin,:) + qmod_in(c_cell,:)
      c_cell = c_cell + 1
    enddo  !end cell loop
  enddo  !end basin loop

  close(UNIT=51)
  close(UNIT=53)
  close(UNIT=54)

  return
end subroutine

!******************************
subroutine read_modflow(outputdir,cellfile,listfile,nTime,model)
!******************************
  implicit none
!input variables
  character(len=strLen), intent(in) :: outputdir
  character(len=strLen), intent(in) :: cellfile
  character(len=strLen), intent(in) :: listfile
  integer(i4b),           intent(in) :: nTime
!output variables
  real(dp),dimension(:),intent(out) :: model
!local variables
  character(len=strLen)             :: filename
  real(dp)                          :: cellfraction
  real(dp)                          :: auxflux(5) ! This is only in case of water balance mode
  integer(i4b)                       :: icell, itime, ivar

!set model to zero
  model = 0.0

!read model flow
  open (UNIT=53,file=trim(listfile),form='formatted',status='old')
  open (UNIT=54,file=trim(cellfile),form='formatted',status='old')

  do icell = 1,Ncells
    read (UNIT=53,*) filename
    read (UNIT=54,*) cellfraction
    filename = trim(outputdir) // trim(filename)
    open (UNIT=55,file= filename,form='formatted',status='old')
    do itime = 1,nTime
      read (UNIT=55,*) (auxflux(ivar), ivar=1,5)
      model(itime) = model(itime) + (auxflux(4) + auxflux(5))*cellfraction
    enddo
    close(UNIT=55)
  enddo

  close(UNIT=53)
  close(UNIT=54)

  return
end subroutine

!******************************
subroutine read_modflow_region(outputdir,cellfile,listfile,infofile,nTime ,model_region)
!******************************
  implicit none
!input variables
  character(len=strLen), intent(in)   :: outputdir
  character(len=strLen), intent(in)   :: cellfile
  character(len=strLen), intent(in)   :: listfile
  character(len=strLen), intent(in)   :: infofile
  integer(i4b),           intent(in)   :: nTime
!output variables
  real(dp),dimension(:,:),intent(out) :: model_region
!local variables
  character(len = 200)                :: filename
  real(dp)                            :: cellfraction,basin_area
  real(dp)                            :: auxflux(5) ! This is only in case of water balance mode
  integer                             :: ibasin, itime, ivar,icell,ncell
  integer                             :: dum,j

!set output variable to zero
  model_region = 0.0

!open a few files
  open (UNIT=53,file=trim(listfile),form='formatted',status='old')
  open (UNIT=54,file=trim(cellfile),form='formatted',status='old')
  open (UNIT=51,file=trim(infofile),form='formatted',status='old')

  do ibasin = 1,nbasin
    read (UNIT=51,*) dum,dum,basin_area,ncell

    do icell = 1,ncell
      read (UNIT=53,*) filename
      read (UNIT=54,*) cellfraction

      filename = TRIM(outputdir) // TRIM(filename)
      open (UNIT=55,file= filename,form='formatted',status='old')
      
      do itime = 1,nTime
        read (UNIT=55,*) (auxflux(ivar), ivar=1,5)
	! The next line is the sum of baseflow + surface runoff
	!qmod(itime) = qmod(itime) + (auxflux(4) + auxflux(5))*cellfraction*BasinArea*( 10 ** 3 )/86400.0 !Metric unit: cms
  !      qmod(itime) = qmod(itime) + 0.0033*(auxflux(4) + auxflux(5))*cellfraction*BasinArea*10.764*( 10 ** 6 )/86400.0 ! English Unit: cfs, basin area in km^2
  !      qmod(icell,itime) = 0.0033*(auxflux(4) + auxflux(5))*cellfraction*BasinArea*10.764/86400.0 ! English Unit: cfs,   basin area in m^2
!	model_region(ibasin,itime) = model_region(ibasin,itime) + 0.0033*(auxflux(4) + auxflux(5))*cellfraction*basin_area*10.764/86400.0 ! English Unit: cfs,   basin area in m^2
        model_region(ibasin,itime) = model_region(ibasin,itime) + (auxflux(4) + auxflux(5))*cellfraction
      enddo
      close(UNIT=55)
    enddo
  enddo
  close(UNIT=51)
  close(UNIT=53)
  close(UNIT=54)

  return
end subroutine

!******************************
subroutine route_modflow_region(qmod_in,route_mod,ushape,uscale)
!******************************

  implicit none

  !input variables
  real(dp),dimension(:,:), intent(in)  :: qmod_in
  real(dp),intent(in)                  :: ushape,uscale
  !output variables
  real(dp),dimension(:,:), intent(out) :: route_mod
  !local variables
  integer(i4b)                          :: i,k,m,ntau,ibasin,j
  real(dp)                             :: dtuh

! route flow for each basin in the region now
  dtuh = real(86400./86400.) 
  ntau = 0
  if (ushape .le. 0.0 .and. uscale .le. 0.0) THEN
    do ibasin=1,nbasin
      route_mod(ibasin,:) = qmod_in(ibasin,:)
    enddo
  else
    do ibasin=1,nbasin
      call duamel(qmod_in(ibasin,1:sim_len-1),ushape,uscale,dtuh,sim_len-1,route_mod(ibasin,1:sim_len-1),ntau)
    enddo
  end if
end subroutine


end module vic_subroutines_parallel
