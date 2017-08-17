module hydroSignature 

  use nrtype 
  use public_var 

  implicit none

  private

  public :: cal_rr
  public :: cal_eqp
  public :: cal_mean_yr
  public :: cal_fms
  public :: cal_qp
  public :: cal_bfi
  public :: cal_events

contains
  
  ! Hydrologic Signature: long-term runoff ratio 
  subroutine cal_rr(daily_q, daily_p, rr, err, message)
    implicit none
    !input variables
    real(dp),               intent(in)  :: daily_q(:)
    real(dp),               intent(in)  :: daily_p(:)
    !output variables
    real(dp),               intent(out) :: rr 
    integer(i4b),           intent(out) :: err                  ! error code
    character(len=strLen),  intent(out) :: message              ! error message
    !local variables
    real(dp)                            :: yr_q
    real(dp)                            :: yr_p
    character(len=strLen)               :: cmessage             ! error message from subroutine
  
    ! initialize error control
    err=0; message='calc_rr/'
    call cal_mean_yr(daily_q, yr_q, err, cmessage)
    call cal_mean_yr(daily_p, yr_p, err, cmessage)
    rr=yr_q/yr_p
    return
  end subroutine

  ! Hydrologic Signature: Q-P elastisity 
  subroutine cal_eqp(daily_q, daily_p, eqp, err, message)
    implicit none
    !input variables
    real(dp),               intent(in)  :: daily_q(:)
    real(dp),               intent(in)  :: daily_p(:)
    !output variables
    real(dp),               intent(out) :: eqp 
    integer(i4b),           intent(out) :: err                  ! error code
    character(len=strLen),  intent(out) :: message              ! error message
    !local variables
    real(dp)                            :: rr 
    real(dp), allocatable               :: yr_q(:)
    real(dp), allocatable               :: yr_p(:)
    real(dp), allocatable               :: diff_yr_q(:)
    real(dp), allocatable               :: diff_yr_p(:)
    real(dp), allocatable               :: yr_eqp(:)
    integer(i4b)                        :: iTime
    integer(i4b)                        :: nTime,nYr            ! number of years and annual time steps 
    integer(i4b)                        :: start_ind, end_ind
    character(len=strLen)               :: cmessage             ! error message from subroutine
  
    ! initialize error control
    err=0; message='calc_eqp/'
    nTime=size(daily_q)
    if (size(daily_p)/=nTime) then; err=10; message=trim(message)//'q and p sizes different'; return; endif
    nYr = floor(nTime/365.0_dp)  !use 365 day per year uniformly to make it easier
    allocate(yr_q(nYr))
    allocate(yr_p(nYr))
    allocate(diff_yr_q(nYr-1))
    allocate(diff_yr_p(nYr-1))
    allocate(yr_eqp(nYr-1))
    call cal_rr(daily_q, daily_p, rr, err, cmessage)
    ! Indices of start and end for first year 
    start_ind = 1 
    end_ind   = start_ind + 364 
    do iTime = 1,nYr
      yr_q(iTime) = sum(daily_q(start_ind:end_ind))
      yr_q(iTime) = yr_q(iTime)/365.0_dp

      yr_p(iTime) = sum(daily_p(start_ind:end_ind))
      yr_p(iTime) = yr_p(iTime)/365.0_dp

      if (iTime >= 2) then
        diff_yr_p(iTime-1) = yr_p(iTime)-yr_p(iTime-1)
        diff_yr_q(iTime-1) = yr_q(iTime)-yr_q(iTime-1)
        yr_eqp(iTime-1)=diff_yr_q(iTime-1)/diff_yr_p(iTime-1)/rr
      endif
      !update starting and ending indice for next year step
      start_ind = end_ind+1
      end_ind   = end_ind+364
    enddo
    eqp=median(yr_eqp)
    return
  end subroutine

  ! Hydrologic Signature: long-term annual mean value 
  subroutine cal_mean_yr(daily_q, mean_yr, err, message)
    implicit none
    !input variables
    real(dp),               intent(in)  :: daily_q(:)           ! daily discharge [mm/day]
    !output variables
    real(dp),               intent(out) :: mean_yr 
    integer(i4b),           intent(out) :: err                  ! error code
    character(len=strLen),  intent(out) :: message              ! error message
    !local variables
    real(dp), allocatable               :: yearly(:) 
    integer(i4b)                        :: iTime
    integer(i4b)                        :: nTime,nYr            ! number of years and annual time steps 
    integer(i4b)                        :: start_ind, end_ind
  
    ! initialize error control
    err=0; message='mean_year/'
    nTime=size(daily_q)
    nYr = floor(nTime/365.0_dp)  !use 365 day per year uniformly to make it easier
    allocate(yearly(nYr))
    ! Indices of start and end for first year 
    start_ind = 1 
    end_ind   = start_ind + 364 
    do itime = 1,nYr
      yearly(itime) = sum(daily_q(start_ind:end_ind))
      yearly(iTime) = yearly(iTime)/365.0_dp
      !update starting and ending indice for next year step
      start_ind = end_ind+1
      end_ind   = end_ind+364
    enddo
    mean_yr=sum(yearly)/real(nYr)
    return
  end subroutine 
  
  ! Hydrologic Signature: Slope of mid section of flow dulation curve(FMS)
  subroutine cal_fms(daily_q, FMS, err, message)
    ! Sawicz et al., 2011. Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA 
    ! HESS 15, 2895-2911, doi:10.5194/hess-15-2895-2011 
    implicit none
    !input variables 
    real(dp),              intent(in)  :: daily_q(:)           ! daily discharge [mm/day]
    !output variables
    real(dp),             intent(out)  :: FMS
    integer(i4b),         intent(out)  :: err           ! error code
    character(len=strLen),intent(out)  :: message       ! error message
    ! local variables
    integer(i4b)                       :: iTime ! loop index
    integer(i4b)                       :: i33,i66
    integer(i4b)                       :: idx(1)
    integer(i4b)                       :: nTime            
    real(dp),dimension(size(daily_q))  :: p               ! probability
    real(dp),dimension(size(daily_q))  :: dailyIn, ascDaily
  
    ! initialize error control
    err=0; message='FMS/'
    nTime=size(daily_q)
    dailyIn=daily_q
    where(dailyIn<verySmall) dailyIn=verySmall
    p=(/ (real(itime)/(real(nTime)+1.0_dp), iTime=1,nTime) /)
    idx=minloc(abs(p-0.66_dp))
    i66=idx(1)
    idx=minloc(abs(p-0.33_dp))
    i33=idx(1)
    ascDaily=dailyIn
    call sort(ascDaily)
    FMS=( log(ascDaily(i66))-log(ascDaily(i33)) )/(0.66_dp - 0.33_dp)
    return
  end subroutine
  
  ! Hydrologic signature: "perc" percentile value of time series (Q"perc")
  subroutine cal_qp(daily_q, QP, perc, err, message) 
    implicit none
    !input variables 
    real(dp),             intent(in)   :: daily_q(:)    ! daily discharge [mm/day]
    real(dp),             intent(in)   :: perc          ! percentile (0-1) vector is ascending order
    !output variables
    real(dp),             intent(out)  :: QP 
    integer(i4b),         intent(out)  :: err           ! error code
    character(len=strLen),intent(out)  :: message       ! error message
    ! local variables
    integer(i4b)                       :: iTime ! loop index
    integer(i4b)                       :: ip
    integer(i4b)                       :: idx(1)
    integer(i4b)                       :: nTime            
    real(dp),dimension(size(daily_q))  :: p               ! probability
    real(dp),dimension(size(daily_q))  :: dailyIn, ascDaily
  
    err=0; message='qp/'
    if (perc < 0 .or. perc > 1)then; err=10; message=trim(message)//'percentile is outside range'; return; endif
    nTime=size(daily_q)
    dailyIn=daily_q
    where(dailyIn<verySmall) dailyIn=verySmall
    p=(/ (real(itime)/(real(nTime)+1.0_dp), iTime=1,nTime) /)
    idx=minloc(abs(p-perc))
    ip=idx(1)
    ascDaily=dailyIn
    call sort(ascDaily)
    QP=ascDaily(ip)
    return
  end subroutine
  
  ! Hydrologic Signature: Baseflow Index (BFI) 
  subroutine cal_bfi(daily_q, BFI, err, message, coef)
    ! Sawicz et al., 2011. Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA 
    ! HESS 15, 2895-2911, doi:10.5194/hess-15-2895-2011 
    implicit none
    !input variables 
    real(dp),             intent(in)   :: daily_q(:)    ! daily discharge [mm/day]
    real(dp), optional,   intent(in)   :: coef          ! digital filter parameter
    !output variables
    real(dp),             intent(out)  :: BFI 
    integer(i4b),         intent(out)  :: err           ! error code
    character(len=strLen),intent(out)  :: message       ! error message
    ! local variables
    integer(i4b)                       :: iTime         ! loop index
    integer(i4b)                       :: nTime            
    real(dp)                           :: co 
    real(dp),dimension(size(daily_q))  :: bDaily        ! baseflow 
  
    err=0; message='qp'
    co=0.925_dp
    if (present(coef))  co=coef
    nTime=size(daily_q)
    bDaily(1)=daily_q(1) 
    do iTime=2,nTime
      bDaily(iTime)=co*bDaily(iTime-1)+((1.0_dp-co)/2.0_dp)*(daily_q(iTime)+daily_q(iTime-1))
      if (bDaily(iTime) > daily_q(iTime)) bDaily(iTime)=daily_q(iTime)
    end do
    BFI=sum(bDaily)/sum(daily_q)
    return
  end subroutine
  
  ! Hydrologic Signature: Average duration (days) and numbers per yr (-) of flow events 
  subroutine cal_events(daily_q, thresh, FRE, DUR, err, message)
    implicit none
    real(dp),               intent(in)  :: daily_q(:)           ! daily discharge [mm/day]
    real(dp),               intent(in)  :: thresh
    !output variables
    real(dp),               intent(out) :: FRE                  ! mean number of events per year [yr-1]
    real(dp),               intent(out) :: DUR                  ! mean event duration [days] 
    integer(i4b),           intent(out) :: err                  ! error code
    character(len=strLen),  intent(out) :: message              ! error message
    !local variables
    integer(i4b)                        :: iTime, jTime
    integer(i4b)                        :: c 
    integer(i4b)                        :: duration 
    integer(i4b)                        :: event 
    real(dp)                            :: med_daily 
    real(dp),    allocatable            :: daily_yr(:) 
    real(dp),    allocatable            :: meanDur(:) 
    integer(i4b),allocatable            :: totEvent(:) 
    integer(i4b)                        :: nTime,nYr            ! number of years and annual time steps 
    integer(i4b)                        :: start_ind, end_ind
    character(len=strLen)               :: cmessage             ! error message from subroutine
  
    ! initialize error control
    err=0; message='highEvents/'

    if (thresh < 0 .or. thresh > 10)then; err=10; message=trim(message)//'multiplier "thresh" is outside range'; return; endif
    nTime=size(daily_q)
    nYr = floor(nTime/365.0_dp)  !use 365 day per year uniformly to make it easier
    ! Indices of start and end for first year 
    start_ind = 1 
    end_ind   = start_ind + 364 
    call cal_qp(daily_q, med_daily, 0.5_dp, err, cmessage)
    allocate(totEvent(nYr))
    allocate(meanDur(nYr))
    do iTime = 1,nYr
      if (allocated(daily_yr)) deallocate(daily_yr)
      allocate(daily_yr(end_ind-start_ind+1))
      daily_yr=daily_q(start_ind:end_ind)
      c=0    
      event=0_i4b
      duration=0_i4b
      do jTime=1,size(daily_yr)
        if ( daily_yr(jTime)>med_daily*thresh ) then
          c=c+1
          if (c==1) event=event+1
        else
          duration=duration+c
          c=0
        endif
      end do
      if (event==0_i4b) then
        meanDur(iTime) = 0.0_dp 
      else
        meanDur(iTime)  = real(duration)/real(event)
      end if
      totEvent(iTime) = event
      !update starting and ending indice for next year step
      start_ind = end_ind+1
      end_ind   = end_ind+364
    enddo
    DUR=sum(meanDur)/real(nYr)
    FRE=real(sum(totEvent))/real(nYr)
    return
  end subroutine 

  function median(X)
     implicit none 
     real(dp),        intent(in) :: X(:)
     real(dp)                     :: median
     real(dp), allocatable        :: temp(:)
     integer(i4b)                 :: N

    N=size(X)
    allocate(temp, source=X)
    call sort(temp)               ! sort the copy
    if (mod(N,2) == 0) then       ! compute the median
      median = (temp(N/2) + temp(N/2+1)) / 2.0_dp
    else
      median = temp(N/2+1)
    end if
    return
  end function 
 
  subroutine pearsn(x,y,r)
    implicit none
    !input variables
    real(dp), dimension(:), intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    !output variables
    real(dp),               intent(out) :: r
    !local variables
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
    r=sxy/(sqrt(sxx*syy)+verySmall)
    return
  end subroutine
  
  subroutine sort(vec)
    implicit none
    real(dp) ,intent(inout) :: vec(:)
    real(dp)                :: buf
    integer(i4b)            :: nsize, i
    integer(i4b)            :: k(1)
  
    nsize = size(vec)
    do i = 1, nsize
      k   = minloc(vec(i:nsize))+i-1
      buf = vec(i)
      vec(i) = vec(k(1))
      vec(k(1)) = buf
    enddo
    return
  end subroutine

end module hydroSignature 
