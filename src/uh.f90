module uh
  
  use nrtype 
  use public_var

  implicit none

  public :: duamel 
  private

contains

  !----------------------------------
  ! perform unit hydrograph convolution 
  subroutine duamel(Q,un1,ut,dt,nq,QB,ntau,inUH)
    implicit none

    ! input 
    real(dp),   dimension(:),          intent(in)  :: Q      ! instantaneous flow
    real(dp),                          intent(in)  :: un1    ! scale parameter
    real(dp),                          intent(in)  :: ut     ! time parameter
    real(dp),                          intent(in)  :: dt     ! time step 
    integer(i4b),                      intent(in)  :: nq     ! size of instantaneous flow series
    integer(i4b),                      intent(in)  :: ntau 
    real(dp),   dimension(:),optional, intent(in)  :: inUH   ! optional input unit hydrograph  
    ! output
    real(dp),dimension(:),             intent(out) :: QB     ! Routed flow
    ! local 
    real(dp),dimension(1000)                       :: uh     ! unit hydrograph (use 1000 time step)
    integer(i4b)                                   :: m      ! size of unit hydrograph
    integer(i4b)                                   :: A,B
    integer(i4b)                                   :: i,j,ij ! loop index 
    integer(i4b)                                   :: ioc    ! total number of time step  
    real(dp)                                       :: top
    real(dp)                                       :: toc
    real(dp)                                       :: tor
    real(dp)                                       :: spv    ! cumulative uh distribution to normalize it to get unit hydrograph
    
    !size of unit hydrograph 
    m=size(uh)
    !initialize unit hydrograph 
    uh=0._dp
    ! Generate unit hydrograph
    if (un1 .lt. 0) then ! if un1 < 0, routed flow = instantaneous flow 
      uh(1)=1.0_dp
      m = 1
    else
      if (present(inUH)) then  !update uh and size of uh
        uh=inUH
        m=size(uh)  
      else 
        spv=0.0_dp
        toc=gf(un1)
        toc=log(toc*ut)
        do i=1,m
          top=i*DT/ut
          tor=(UN1-1)*log(top)-top-toc
          uh(i)=0.0_dp
          if(tor.GT.-8.0_dp) then 
            uh(i)=exp(tor)
          else 
            if (i .GT. 1) then 
              uh(i) = 0.0_dp
            end if 
          end if 
          spv=spv+uh(i) ! accumulate uh each uh time step
        end do
        if (spv .EQ. 0) spv=1.0E-5
        spv=1.0_dp/spv  
        do i=1,m
          uh(I)=uh(i)*spv  ! normalize uh so cumulative uh = 1
        end do
      endif
    endif
      
    ! do unit hydrograph convolution
    IOC=nq+NTAU
    if (nq.LE.m) then
      do i=1,IOC
        QB(i)=0.0_dp
        A=1
        if(i.GT.m) A=I-m+1
        B=I
        if(i.GT.nq) B=nq
        do j=A,B
          ij=i-j+1
          QB(i)=QB(i)+Q(J)*uh(ij)
        end do
      end do
    else
      do i=1,IOC
        QB(i)=0.0_dp
        A=1
        if(i.GT.nq) A=i-nq+1
        B=i
        if(i.GT.M) B=M 
        do j=A,B
          ij=i-j+1
          QB(i)=QB(i)+uh(J)*Q(ij)
        end do
      end do 
    end if
  
  end subroutine duamel
  
  !=================================================================
  function gf(Y)
  
    implicit none
  
    real(dp),intent(in)  :: y
    real(dp)             :: gf 
    real(dp)             :: x
    real(dp)             :: h
  
    H=1_dp
    x=y

    do 
      if(x.le.0_dp) exit
      if(x.lt.2_dp .and. x.gt.2_dp) then
        gf=H
        exit
      end if
      if(x.gt.2_dp) then
        if(x.le.3_dp) then
          x=x-2_dp
          h=(((((((.0016063118_dp*x+0.0051589951_dp)*x+0.0044511400_dp)*x+.0721101567_dp)*x  &
            +.0821117404_dp)*x+.4117741955_dp)*x+.4227874605_dp)*x+.9999999758_dp)*h
          gf=H 
          exit
        else
          x=x-1_dp
          H=H*x
          cycle
        end if
      else
        H=H/x
        x=x+1_dp
        cycle
      end if
    end do
  
  end function gf

end module uh
