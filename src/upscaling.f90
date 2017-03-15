module upscaling
! This module contain (up)scale operators from native geophysical data (geo-poly) resolution to model resolution (model-poly)
! For each model-poly, scaling operator takes areal weights (weight vector) and data values (data vector) of each geo-poly 
! and compute a single value 
!
! Scale operators and subroutines included here are 
! 1. powermean: weighted power mean
! 2. wamean: weighted arithmatic mean
! 3. wgmean: weighted geometric mean
! 4. whmean: weighted harmonic mean
! 5. wmedi:  weighted median 
! Also subroutine that allows to switch scaling method by specifying the method in input 

use nrtype                            ! variable types, etc.
use data_type                         ! Including custum data structure definition
use public_var                        ! Including common constant (physical constant, other e.g., dmiss, etc.)

implicit none

private

public::aggreg

contains

! *********************************************************************
! subroutine: Spatical aggregation with selected method 
! *********************************************************************
subroutine aggreg(wgtval, wgtvec, datvec, method, pnorm, err, message )
  implicit none
  ! input
  real(dp),              intent(in)      :: wgtvec(:)        ! weight vector
  real(dp),              intent(in)      :: datvec(:)        ! data value vector
  character(*),          intent(in)      :: method
  real(dp),              intent(in)      :: pnorm            ! p value for pnorm 
  ! output 
  real(dp),              intent(out)     :: wgtval           ! weighted value
  integer(i4b),          intent(out)     :: err              ! error code
  character(len=strLen), intent(out)     :: message          ! error message for current routine
  ! local 
  character(len=strLen)                  :: cmessage         ! error message from subroutine

  err=0; message="aggreg/"
  select case(trim(method))
    case('pnorm'); call wpnorm(wgtval, wgtvec, datvec, err, cmessage, p_exp_in=pnorm )
    case('wmedi'); print*,'No pnorm value required';call wmedi (wgtval, wgtvec, datvec, err, cmessage)
    case('wmode'); print*,'No pnorm value required';call wmode (wgtval, wgtvec, datvec, err, cmessage)
    case('asum');  print*,'No pnorm value required';call asum  (wgtval, wgtvec, datvec, err, cmessage)
    case default;  err=10; message=trim(message)//'aggre method not avaiable'; return
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage);return;endif
  return
end subroutine

! *********************************************************************
! Subroutine: Weighted p-norm aggregation 
!
! purpose: aggregating values according to a p-norm
!
!          result = sum(w * abs(x)^p)^(1/p),
!
!          where,
!          w = weights
!          x = values to aggregate
!          p = real number
!
! author: Stephan Thober (stephan.thober@ufz.de)
! created: Nov 2016
!
! *********************************************************************
subroutine wpnorm(wgtval, wgtvec, datvec, err, message, p_exp_in )
  implicit none
  ! input
  real(dp),               intent(in)  :: wgtvec(:)          ! weight vector
  real(dp),               intent(in)  :: datvec(:)          ! data value vector
  real(dp),    optional,  intent(in)  :: p_exp_in           ! p used in p-norm
  ! output 
  real(dp),               intent(out) :: wgtval             ! weighted (scaled) value
  integer(i4b),           intent(out) :: err                ! error code
  character(len=strLen),  intent(out) :: message            ! error message for current routine
  ! local 
  real(dp)                            :: wgtval_sum         ! Sum of weight (should be one)
  real(dp),parameter                  :: wgtMin=1.e-50_dp   ! minimum value for weight 
  real(dp),parameter                  :: paramMin=1.e-50_dp ! minimum value for parameters (for now exclued missing value, -999) 
  logical(lgc),allocatable            :: mask(:)            ! maks
  real(dp),allocatable                :: wgtvec_packed(:)   ! packed weight vector
  real(dp),allocatable                :: datvec_packed(:)   ! packed data value vector 
  real(dp)                            :: wgtvec_sum         ! packed weight vector
  real(dp)                            :: p_exp              ! local p used in p-norm
  integer(i4b)                        :: nElm_org           ! number of vector elements -original vector
  integer(i4b)                        :: nElm               ! number of vector elements -packed vector 
  integer(i4b)                        :: iElm               ! index of vector 

  err=0; message="wpnorm/"
  ! initialize p
  p_exp = 1._dp ! Default is weighted arithmetic mean
  if (present(p_exp_in)) p_exp = p_exp_in
  ! Compute size of dimension
  nElm_org=size(wgtvec)
  if (nElm_org /= size(datvec))then; err=20;message=trim(message)//'data and wgtvec vector:different size';return;endif
  ! Create mask 
  allocate(mask(nElm_org),stat=err); if(err/=0)then;message=trim(message)//'problem allocating mask';return;endif
  mask=(wgtvec > wgtMin .and. datvec > paramMin)
  ! Pack vector
  allocate(wgtvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'problem allocating wgtvec_packed';return;endif
  allocate(datvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'problem allocating datvec_packed';return;endif
  wgtvec_packed=pack(wgtvec,mask)
  datvec_packed=pack(datvec,mask)
  ! Re-check size of packed vector
  nElm=size(wgtvec_packed)
  if (nElm /= size(datvec_packed))then;err=20;message=trim(message)//'data and wgtvec vectors: different size';return;endif
  if (nElm > 0) then
    ! Recompute weight
    wgtvec_sum = sum(wgtvec_packed)
    if (wgtvec_sum /= 1.0) wgtvec_packed = wgtvec_packed / wgtvec_sum 
    wgtval_sum = 0._dp
    if (p_exp .gt. 100._dp) then
!      print*, '    using maximum norm'
      wgtval = maxval(datvec_packed(1:nElm))
    else if (abs(p_exp) .lt. epsilon(1._dp)) then
!      print*, '    using geometric mean'
      do iElm = 1, nElm
        wgtval_sum = wgtval_sum + log(datvec_packed(iElm)) * wgtvec_packed(iElm)
      end do 
      wgtval = exp(wgtval_sum)
    else if (p_exp .lt. -100._dp) then
!      print*, '    using minimum '
      wgtval = minval(datvec_packed(1:nElm))
    else
!      print*, '    using p-norm with exponent: ', p_exp
      do iElm = 1, nElm
        !wgtval_sum = wgtval_sum + exp(log(abs(datvec_packed(iElm))) * p_exp) * wgtvec_packed(iElm)
        wgtval_sum = wgtval_sum + (abs(datvec_packed(iElm))**p_exp)*wgtvec_packed(iElm)
      end do
      !wgtval = exp(log(wgtval_sum) * (1._dp / p_exp))
      wgtval = wgtval_sum**(1._dp / p_exp)
    end if
  else
    wgtval=dmiss
  end if

end subroutine

! *********************************************************************
! subroutine: computing unweighted sum 
! *********************************************************************
subroutine asum(wgtval, wgtvec, datvec, err, message )
  implicit none
  ! input
  real(dp),               intent(in)      :: wgtvec(:)           ! weight vector
  real(dp),               intent(in)      :: datvec(:)           ! data value vector
  ! output 
  real(dp),               intent(out)     :: wgtval              ! weighted (scaled) value
  integer(i4b),           intent(out)     :: err                 ! error code
  character(len=strLen),  intent(out)     :: message             ! error message for current routine
  ! local 
  real(dp),parameter                      :: wgtMin=1.e-50_dp    ! minimum value for weight 
  real(dp),parameter                      :: paramMin=1.e-50_dp  ! minimum value for parameters (for now exclued missing value, -999) 
  logical(lgc),allocatable                :: mask(:)             ! maks
  real(dp),allocatable                    :: wgtvec_packed(:)    ! packed weight vector
  real(dp),allocatable                    :: datvec_packed(:)    ! packed data value vector 
  integer(i4b)                            :: nElm_org            ! number of vector elements -original vector
  integer(i4b)                            :: nElm                ! number of vector elements -packed vector 

  err=0; message="asum/"
  ! check : sum of weight should be one
  nElm_org=size(wgtvec)
  if (nElm_org /= size(datvec))then; err=20;message=trim(message)//'data and wgtvec vector:different size';return;endif
  ! Create mask 
  allocate(mask(nElm_org),stat=err); if(err/=0)then;message=trim(message)//'problem allocating mask';return;endif
  mask=(wgtvec > wgtMin .and. datvec > paramMin)
  ! Pack vector
  allocate(wgtvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'problem allocating wgtvec_packed';return;endif
  allocate(datvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'problem allocating datvec_packed';return;endif
  wgtvec_packed=pack(wgtvec,mask)
  datvec_packed=pack(datvec,mask)
  ! Re-check size of packed vector
  nElm=size(wgtvec_packed)
  if (nElm /= size(datvec_packed))then;err=20;message=trim(message)//'data and wgtvec vectors: different size';return;endif
  if (nElm > 0) then
    wgtval=sum(datvec_packed)
  else
    wgtval=dmiss
  endif
  return 
end subroutine

! *********************************************************************
! subroutine: computing weighted median  
! *********************************************************************
subroutine wmedi(wgtval, wgtvec, datvec, err, message )
  implicit none
  ! input
  real(dp),dimension(:), intent(in)    :: wgtvec             ! original weight vector
  real(dp),dimension(:), intent(in)    :: datvec             ! original data vector 
  ! output 
  real(dp),              intent(out)   :: wgtval             ! weighted (scaled) value
  integer(i4b),          intent(out)   :: err             ! error code
  character(len=strLen), intent(out)   :: message          ! error message for current routine
  ! local 
  real(dp),parameter                   :: wgtMin=1.e-50_dp   ! minimum value for weight 
  real(dp)                             :: swapDat 
  real(dp)                             :: swapWgt 
  real(dp)                             :: wsum               ! sum of weight at ith sorted weight vector 
  integer(i4b)                         :: k                  ! counter 
  real(dp)                             :: wgtvec_sum         ! Sum of weight (should be one)
  logical(lgc),allocatable             :: mask(:)            ! maks
  real(dp),allocatable                 :: wgtvec_packed(:)   ! packed weight vector
  real(dp),allocatable                 :: datvec_packed(:)   ! packed data value vector 
  integer(i4b)                         :: nElm_org           ! number of vector elements -original vector
  integer(i4b)                         :: nElm 
  integer(i4b)                         :: iElm 
  integer(i4b)                         :: jElm 

  err=0; message="whmean/"
  nElm_org =size(wgtvec)
  if (nElm_org /= size(datvec))then; err=20;message=trim(message)//'data and wgtvec vector:different size';return;endif
  ! Create mask 
  allocate(mask(nElm_org),stat=err); if(err/=0)then;message=trim(message)//'error allocating mask';return;endif
  mask = (wgtvec > wgtMin .and. (datvec > dmiss .or. datvec < dmiss) )
  ! Pack vector
  allocate(wgtvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'error allocating wgtvec_packed';return;endif
  allocate(datvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'error allocating datvec_packed';return;endif
  wgtvec_packed = pack(wgtvec,mask)
  datvec_packed = pack(datvec,mask)
  ! Re-check size of packed vector
  nElm=size(wgtvec_packed)
  if (nElm /= size(datvec_packed))then;err=20;message=trim(message)//'data and wgtvec: different size';return;endif
  
  if (nElm > 0) then
    ! Recompute weight
    wgtvec_sum = sum(wgtvec_packed)
    if ( wgtvec_sum > 1.0_dp .and. wgtvec_sum < 1.0_dp ) &
      wgtvec_packed = wgtvec_packed/wgtvec_sum    
    !sorting data and wgtvec based on data values
    do iElm=1,nElm
      do jElm=iElm,nElm
        if(datvec_packed(jElm) < datvec_packed(iElm)) then
          swapDat = datvec_packed(iElm)
          datvec_packed(iElm) = datvec_packed(jElm)
          datvec_packed(jElm) = swapDat
          
          swapWgt = wgtvec_packed(iElm)
          wgtvec_packed(iElm) = wgtvec_packed(jElm)
          wgtvec_packed(jElm) = swapWgt
         end if 
      end do
    end do
    ! Initialize
    k = 1
    wsum = wgtvec_packed(1)
    ! Find weighted Median
    do 
      if (wsum > wgtvec_sum/2) exit 
      k = k+1
      wsum = wsum + wgtvec_packed(k)
    enddo
    wgtval = datvec_packed(k)
  else
    wgtval = dmiss 
  endif
  return
end subroutine

! *********************************************************************
! subroutine: computing weighted mode 
! *********************************************************************
subroutine wmode(wgtval, wgtvec, datvec, err, message )
  ! Define variables
  implicit none
  ! input
  real(dp),              intent(in)    :: wgtvec(:)             ! original weight vector
  real(dp),              intent(in)    :: datvec(:)             ! original data vector 
  ! output 
  real(dp),              intent(out)   :: wgtval             ! weighted (scaled) value
  integer(i4b),          intent(out)   :: err             ! error code
  character(len=strLen), intent(out)   :: message          ! error message for current routine
  ! local 
  real(dp),parameter                   :: wgtMin=1.e-50_dp   ! minimum value for weight 
  real(dp)                             :: wgtvec_sum         ! Sum of weight (should be one)
  real(dp),allocatable                 :: wsum(:)            ! sum of weight where data values are identical 
  logical(lgc),allocatable             :: mask(:)            ! masks
  logical(lgc),allocatable             :: maskunq(:)         ! mask to identify unique element
  real(dp),allocatable                 :: wgtvec_packed(:)   ! packed weight vector
  real(dp),allocatable                 :: datvec_packed(:)   ! packed data value vector 
  real(dp),allocatable                 :: datvec_unq(:)      ! unique data value vector 
  integer(i4b),allocatable             :: idxdat(:)          ! index vector indicating location of unique element
  integer(i4b)                         :: nElm_org           ! number of vector elements -original vector
  integer(i4b)                         :: nElm               ! number of vector elements masked
  integer(i4b)                         :: iElm               ! loop index of vector element
  integer(i4b)                         :: jElm               ! loop index of vector element
  integer(i4b)                         :: iMax(1)            ! index where element is max in vector

  err=0; message="wmode/"
  nElm_org =size(wgtvec)
  if (nElm_org /= size(datvec))then;err=20;message=trim(message)//'data and wgtvec vectors: different size';return;endif
  ! Create mask to eliminate element with no weight values
  allocate(mask(nElm_org),stat=err); if(err/=0)then;message=trim(message)//'error allocating mask';return;endif
  !mask = (wgtvec > wgtMin .and. datvec > paramMin )
  mask = (wgtvec > wgtMin )
  ! Pack vector
  allocate(wgtvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'error allocating wgtvec_packed';return;endif
  allocate(datvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'error allocating datvec_packed';return;endif
  wgtvec_packed = pack(wgtvec,mask)
  datvec_packed = pack(datvec,mask)
  ! Re-check size of packed vector
  nElm=size(wgtvec_packed)
  if (nElm /= size(datvec_packed))then;err=20;message=trim(message)//'data and wgtvec: different size';return;endif
    
  if (nElm > 0) then
    ! Recompute weight
    wgtvec_sum = sum(wgtvec_packed)
    if (wgtvec_sum /= 1.0) &
      wgtvec_packed = wgtvec_packed/wgtvec_sum

    ! Find unique element in data vector
    allocate(maskunq(nElm),stat=err); if(err/=0)then;message=trim(message)//'error allocating maskunq';return;endif
    maskunq = .true.
    do iElm = nElm,2,-1
      maskunq(iElm) = .not.(any(int(datvec_packed(1:iElm-1))==int(datvec_packed(iElm))))
    end do
    ! make index vector
    allocate(idxdat, source=pack([(iElm,iElm=1,nElm)],maskunq))
    ! copy the unique element in datvec_unq
    allocate(datvec_unq, source=datvec_packed(idxdat))
    allocate(wsum(size(datvec_unq)),stat=err);if(err/=0)then;message=trim(message)//'error allocating wsum';return;endif
    wsum = 0._dp
    do iElm=1,size(datvec_unq)
      do jElm=1,nElm
        if (datvec_packed(jElm) == datvec_unq(iElm)) & 
          wsum(iElm) =  wsum(iElm) + wgtvec_packed(jElm)
      end do
    end do
    iMax=maxloc(wsum)
    wgtval = datvec_unq(imax(1))
  else
    wgtval = int(dmiss)
  endif
  return
end subroutine

! *********************************************************************
! subroutine: computing weighted harmonic mean 
! *********************************************************************
subroutine whmean(wgtval, wgtvec, datvec, err, message )
  implicit none
  ! input
  real(dp),               intent(in)     :: wgtvec(:)           ! Original weight vector
  real(dp),               intent(in)     :: datvec(:)           ! Original data value vector
  ! output 
  real(dp),               intent(out)    :: wgtval              ! weighted (scaled) value
  integer(i2b),           intent(out)    :: err             ! error code
  character(len=strLen),  intent(out)    :: message          ! error message for current routine
  ! local 
  real(dp),parameter                     :: wgtMin=1.e-50_dp    ! minimum value for weight 
  real(dp)                               :: wgtval_sum          ! Sum of weight (should be one)
  logical(lgc),allocatable               :: mask(:)             ! maks
  real(dp),allocatable                   :: wgtvec_packed(:)    ! packed weight vector
  real(dp),allocatable                   :: datvec_packed(:)    ! packed data value vector 
  real(dp)                               :: wgtvec_sum          ! Sum of weight (should be one)
  integer(i4b)                           :: nElm_org            ! number of vector elements -original vector
  integer(i4b)                           :: nElm                ! Number of vector elements 
  integer(i4b)                           :: iElm                ! Index of vector 

  err=0; message="whmean/"
  ! check : sum of weight should be one
  ! check : sum of datvec is positive real 
  nElm_org =size(wgtvec)
  if (nElm_org /= size(datvec))then; err=20;message=trim(message)//'data and wgtvec vector:different size';return;endif
  ! Create mask 
  allocate(mask(nElm_org),stat=err); if(err/=0)then;message=trim(message)//'error allocating mask';return;endif
  mask = (wgtvec > wgtMin .and. (datvec > dmiss .or. datvec < dmiss))
  ! Pack vector
  allocate(wgtvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'error allocating wgtvec_packed';return;endif
  allocate(datvec_packed(count(mask)),stat=err); if(err/=0)then;message=trim(message)//'error allocating datvec_packed';return;endif
  wgtvec_packed = pack(wgtvec,mask)
  datvec_packed = pack(datvec,mask)
  ! Re-check size of packed vector
  nElm=size(wgtvec_packed)
  if (nElm /= size(datvec_packed))then;err=20;message=trim(message)//'data and wgtvec: different size';return;endif

  if (nElm > 0) then
    ! Recompute weight
    wgtvec_sum = sum(wgtvec_packed)
    if (wgtvec_sum /= 1.0) &
      wgtvec_packed = wgtvec_packed/wgtvec_sum    
    if ( any(datvec_packed == 0._dp) ) then
      wgtval = 0._dp
    else
      wgtval_sum =0._dp
      overlap_loop: do iElm=1,nElm
        wgtval_sum=wgtval_sum+wgtvec_packed(iElm)/datvec_packed(iElm)
      enddo overlap_loop
      wgtval = 1._dp/wgtval_sum
    endif
  else
    wgtval = dmiss 
  endif
  return
end subroutine


end module upscaling
