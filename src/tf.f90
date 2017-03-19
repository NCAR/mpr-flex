module tf
! Compute model soil parameter transfer function 
use nrtype                                        ! variable types, etc.
use data_type                                     ! Including custum data structure definition
use public_var                                     ! Including common constant (physical constant, other e.g., missingVal, etc.)
use var_lookup, only:ixVarSoilData, ixVarVegData, ixBeta, ixGamma, nBeta

implicit none

private

public::comp_model_param

contains

! ********************************************************************************************
! Public subroutine: Execute computation of model-dependant soil parameter transfer function  
! *********************************************************************************************
subroutine comp_model_param(parSxySz,          &  ! in/output: soil parameter values for all L0 polygons that are included
                            parVxy,            &  ! in/output: veg parameter values for all L0 polygons that are included 
                            sdata,             &  ! input: soil data
                            vdata,             &  ! input: vege data
                            gammaParMasterMeta,&  ! input: gamma parameter meta file - val of calibrating parammeter is adjusted via calibration 
                            nSLyr,             &  ! input: number of soil layers
                            nSPoly,            &  ! input: number of soil polygons
                            nVPoly,            &  ! input: number of vege polygons (to be deleted as nVpoly=nSpoly)
                            err,message) 

  use globalData, only:betaMaster, soilBetaInGamma, vegBetaInGamma, betaNeeded
  use get_ixname, only:get_ixBeta
  implicit none
  ! in/out
  type(namedvar2),      intent(inout) :: parSxySz(:)            ! soil parameter values for ParSxySz(:)%varDat(lyr,poly) 
  type(namedvar2),      intent(inout) :: parVxy(:)              ! veg parameter values for ParVxy(:)%varDat(lyr,poly) 
  ! input
  type(namevar),        intent(in)    :: sdata(:)               ! storage of soil data strucuture
  type(namevar),        intent(in)    :: vdata(:)               ! storage of veg data strucuture
  type(par_meta)                      :: gammaParMasterMeta(:)
  integer(i4b),         intent(in)    :: nSLyr                  ! number of soil layer
  integer(i4b),         intent(in)    :: nSPoly                 ! number of soil polygons
  integer(i4b),         intent(in)    :: nVPoly                 ! number of vege polygons (to be deleted)
  ! output
  integer(i4b),         intent(out)   :: err                   ! error code
  character(len=strLen),intent(out)   :: message                ! error message for current routine
  ! Local 
  type(namedvar2)                     :: parTemp(nBeta)          ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  integer(i4b)                        :: ix                     ! index of gamma parameter 
  integer(i4b)                        :: idBeta                 ! id of beta parameter array 
  integer(i4b)                        :: iParm                  ! Loop index of model parameters (e.g., VIC)
  logical(lgc)                        :: checkDone(nBeta)        ! used to check if the VIC parameter is processed

  err=0; message="comp_model_param/"
  first: associate(sclass  => sdata(ixVarSoilData%soilclass)%ivar2,   &
                   hslyrs  => sdata(ixVarSoilData%hslyrs)%dvar2,      &
                   sand    => sdata(ixVarSoilData%sand_frc)%dvar2,    & 
                   silt    => sdata(ixVarSoilData%silt_frc)%dvar2,    & 
                   clay    => sdata(ixVarSoilData%clay_frc)%dvar2,    & 
                   bulkd   => sdata(ixVarSoilData%bulk_density)%dvar2,& 
                   elestd  => sdata(ixVarSoilData%ele_std)%dvar1,     & 
                   elemean => sdata(ixVarSoilData%ele_mean)%dvar1,    & 
                   slpmean => sdata(ixVarSoilData%slp_mean)%dvar1,    &
                   monLai  => vdata(ixVarVegData%lai)%dvar2,          &
                   gammaPar=> gammaParMasterMeta(:)%val)
  do iParm = 1,size(betaNeeded)
    ix = get_ixBeta(betaNeeded(iParm)) 
    if (trim(betaMaster(ix)%ptype)=='soil')then
      allocate(parTemp(ix)%varData(nSLyr,nSPoly) ,stat=err); if(err/=0)then;message=trim(message)//'error allocating parTemp';stop;endif
    elseif (betaMaster(ix)%ptype=='veg')then 
      allocate(parTemp(ix)%varData(nMonth,nVPoly) ,stat=err); if(err/=0)then;message=trim(message)//'error allocating parTemp';stop;endif
    endif
    second: associate (xPar => parTemp(ix)%varData, &
                       tfid => betaMaster(ix)%tftype)
    if (tfid==-999_i2b) tfid=1_i2b
    select case(ix)
      case(ixBeta%ks)
        checkDone(ix)=.true.
        xPar = ks( sand, clay, gammaPar, tfid)
      case(ixBeta%bd)
        checkDone(ix)=.true.
        xPar  = bd( bulkd, gammaPar,tfid )
      case(ixBeta%phi)
        if(.not.checkDone(ixBeta%bd)) then;err=10;message=trim(message)//'need to process bd before phi';return;endif
        checkDone(ix)=.true.
        xPar = phi( sand, clay, parTemp(ixBeta%bd)%varData,gammaPar, tfid)
      case(ixBeta%b)
        checkDone(ix)=.true.
        xPar = ret_curve( sand, clay, gammaPar, tfid)
      case(ixBeta%psis)
        checkDone(ix)=.true.
        xPar = psis( sand, silt, gammaPar, tfid)
      case(ixBeta%fc)
        if(.not.checkDone(ixBeta%psis)) then;err=10;message=trim(message)//'need to process psis before fc';return;endif
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before fc';return;endif
        if(.not.checkDone(ixBeta%b))    then;err=10;message=trim(message)//'need to process b before fc';return;endif
        checkDone(ix)=.true.
        xPar = fc(sand, parTemp(ixBeta%phi)%varData, parTemp(ixBeta%psis)%varData, parTemp(ixBeta%b)%varData,gammaPar, tfid)
      case(ixBeta%wp)
        if(.not.checkDone(ixBeta%psis)) then;err=10;message=trim(message)//'need to process psis before wp';return;endif
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before wp';return;endif
        if(.not.checkDone(ixBeta%b))    then;err=10;message=trim(message)//'need to process b before wp';return;endif
        checkDone(ix)=.true.
        xPar = wp( parTemp(ixBeta%phi)%varData, parTemp(ixBeta%psis)%varData, parTemp(ixBeta%b)%varData, gammaPar, tfid) 
      case(ixBeta%myu)
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before myu';return;endif
        if(.not.checkDone(ixBeta%fc))   then;err=10;message=trim(message)//'need to process fc before myu';return;endif
        checkDone(ix)=.true.
        xPar= myu( parTemp(ixBeta%phi)%varData, parTemp(ixBeta%fc)%varData, gammaPar, tfid)
      case(ixBeta%binfilt)
        checkDone(ix)=.true.
        xPar=spread( infilt( elestd, gammaPar, tfid), 1, nSLyr)
      case(ixBeta%D1)
        if(.not.checkDone(ixBeta%ks))  then;err=10;message=trim(message)//'need to process "ks" before "D1"';return;endif
        if(.not.checkDone(ixBeta%phi)) then;err=10;message=trim(message)//'need to process "phi" before "D1"';return;endif
        checkDone(ix)=.true. 
        xPar=D1( slpmean,                    &
                 parTemp(ixBeta%ks)%varData,  &
                 parTemp(ixBeta%phi)%varData, & 
                 hslyrs,                     &
                 gammaPar,                   &
                 tfid)
      case(ixBeta%Ds)
        if(.not.checkDone(ixBeta%D1))    then;err=10;message=trim(message)//'need to process "D1" before "Ds"';return;endif
        if(.not.checkDone(ixBeta%D3))    then;err=10;message=trim(message)//'need to process "D3" before "Ds"';return;endif
        if(.not.checkDone(ixBeta%Dsmax)) then;err=10;message=trim(message)//'need to process "Dsmax" before "Ds"';return;endif
        checkDone(ix)=.true. 
        xPar=Ds( parTemp(ixBeta%D1)%varData, parTemp(ixBeta%D3)%varData, parTemp(ixBeta%Dsmax)%varData, tfid)
      case(ixBeta%D4)
        checkDone(ix)=.true.
        xPar=D4(gammaPar, tfid)
      case(ixBeta%c)
        if(.not.checkDone(ixBeta%D4)) then;err=10;message=trim(message)//'need to process "D4" before "c"';return;endif
        checkDone(ix)=.true.
        xPar=cexpt(parTemp(ixBeta%D4)%varData, tfid)
      case(ixBeta%SD)
        checkDone(ix)=.true.
        xPar=soilDensity(parTemp(ixBeta%sd)%varData, gammaPar ,tfid)
      case(ixBeta%expt)
        if(.not.checkDone(ixBeta%b)) then;err=10;message=trim(message)//'need to process "b" before "expt"';return;endif
        checkDone(ix)=.true.
        xPar=expt( parTemp(ixBeta%b)%varData, gammaPar ,tfid)
      case(ixBeta%D2)
        if(.not.checkDone(ixBeta%ks)) then;err=10;message=trim(message)//'need to process "ksat" before "D2"';return;endif
        if(.not.checkDone(ixBeta%D4)) then;err=10;message=trim(message)//'need to process "D4" before "D2"';return;endif
        checkDone(ix)=.true.
        xPar=D2( slpmean, parTemp(ixBeta%ks)%varData, parTemp(ixBeta%D4)%varData, gammaPar, tfid )
      case(ixBeta%Dsmax)
        if(.not.checkDone(ixBeta%D1)) then;err=10;message=trim(message)//'need to process "D1" before "Dsmax"';return;endif
        if(.not.checkDone(ixBeta%D2)) then;err=10;message=trim(message)//'need to process "D2" before "Dsmax"';return;endif
        if(.not.checkDone(ixBeta%D3)) then;err=10;message=trim(message)//'need to process "D3" before "Dsmax"';return;endif
        if(.not.checkDone(ixBeta%c))  then;err=10;message=trim(message)//'need to process "c" before "Dsmax"';return;endif
        checkDone(ix)=.true. 
        xPar=Dsmax( parTemp(ixBeta%D1)%varData,   &
                    parTemp(ixBeta%D2)%varData,   & 
                    parTemp(ixBeta%D3)%varData,   & 
                    parTemp(ixBeta%c)%varData,    & 
                    parTemp(ixBeta%phi)%varData,  &
                    hslyrs,                      &
                    tfid)
      case(ixBeta%bbl)
        if(.not.checkDone(ixBeta%expt)) then;err=10;message=trim(message)//'need to process "expt" before "bubble"';return;endif
        checkDone(ix)=.true.
        xPar=bubble( parTemp(ixBeta%expt)%varData, gammaPar, tfid)
      case(ixBeta%WcrFrac)
        if(.not.checkDone(ixBeta%fc))then;err=10;message=trim(message)//'need to process "fc" before "WcrFrac"';return;endif 
        if(.not.checkDone(ixBeta%phi))then;err=10;message=trim(message)//'need to process "phi" before "WcrFrac"';return;endif 
        checkDone(ix)=.true.
        xPar=WcrFrac( parTemp(ixBeta%fc)%varData,parTemp(ixBeta%phi)%varData,gammaPar , tfid)
      case(ixBeta%WpwpFrac)
        if(.not.checkDone(ixBeta%wp))then;err=10;message=trim(message)//'need to process "wp" before "WpwpFrac"';return;endif 
        if(.not.checkDone(ixBeta%phi))then;err=10;message=trim(message)//'need to process "phi" before "WpwpFrac"';return;endif 
        checkDone(ix)=.true.
        xPar=WpwpFrac( parTemp(ixBeta%wp)%varData, parTemp(ixBeta%phi)%varData, gammaPar, tfid)
      case(ixBeta%D3)
        if(.not.checkDone(ixBeta%fc))   then;err=10;message=trim(message)//'need to process fc before D3';return;endif
        checkDone(ix)=.true.
        xPar=D3( parTemp(ixBeta%fc)%varData, hslyrs, gammaPar , tfid)
      case(ixBeta%Ws)
        if(.not.checkDone(ixBeta%D3))then;err=10;message=trim(message)//'need to process "D3" before "Dsmax"';return;endif 
        checkDone(ix)=.true. 
        xPar=Ws( parTemp(ixBeta%D3)%varData, parTemp(ixBeta%phi)%varData,hslyrs , tfid)
      case(ixBeta%twm)
        if(.not.checkDone(ixBeta%fc))   then;err=10;message=trim(message)//'need to process fc before twm';return;endif
        if(.not.checkDone(ixBeta%wp))   then;err=10;message=trim(message)//'need to process wp before twm';return;endif
        checkDone(ix)=.true. 
        xPar= twm( parTemp(ixBeta%fc)%varData,parTemp(ixBeta%wp)%varData,hslyrs , tfid)
      case(ixBeta%fwm)
        if(.not.checkDone(ixBeta%fc))   then;err=10;message=trim(message)//'need to process fc before fwm';return;endif
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before fwm';return;endif
        checkDone(ix)=.true. 
        xPar= fwm( parTemp(ixBeta%phi)%varData,parTemp(ixBeta%fc)%varData,hslyrs, tfid )
      case(ixBeta%fsm)
        if(.not.checkDone(ixBeta%fwm))then;err=10;message=trim(message)//'need to process "fwm" before "fsm"';return;endif 
        checkDone(ix)=.true. 
        xPar= fsm( parTemp(ixBeta%fwm)%varData,parTemp(ixBeta%phi)%varData,parTemp(ixBeta%wp)%varData,gammaPar, tfid )
      case(ixBeta%fpm)
        if(.not.checkDone(ixBeta%fwm))then;err=10;message=trim(message)//'need to process "fwm" before "fpm"';return;endif 
        if(.not.checkDone(ixBeta%fsm))then;err=10;message=trim(message)//'need to process "fsm" before "fpm"';return;endif 
        checkDone(ix)=.true. 
        xPar= fpm( parTemp(ixBeta%fwm)%varData, parTemp(ixBeta%fsm)%varData, tfid )
      case(ixBeta%zk)
        if(.not.checkDone(ixBeta%fc))   then;err=10;message=trim(message)//'need to process fc before zk';return;endif
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before zk';return;endif
        checkDone(ix)=.true. 
        xPar= zk( parTemp(ixBeta%phi)%varData,parTemp(ixBeta%fc)%varData,gammaPar, tfid )
      case(ixBeta%zsk)
        if(.not.checkDone(ixBeta%fc))   then;err=10;message=trim(message)//'need to process fc before zsk';return;endif
        if(.not.checkDone(ixBeta%wp))   then;err=10;message=trim(message)//'need to process wp before zsk';return;endif
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before zsk';return;endif
        checkDone(ix)=.true. 
        xPar= zsk( parTemp(ixBeta%phi)%varData,parTemp(ixBeta%fc)%varData,parTemp(ixBeta%wp)%varData,gammaPar, tfid )
      case(ixBeta%zpk)
        if(.not.checkDone(ixBeta%ks))   then;err=10;message=trim(message)//'need to process ks before zpk';return;endif
        if(.not.checkDone(ixBeta%myu))  then;err=10;message=trim(message)//'need to process myu before zpk';return;endif
        checkDone(ix)=.true. 
        xPar= zpk( parTemp(ixBeta%ks)%varData,parTemp(ixBeta%myu)%varData,hslyrs,gammaPar, tfid )
      case(ixBeta%pfree)
        if(.not.checkDone(ixBeta%phi))  then;err=10;message=trim(message)//'need to process phi before pfree';return;endif
        if(.not.checkDone(ixBeta%wp))   then;err=10;message=trim(message)//'need to process wp before pfree';return;endif
        checkDone(ix)=.true. 
        xPar= pfree( parTemp(ixBeta%phi)%varData,parTemp(ixBeta%wp)%varData,gammaPar, tfid )
      case(ixBeta%zperc)
        if(.not.checkDone(ixBeta%twm))then;err=10;message=trim(message)//'need to process "twm" before "pfree"';return;endif 
        if(.not.checkDone(ixBeta%fsm))then;err=10;message=trim(message)//'need to process "fsm" before "pfree"';return;endif 
        if(.not.checkDone(ixBeta%zsk))then;err=10;message=trim(message)//'need to process "zsk" before "pfree"';return;endif 
        if(.not.checkDone(ixBeta%fpm))then;err=10;message=trim(message)//'need to process "fpm" before "pfree"';return;endif 
        if(.not.checkDone(ixBeta%zpk))then;err=10;message=trim(message)//'need to process "zpk" before "pfree"';return;endif 
        checkDone(ix)=.true. 
        xPar= zperc(parTemp(ixBeta%twm)%varData,parTemp(ixBeta%fsm)%varData,parTemp(ixBeta%zsk)%varData,parTemp(ixBeta%fpm)%varData,parTemp(ixBeta%zsk)%varData, tfid)
      case(ixBeta%rexp)
        if(.not.checkDone(ixBeta%wp))then;err=10;message=trim(message)//'need to process "wp" before "rexp"';return;endif 
        checkDone(ix)=.true. 
        xPar= rexp( parTemp(ixBeta%wp)%varData,gammaPar, tfid )
      case(ixBeta%lai)
        xPar= lai( monLai, gammaPar, tfid )
    end select ! end of parameter case
    end associate second
  end do ! end of parameter loop
  end associate first
  ! extract beta parameters in 'CalPar' list
  do iParm=1,size(soilBetaInGamma)
    idBeta=get_ixBeta(trim(soilBetaInGamma(iParm)))
    parSxySz(iParm)%varData=parTemp(idBeta)%varData 
  enddo
  do iParm=1,size(vegBetaInGamma)
    idBeta=get_ixBeta(trim(vegBetaInGamma(iParm)))
    parVxy(iParm)%varData=parTemp(idBeta)%varData 
  enddo
  return
end subroutine 

! *********************************************************************
! Library of hydrologic model transfer functions
! *********************************************************************
! infilt parameter 
! *********************************************************************
function infilt(elestd_in, gammaPar, opt)
  ! Use Arno scheme 
  ! b = (sigma_elev-sigma_min)/(sigma_ele+sigma_max)
  ! where
  ! sigma_ele: standard deviation of elevation in polygon
  ! sigma_min: parameter - depending on spatial scale of operation 
  ! sigma_max: parameter - depending on spatial scale of operation 
  !
  ! Reference: 
  ! van den Hurk, B. J. J. M., and P. Viterbo, 2003: The Torne-Kalix PILPS 2(e) experiment  
  ! as a test bedfor modifications to the ECMWF land surface scheme. Global Planet. Change, 38, 165?173
  !
  ! Dümenil, L., and E. Todini, 1992: A rainfall-runoff scheme for use in the Hamburg climate model. 
  ! Advances in Theoretical Hydrology: A Tribute to James Dooge. J. P. O?Kane, Ed., 
  ! European Geophysical Society Series on Hydrological Sciences, Vol. 1, Elsevier, 129-157. 

  implicit none
  ! input
  real(dp), intent(in)    :: elestd_in(:)
  real(dp), intent(in)    :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: infilt(size(elestd_in)) 
  ! local 
  real(dp),parameter      :: infilt_min=0.03_dp
  real(dp),parameter      :: infilt_max=0.50_dp
  character(len=strLen)   :: message                           
  
  message="infilt/"
  associate(g1=>gammaPar(ixGamma%binfilt1gamma1), &
            g2=>gammaPar(ixGamma%binfilt1gamma2))
  select case(opt)
    case(1); 
      where ( elestd_in /= dmiss ) 
        infilt = (log(elestd_in+verySmall)-g1)/(log(elestd_in+verySmall)+g2*10) !transfer function 
      else where
        infilt = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  !cap value with upper and lower bounds 
  where ( infilt > infilt_max ) infilt=infilt_max 
  where ( infilt /= dmiss .and. infilt < infilt_min ) infilt=infilt_min 
  return
end function

! *********************************************************************
! residual_moist parameter 
! *********************************************************************
function residMoist(opt)
  implicit none
  ! input
   integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp) :: residMoist
  ! local 
  character(len=strLen)   :: message                           
 
  message="residMoist/"
  select case(opt)
    case(1); 
      residMoist = 0._dp
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
 return
end function
         
! ***********
!  Nijssen basiflow D1 parameter
! *********************************************************************
function D1(slope_in, ks_in, phi_in, h_in, gammaPar, opt)
  implicit none
  ! input
  real(dp), intent(in)    :: slope_in(:)   ! mean slope [%]
  real(dp), intent(in)    :: ks_in(:,:)    ! Ksat [mm/s]
  real(dp), intent(in)    :: phi_in(:,:)   ! porosity [-]
  real(dp), intent(in)    :: h_in(:,:)     ! layer thickness [m]
  real(dp), intent(in)    :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: D1(size(ks_in,1),size(ks_in,2)) !output: [day-1]
  ! local 
  real(dp),allocatable    :: slope2d(:,:)
  integer(i4b)            :: n1            ! number of 1st dimension 
  integer(i4b)            :: n2            ! number of 1st dimension 
  real(dp), allocatable   :: S(:,:)              ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp), parameter     :: D1_min=0.0001_dp
  real(dp), parameter     :: D1_max=1.0_dp
  character(len=strLen)   :: message                           
   
  message="D1/"
  ! local variable allocation
  n1=size(D1,1)
  n2=size(D1,2)
  allocate(S(n1,n2)) 
  allocate(slope2d(n1,n2)) 
  S=phi_in*h_in
  S=1.0_dp
  slope2d=spread(slope_in,1,n1)
  associate(g1=>gammaPar(ixGamma%D11gamma1))
  ! compute parameters 
  select case(opt)
    case(1); 
      where ( slope2d /= dmiss .and. Ks_in /= dmiss )
        D1 = S**(-1)*10**(-1*g1)*Ks_in*(60*60*24)*(slope2d*0.01)
      else where
        D1 = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate 
  ! cap value with upper and lower bounds 
  where ( D1 > D1_max ) D1=D1_max 
  where ( D1 > 0._dp .and. D1 < D1_min ) D1=D1_min 
  return
end function

! ***********
!  Arno basiflow Ds parameter
! *********************************************************************
function Ds( D1, D3, Dsmax, opt)
 implicit none
  ! input
  real(dp), intent(in)     :: D1(:,:)       ! nijssen baseflow D1 parameter [day-1]
  real(dp), intent(in)     :: D3(:,:)       ! nijssen baseflow D3 parameter [mm]
  real(dp), intent(in)     :: Dsmax(:,:)    ! ARNO Dsmax parameter [mm/day]
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                 :: Ds(size(D1,1),size(D1,2)) !output: [-]
  ! local 
  real(dp), parameter      :: Ds_min=0.0001_dp
  real(dp), parameter      :: Ds_max=1.0_dp
  character(len=strLen)   :: message                           

  message="Ds/"
  select case(opt)
    case(1); 
      where ( D1 /= dmiss .and. D3 /= dmiss .and. Dsmax /= dmiss )
        Ds = D1 * D3 / Dsmax
      else where
        Ds = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
 ! cap value with upper and lower bounds 
 where ( Ds > Ds_max ) Ds=Ds_max 
 where ( Ds > 0._dp .and. Ds < Ds_min ) Ds=Ds_min 
 return
end function

! ***********
! Nijssen baseflow D2 parameter
! *********************************************************************
function D2(slope_in, ks_in, D4_in, gammaPar, opt) 
  implicit none
  ! input
  real(dp), intent(in)    :: slope_in(:)   ! slope percent
  real(dp), intent(in)    :: Ks_in(:,:)    ! ksat [mm/s]
  real(dp), intent(in)    :: D4_in(:,:)    ! VIC D4 parameter [-]
  real(dp), intent(in)    :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)               :: D2(size(Ks_in,1),size(Ks_in,2)) ! output [day^-D4]
  ! local 
  real(dp),allocatable   :: slope2d(:,:)
  integer(i4b)           :: n1           ! number of element for 1st dimension
  integer(i4b)           :: n2           ! number of element for 2nd dimension
  real(dp), allocatable  :: S(:,:)             ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp), parameter    :: D2_min=0.0001_dp
  real(dp), parameter    :: D2_max=1.0_dp
  character(len=strLen)   :: message                           
  
  message="D2/"
 ! local variable allocation
  n1=size(D2,1)
  n2=size(D2,2)
  allocate(S(n1,n2))
  allocate(slope2d(n1,n2)) 
  S=1.0_dp
  slope2d=spread(slope_in,1,n1)
 ! compute parameters 
  associate(g1=>gammaPar(ixGamma%D21gamma1))
  select case(opt)
    case(1); 
       where ( slope2d /= dmiss .and. Ks_in /= dmiss )
         D2 = S**(-1*D4_in)*10**(-1*g1)*Ks_in*(60*60*24)*(slope2d*0.01)
       else where
         D2 = dmiss
       end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
 ! cap value with upper and lower bounds 
  where ( D2 > D2_max ) D2=D2_max
  where ( D2 > 0._dp .and. D2 < D2_min ) D2=D2_min
  return
end function

! ***********
! Arno baseflow Dsmax parameter
! *********************************************************************
function Dsmax( D1,           & ! input:  Nijssen baseflow D1 parameter [day^-1]
                D2,           & ! input:  Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
                D3,           & ! input:  Nijssen baseflow D3 parameter [mm]
                c,            & ! input:  c parameter [mm]
                phi_in,       & ! input:  porosity [cm^3/cm^-3]
                h_in,         & ! input:  Soil layer thickness [m]
                opt)          
  implicit none
  ! input
  real(dp), intent(in)    :: D1(:,:)        ! Nijssen baseflow D1 parameter [day^-1]
  real(dp), intent(in)    :: D2(:,:)        ! Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
  real(dp), intent(in)    :: D3(:,:)        ! Nijssen baseflow D3 parameter [mm]
  real(dp), intent(in)    :: c(:,:)         ! c parameter [mm]
  real(dp), intent(in)    :: phi_in(:,:)    ! porosity [cm^3/cm^-3]
  real(dp), intent(in)    :: h_in(:,:)      ! Soil layer thickness [m]
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: Dsmax(size(D1,1),size(D1,2))     ! Dsmax parameter for Arno baseflow [mm day-1]
  ! local 
  real(dp), parameter     :: Dsmax_min=0.1_dp
  real(dp), parameter     :: Dsmax_max=30.0_dp
  character(len=strLen)   :: message                           

  message="Dsmax/"
  select case(opt)
    case(1); 
      where ( phi_in /= dmiss .and. h_in /= dmiss )
        Dsmax = D2*(phi_in*h_in*1000-D3)**c+D1*(phi_in*h_in*1000)
      else where
        Dsmax = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
 ! cap value with upper and lower bounds 
  where ( Dsmax > Dsmax_max ) Dsmax=Dsmax_max
  where ( Dsmax > 0._dp .and. Dsmax < Dsmax_min ) Dsmax=Dsmax_min
  return
end function

! ***********
!  Nijssen baseflow D3 parameter
! *********************************************************************
function D3( fc_in, h_in, gammaPar, opt) 
  implicit none
  ! input
  real(dp), intent(in)    :: fc_in(:,:)    ! input: field capacity [frac]
  real(dp), intent(in)    :: h_in(:,:)     ! input: layer thickness [m]
  real(dp), intent(in)    :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)              :: D3(size(fc_in,1),size(fc_in,2)) !output: [mm]
  ! local 
  real(dp), parameter   :: D3_min=0.0001_dp
  real(dp), parameter   :: D3_max=1000.0_dp
  character(len=strLen)   :: message                           

  message="D3/"
  associate(g1=>gammaPar(ixGamma%D31gamma1))
  select case(opt)
    case(1); 
      where ( fc_in /= dmiss .and. h_in /= dmiss ) 
        D3 = g1* fc_in * (h_in*1000)
      else where
        D3 = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
 ! cap value with upper and lower bounds 
  where ( D3 > D3_max ) D3=D3_max
  where ( D3 > 0._dp .and. D3 < D3_min ) D3=D3_min 
  return
end function

! ***********
!  Arno baseflow Ws parameter (conversion equation)
! *********************************************************************
function Ws( D3,         & ! input:  D3 parameter [mm]
             phi_in,     & ! input:  porosity [cm^3/cm^-3]
             h_in,       & ! input:  Soil layer thickness [m]
             opt)          
  implicit none
  ! input
  real(dp),    intent(in) :: D3(:,:)
  real(dp),    intent(in) :: phi_in(:,:)
  real(dp),    intent(in) :: h_in(:,:)
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: Ws(size(D3,1),size(D3,2))   ! output [-]
  ! local 
  real(dp), parameter     :: Ws_min=0.05_dp
  real(dp), parameter     :: Ws_max=1.0_dp
  character(len=strLen)   :: message                           

  message="Ws/"
  select case(opt)
    case(1); 
      where ( phi_in /= dmiss .and. h_in /= dmiss .and. D3 /= dmiss ) 
        Ws = D3 / phi_in / (h_in*1000)
      else where
        Ws = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
 ! cap value with upper and lower bounds 
  where ( Ws > Ws_max ) Ws=Ws_max
  where ( Ws > 0._dp .and. Ws < Ws_min ) Ws=Ws_min 
  return
end function

! ***********
!  Nijssen baseflow D4 parameter
! *********************************************************************
function D4( gammaPar,opt )
  implicit none
  ! input
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: D4
  ! local 
  character(len=strLen)   :: message                           
  
  message="D4/"
  associate(g1=>gammaPar(ixGamma%D41gamma1))
  select case(opt)
    case(1); 
      D4 = g1 
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function D4 

! ***********
!  c parameter
! *********************************************************************
function cexpt( D4, opt )
  implicit none
  ! input
  real(dp),    intent(in) :: D4(:,:)
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: cexpt(size(D4,1),size(D4,2))
  ! local 
  character(len=strLen)   :: message                           
  
  message="cexpt/"
  select case(opt)
    case(1); 
      cexpt = D4 
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  return
end function cexpt 

! ***********
! computing expt parameter 
! *********************************************************************
function expt( b_in, gammaPar, opt )
  implicit none
  ! input
  real(dp),    intent(in) :: b_in(:,:)     ! inpuy [-]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)               :: expt(size(b_in,1),size(b_in,2)) ! exponent in campbel equation [-]
  ! local 
  character(len=strLen)   :: message                           

  message="expt/"
  associate(g1=>gammaPar(ixGamma%exp1gamma1), &
            g2=>gammaPar(ixGamma%exp1gamma2))
  select case(opt)
    case(1); 
      where ( b_in /= dmiss ) 
        expt = g1+g2*b_in
      else where
        expt = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function expt 

! ************
! computing init_moist parameter  
! *********************************************************************
function initMoist( phi_in, h_in, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: phi_in(:,:)   ! porosity [-]
  real(dp),    intent(in) :: h_in(:,:)     ! thickness [m]
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: initMoist(size(phi_in,1),size(phi_in,2))
  ! local  
  character(len=strLen)   :: message                           
  
  message="initMoist/"
  select case(opt)
    case(1); 
      where ( phi_in /= dmiss ) 
        initMoist = phi_in*(h_in*1000.0_dp)
      else where
        initMoist = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  return
end function

! ***********
! bubble parameter 
! *********************************************************************
function bubble( expt_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: expt_in(:,:) 
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: bubble(size(expt_in,1),size(expt_in,2))
  ! local 
  character(len=strLen)   :: message                           
  
  message="bubble/"
  associate(g1=>gammaPar(ixGamma%bbl1gamma1), &
            g2=>gammaPar(ixGamma%bbl1gamma2))
  select case(opt)
    case(1); 
      where ( expt_in /= dmiss ) 
        bubble = g1*expt_in+g2
      else where
        bubble = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! ***********
! soil_density parameter 
! *********************************************************************
function soilDensity( srho_in, gammaPar, opt )
  implicit none
  ! input
  real(dp),    intent(in) :: srho_in(:,:)  ! input: kg/m^3]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: soilDensity(size(srho_in,1),size(srho_in,2))
  ! local 
  character(len=strLen)   :: message                           
  
  message="soilDensity/"
  associate(g1=>gammaPar(ixGamma%sd1gamma1))
  select case(opt)
    case(1); 
      where ( srho_in /= dmiss ) 
        soilDensity = g1*srho_in
      else where
        soilDensity = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! ***********
! WcrFrac parameter  
! *********************************************************************
function WcrFrac(fc_in, phi_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: fc_in(:,:)    ! input: field capacity [frac]
  real(dp),    intent(in) :: phi_in(:,:)   ! input: porosity [frac]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: WcrFrac(size(fc_in,1),size(fc_in,2)) !output: [frac]
  ! local 
  character(len=strLen)   :: message                           
  
  message="WcrFrac/"
  associate(g1=>gammaPar(ixGamma%WcrFrac1gamma1))
  select case(opt)
    case(1); 
      where ( fc_in /= dmiss .and. phi_in /= dmiss ) 
        WcrFrac = g1*fc_in/phi_in
      else where
        wcrFrac = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! ************
! WpwpFrac parameter  
! *********************************************************************
function WpwpFrac( wp_in, phi_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: wp_in(:,:)    ! Wilting point [frac]
  real(dp),    intent(in) :: phi_in(:,:)   ! Porosity [frac]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                :: WpwpFrac(size(wp_in,1),size(wp_in,2)) !output: [frac] 
  ! local 
  character(len=strLen)   :: message                           
  
  message="WpwFrac/"
  associate(g1=>gammaPar(ixGamma%WpwpFrac1gamma1))
  select case(opt)
    case(1); 
      where ( wp_in /= dmiss .and. phi_in /= dmiss ) 
        WpwpFrac = g1*wp_in/phi_in
      else where
        wpwpFrac = dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! ************
! TWM parameter (UZTWM for top layer and LZTWM for bottom layer)
! *********************************************************************
function twm( fc_in, wp_in, h_in, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: fc_in(:,:)    ! input: field capacity [-]           
  real(dp),    intent(in) :: wp_in(:,:)    ! input: wilting point [-]
  real(dp),    intent(in) :: h_in(:,:)     ! input: thickness of layer [mm]
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: twm(size(wp_in,1),size(wp_in,2)) ! output: tension water maximum [mm]
  character(len=strLen)   :: message                           

  message="twm/"
  select case(opt)
    case(1); 
      where ( wp_in/=dmiss .and. fc_in/=dmiss ) 
        twm=(fc_in-wp_in)*h_in*1000.0_dp  ! convert m to mm 
      else where
        twm=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  return
end function 

! *********************************************************************
! FWM parameter (UZFWM for top layer and LZFWM for bottom layer)
! *********************************************************************
function fwm( phi_in, fc_in, h_in, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: phi_in(:,:)   ! input: porosity [-]
  real(dp),    intent(in) :: fc_in(:,:)    ! input: field capacity [-]           
  real(dp),    intent(in) :: h_in(:,:)     ! input: thickness of layer [mm]
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: fwm(size(fc_in,1),size(fc_in,2)) ! output: free water maximum [mm]
  character(len=strLen)   :: message                           

  message="fwm/"
  select case(opt)
    case(1); 
      where ( phi_in/=dmiss .and. fc_in/=dmiss .and. h_in/=dmiss ) 
        fwm=(phi_in-fc_in)*h_in*1000.0_dp ! convert m to mm
      else where
        fwm=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  return
end function 

! *********************************************************************
! FSM parameter (LZFSM for bottom layer)
! *********************************************************************
function fsm( fwm_in, phi_in, wp_in, gammaPar ,opt)
  implicit none
  ! input
  real(dp),    intent(in) :: fwm_in(:,:)   ! Free water maximum 
  real(dp),    intent(in) :: phi_in(:,:)   ! Porosity
  real(dp),    intent(in) :: wp_in(:,:)    ! Wilting point
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: fsm(size(wp_in,1),size(wp_in,2)) ! output: supplementary free water maximum [mm]
  character(len=strLen)   :: message                           

  message="fsm/"
  associate(g1=>gammaPar(ixGamma%fsm1gamma1))
  select case(opt)
    case(1); 
      where ( phi_in/=dmiss .and. wp_in/=dmiss .and. fwm_in/=dmiss ) 
        fsm=fwm_in*(wp_in/phi_in)**g1
      else where
        fsm=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function 

! *********************************************************************
! FPM parameter (LZFPM for bottom layer)
! *********************************************************************
function fpm( fwm_in, fsm_in, opt )
  implicit none
  ! input
  real(dp),    intent(in) :: fwm_in(:,:)   ! Free water maximum 
  real(dp),    intent(in) :: fsm_in(:,:)   ! supplementary Free water maximum 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local  
  real(dp)                :: fpm(size(fwm_in,1),size(fwm_in,2)) ! output: primary free water maximum [mm]
  character(len=strLen)   :: message                           

  message="fpm/"
  select case(opt)
    case(1); 
      where ( fwm_in/=dmiss .and. fsm_in/=dmiss )
        fpm=fwm_in-fsm_in
      else where
        fpm=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  return
end function 

! *********************************************************************
! ZK parameter (UZK for top layer)
! *********************************************************************
function zk( phi_in, fc_in, gammaPar,opt)
  implicit none
  ! input
  real(dp),    intent(in) :: phi_in(:,:)   ! Porosity [frac]
  real(dp),    intent(in) :: fc_in(:,:)    ! field capacity [frac]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: zk(size(fc_in,1),size(fc_in,2)) ! output: draw coefficient from free water content [/day] 
  character(len=strLen)   :: message                           

  message="zk/"
  associate(g1=>gammaPar(ixGamma%zk1gamma1))
  select case(opt)
    case(1); 
      where ( phi_in/=dmiss .and. fc_in/=dmiss ) 
        zk=1.0_dp-(fc_in/phi_in)**g1
      else where
        zk=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function 

! *********************************************************************
! ZSK parameter (LZSK for top layer)
! *********************************************************************
function zsk( phi_in, fc_in, wp_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: phi_in(:,:)   ! Porosity [frac]
  real(dp),    intent(in) :: fc_in(:,:)    ! field capacity [frac]
  real(dp),    intent(in) :: wp_in(:,:)    ! Wilting point [frac]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: zsk(size(fc_in,1),size(fc_in,2)) ! output: draw coefficient from supplementary free water content [/day] 
  character(len=strLen)   :: message                           

  message="zsk/"
  associate(g1=>gammaPar(ixGamma%zsk1gamma1),&
            g2=>gammaPar(ixGamma%zsk1gamma2))
  select case(opt)
    case(1); 
      where ( phi_in/=dmiss .and. fc_in/=dmiss .and. wp_in/=dmiss ) 
        zsk=(1.0_dp-(fc_in/phi_in)**g1)/(1.0_dp+g2*(1.0_dp-wp_in)) 
      else where
        zsk=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function 

! *********************************************************************
! ZPK parameter (LZPK for top layer)
! *********************************************************************
function zpk( ks_in, h_in, myu_in, gammaPar,opt)
  implicit none
  ! input
  real(dp),    intent(in) :: ks_in(:,:)    ! Ksat [mm/s]
  real(dp),    intent(in) :: myu_in(:,:)   ! specific yield [-]
  real(dp),    intent(in) :: h_in(:,:)     ! input: thickness of layer [m]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: zpk(size(ks_in,1),size(ks_in,2)) ! output: draw coefficient from primary free water content [/day] 
  real(dp)                :: dt                               ! time step of the simulation [hr]
  character(len=strLen)   :: message                           

  message="zpk/"
  dt=24.0_dp ! unit: hr
  associate(g1=>gammaPar(ixGamma%zpk1gamma1))
  select case(opt)
    case(1); 
      where ( ks_in/=dmiss .and. myu_in/=dmiss .and. h_in/= dmiss ) 
        zpk=1.0_dp-exp(-1.0_dp*g1**2.0_dp*pi*ks_in*60.0_dp*h_in*1000.0_dp*dt/myu_in)
      else where
        zpk=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function 

! *********************************************************************
! PFREE parameter 
! *********************************************************************
function pfree( phi_in, wp_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in) :: phi_in(:,:)   ! input: porosity [frac]
  real(dp),    intent(in) :: wp_in(:,:)    ! input: wilting point [frac]
  real(dp),    intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b),intent(in) :: opt           ! input: option for transfer function form
  ! local 
  real(dp)                :: pfree(size(wp_in,1),size(wp_in,2)) ! output: tension water maximum [mm]
  character(len=strLen)   :: message                           

  message="pfree/"
  associate(g1=>gammaPar(ixGamma%pfree1gamma1))
  select case(opt)
    case(1); 
      where ( phi_in/=dmiss .and. wp_in/=dmiss ) 
        pfree=(wp_in/phi_in)**g1
      else where 
        pfree=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! *********************************************************************
! ZPERC parameter 
! *********************************************************************
function zperc( twm_in, fsm_in, zsk_in, fpm_in, zpk_in, opt)
  implicit none
  ! input
  real(dp),    intent(in)  :: twm_in(:,:)   ! input: tension water maximum [mm] 
  real(dp),    intent(in)  :: fsm_in(:,:)   ! input: supplemental free water maximum [mm] 
  real(dp),    intent(in)  :: zsk_in(:,:)   ! input: flow rate from supplemental free water maximum [day-1] 
  real(dp),    intent(in)  :: fpm_in(:,:)   ! input: primary free water maximum [mm] 
  real(dp),    intent(in)  :: zpk_in(:,:)   ! input: flow rate from primary free water maximum [day-1] 
  integer(i2b),intent(in)  :: opt           ! input: option for transfer function form
  ! output 
  real(dp)                 :: zperc(size(twm_in,1),size(twm_in,2))    ! output: ratio of max and min percolation rates [day-1]
  ! local 
  character(len=strLen)   :: message                           
                       
  message="zperc/"
  select case(opt)
    case(1); 
      where ( twm_in/=dmiss .and. fsm_in/=dmiss ) 
        zperc=(twm_in+fsm_in*(1.0_dp-zsk_in)+fsm_in*(1.0_dp-zpk_in))/ &
              (fsm_in*zsk_in+fpm_in*zpk_in)
      else where
        zperc=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  return
end function

! *********************************************************************
! REXP parameter 
! *********************************************************************
function rexp( wp_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in)  :: wp_in(:,:)      ! input: wilting point [frac]
  real(dp),    intent(in)  :: gammaPar(:)     ! input: gamma parameter array 
  integer(i2b),intent(in)  :: opt             ! input: option for transfer function form
  ! local 
  real(dp)                 :: rexp(size(wp_in,1),size(wp_in,2)) ! output: tension water maximum [mm]
  character(len=strLen)    :: message                           
  
  message="rexp/"
  associate(g1=>gammaPar(ixGamma%rexp1gamma1))
  select case(opt)
    case(1); 
      where ( wp_in/=dmiss ) 
        rexp=sqrt(wp_in/(g1-0.001))
      else where
        rexp=dmiss
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for saturated hydraulic conductivity (ks)
! *********************************************************************
function ks( sand_in, clay_in, gammaPar, opt)
  implicit none
  ! input
  real(dp), intent(in)  :: sand_in(:,:)                        ! input: sand [percent] 
  real(dp), intent(in)  :: clay_in(:,:)                        ! input: clay [percent] 
  real(dp), intent(in)  :: gammaPar(:)                         ! input: gamma parameter array 
  integer(i2b)          :: opt                                 ! input: option for transfer function form
  ! local 
  real(dp)              :: ks(size(sand_in,1),size(sand_in,2)) ! output: mm/s
  character(len=strLen) :: message                           

  message="ks/"
  ! opt 1: Cosby et al. WRR 1984
  ! opt 2: campbell & shiozawa 1994 
  associate(g1=>gammaPar(ixGamma%ks1gamma1), &
            g2=>gammaPar(ixGamma%ks1gamma2), &
            g3=>gammaPar(ixGamma%ks1gamma3), &
            g4=>gammaPar(ixGamma%ks2gamma1), &
            g5=>gammaPar(ixGamma%ks2gamma2), &
            g6=>gammaPar(ixGamma%ks2gamma3) )
  select case(opt)
    case(1); 
      where ( sand_in /= dmiss .or. clay_in /= dmiss ) 
        ks = g1 + g2*sand_in + g3*clay_in
        ks = (10**ks)*25.4/60/60   ! 25.4 mm/inch. Cosby give Ksat in inch/hr 
      else where
        ks = dmiss 
      end where
    case(2); 
      where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
        ks = g4*exp(g5*sand_in+g6*clay_in)
      else where
        ks = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized'; stop
  end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for bulk density 
! *********************************************************************
function bd( bd_in, gammaPar,opt )
  implicit none
  ! input
  real(dp), intent(in)  :: bd_in(:,:)    ! input: bd from dataset [kg/m3]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  integer(i2b)          :: opt            ! input: option for transfer function form
  ! output 
  ! local 
  real(dp)              :: bd(size(bd_in,1),size(bd_in,2))
  real(dp),parameter    :: bd_min=805.0_dp
  real(dp),parameter    :: bd_max=1880.0_dp
  real(dp),allocatable  :: bdslope(:,:)
  real(dp),allocatable  :: bd_temp(:,:)
  integer(i4b)          :: n1                  ! number of 1st dimension 
  integer(i4b)          :: n2                  ! number of 2nd dimension 
  character(len=strLen) :: message                           

  message="bd/"
  n1=size(bd_in,1)
  n2=size(bd_in,2)
  allocate(bdslope(n1,n2))
  allocate(bd_temp(n1,n2))
  bdslope=0.0_dp
  bd_temp=0.0_dp
  associate(g1=>gammaPar(ixGamma%bd1gamma1))
  select case(opt)
    case(1);  ! 
      where ( bd_in /= dmiss ) 
        bd_temp = g1*bd_in
        bdslope=(bd_temp-bd_min)/(bd_max-bd_min)
        where ( bdslope > 1.0_dp) bdslope=1.0_dp
        where ( bdslope < 0.0_dp) bdslope=0.0_dp
        bd = bdslope*(bd_max-bd_min)+bd_min
      else where
        bd = dmiss 
      end where
    case default;print*,trim(message)//'opt not recognized';stop
  end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for porosity 
! *********************************************************************
function phi(sand_in, clay_in, db_in, gammaPar, opt)      
  implicit none
  ! input
  real(dp), intent(in)       :: sand_in(:,:)   ! input: sand [percent] 
  real(dp), intent(in)       :: clay_in(:,:)   ! input: clay [percent]
  real(dp), intent(in)       :: db_in(:,:)     ! input: bulk density [kg/m^3]
  real(dp), intent(in)       :: gammaPar(:)    ! input: gamma parameter array 
  integer(i2b)               :: opt            ! option for transfer function form
  ! local 
  real(dp)                   :: phi(size(db_in,1),size(db_in,2))  ! estimated porosity [fraction]
  character(len=strLen)      :: message        ! error message

  ! opt 1: Cosby et al. WRR 1984
  ! opt 2: Zacharias & Wessolek 2007
  message="comp_porosity/"

  associate(g1=>gammaPar(ixGamma%phi1gamma1), &
            g2=>gammaPar(ixGamma%phi1gamma2), &
            g3=>gammaPar(ixGamma%phi1gamma3), &
            g4=>gammaPar(ixGamma%phi2gamma1), &
            g5=>gammaPar(ixGamma%phi2gamma2), &
            g6=>gammaPar(ixGamma%phi2gamma3), &
            g7=>gammaPar(ixGamma%phi2gamma4), &
            g8=>gammaPar(ixGamma%phi2gamma5), &
            g9=>gammaPar(ixGamma%phi2gamma6))
    select case(opt)
      case(1);  ! Cosby
        where ( sand_in /= dmiss .or. clay_in /= dmiss ) 
          phi = ( g1+g2*sand_in+g3*clay_in )/100.0_dp
        else where
          phi = dmiss 
        end where
      case(2);  ! Zacharias & Wessolek 2007  
        where ( sand_in /= dmiss .and. clay_in /= dmiss .and. Db_in /= dmiss ) 
          where ( sand_in < 66.5_dp) 
            phi = g4+g5*clay_in+g6*Db_in/1000._dp
          else where
            phi = g7+g8*clay_in+g9*Db_in/1000._dp
          end where 
        else where
          phi = dmiss 
        end where
      case default;print*,trim(message)//'opt not recognized';stop
    end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for field capacity 
! *********************************************************************
function fc(sand_in, phi_in, psis_in, b_in, gammaPar, opt)
  implicit none
  ! input
  real(dp), intent(in)       :: sand_in(:,:)   ! input: sand  
  real(dp), intent(in)       :: phi_in(:,:)    ! input: porosity [fraction]  
  real(dp), intent(in)       :: psis_in(:,:)   ! input: saturation matric potential [kPa] 
  real(dp), intent(in)       :: b_in(:,:)      ! input: slope of cambell retention curve 
  real(dp), intent(in)       :: gammaPar(:)    ! input: gamma parameter array 
  integer(i2b)               :: opt            ! id for transfer function form
  ! output 
  character(len=strLen)      :: message        ! error message
  ! local 
  real(dp)                   :: fc(size(b_in,1),size(b_in,2))        ! estimated field capacity [fraction]
  real(dp),allocatable       :: psi_fc(:,:)    ! matric potential at field capacity [kPa]  
  integer(i4b)               :: nSpoly         ! number of soil polygon 
  integer(i4b)               :: nSLyr          ! number of soil layer 

  ! opt 1: Campbell 1974 
  message="fc/"

  nSpoly=size(phi_in,1)
  nSLyr=size(phi_in,2)
  allocate(psi_fc(nSpoly,nSLyr))
  psi_fc(:,:)=-20
  associate(g1=>gammaPar(ixGamma%fc1gamma1))
  where (sand_in > 69) psi_fc=-10
  select case(opt)
    case(1);  !campbell
      where ( phi_in /= dmiss .and. sand_in /= dmiss .and. b_in /= dmiss .and. psis_in /= dmiss ) 
        fc = g1*phi_in*(psi_fc/psis_in)**(-1/b_in)
      else where
        fc = dmiss 
      end where
    case default; message=trim(message)//'opt not recognized' 
  end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for wilting point 
! *********************************************************************
function wp( phi_in, psis_in, b_in, gammaPar, opt) 
  implicit none
  ! input
  real(dp), intent(in)       :: phi_in(:,:)    ! input: porosity [fraction]  
  real(dp), intent(in)       :: psis_in(:,:)   ! input: saturation matric potential [kPa]  
  real(dp), intent(in)       :: b_in(:,:)      ! input: slope of cambell retention curve [-]
  real(dp), intent(in)       :: gammaPar(:)    ! input: gamma parameter array 
  integer(i2b)               :: opt            ! input: option for transfer function form
  ! local 
  real(dp)                   :: wp(size(b_in,1),size(b_in,2))        ! estimated field capacity [frac]
  real(dp),allocatable       :: psi_wp(:,:)    ! matric potential at wilting point [kPa]  
  integer(i4b)               :: nSpoly         ! number of soil polygon 
  integer(i4b)               :: nSLyr          ! number of soil layer 
  character(len=strLen)      :: message        ! error message

  ! opt 1: Campbell 1974 
  message="wp/"

  nSpoly=size(phi_in,1)
  nSLyr=size(phi_in,2)
  allocate(psi_wp(nSpoly,nSLyr))
  psi_wp(:,:)=-1500
  associate(g1=>gammaPar(ixGamma%wp1gamma1))
  select case(opt)
    case(1);  !Cosby et al. 
      where ( phi_in /= dmiss .and. b_in /= dmiss .and. psis_in /= dmiss ) 
        wp = g1*phi_in*(psi_wp/psis_in)**(-1/b_in)
      else where
        wp = dmiss 
      end where
    case default; print*,trim(message)//'opt not recognized';stop
  end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for b (slope of retention curve in log space)
! *********************************************************************
function ret_curve(sand_in, clay_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in)   :: sand_in(:,:)   ! input: sand [percent]  
  real(dp),    intent(in)   :: clay_in(:,:)   ! input: clay [percent]
  real(dp),    intent(in)   :: gammaPar(:)    ! input: gamma parameter array 
  integer(i2b)              :: opt            ! input: option for transfer function form
  ! local 
  real(dp)                  :: ret_curve(size(sand_in,1),size(sand_in,2)) ! computed [-] 
  character(len=strLen)     :: message        ! error message

  ! opt 1: Cosby et al. WRR 1984
  message="ret_curve/"
  associate(g1=>gammaPar(ixGamma%b1gamma1), &
            g2=>gammaPar(ixGamma%b1gamma2), &
            g3=>gammaPar(ixGamma%b1gamma2))
    select case(opt)
      case(1); 
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          ret_curve = g1+g2*sand_in+g3*clay_in
        else where
          ret_curve = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized';stop
    end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for saturation matric potential 
! *********************************************************************
function psis( sand_in, silt_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),     intent(in)   :: sand_in(:,:)   ! input: sand percentage [percentage]
  real(dp),     intent(in)   :: silt_in(:,:)   ! input: silt percentage [percentage]
  real(dp),     intent(in)   :: gammaPar(:)    ! input: gamma parameter array 
  integer(i2b), intent(in)   :: opt            ! input: option for transfer function form
  ! local 
  real(dp)                   :: psis(size(sand_in,1),size(sand_in,2))      ! output: saturation matric potential [kPa]  
  character(len=strLen)      :: message        ! error message

  ! opt 1: Cosby et al. WRR 1984
  message="psis/"
  associate(g1=>gammaPar(ixGamma%psis1gamma1), &
            g2=>gammaPar(ixGamma%psis1gamma2), &
            g3=>gammaPar(ixGamma%psis1gamma2))
    select case(opt)
      case(1);  !Cosby et al. 
        where ( sand_in /= dmiss .and. silt_in /= dmiss ) 
          psis = g1 + g2*sand_in + g3*silt_in
          psis = -1.0_dp*10**psis*0.0980665_dp        ! 0.0980665 kPa/cm-H2O. Cosby give psi_sat in cm of water (cm-H2O)
        else where
          psis = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized';stop
    end select
  end associate
  return
end function

! *********************************************************************
! pedo-transfer function for specific yield  
! *********************************************************************
function myu(phi_in, fc_in, gammaPar, opt)
  implicit none
  ! input
  real(dp),    intent(in)    :: phi_in(:,:)     ! input: porosity [frac]
  real(dp),    intent(in)    :: fc_in(:,:)      ! input: field capacity [frac]
  real(dp),    intent(in)    :: gammaPar(:)     ! input: gamma parameter array 
  integer(i2b),intent(in)    :: opt             ! input: option for transfer function form
  ! local 
  real(dp)                   :: myu(size(fc_in,1),size(fc_in,2))        ! output: specific yield [-]  
  character(len=strLen)      :: message        ! error message

  ! opt 1: Koren et al. 2003
  message="myu/"
  associate(g1=>gammaPar(ixGamma%myu1gamma1), &
            g2=>gammaPar(ixGamma%myu1gamma2))
  select case(opt)
    case(1);  ! koren
      where ( phi_in /= dmiss .and. fc_in /= dmiss ) 
        myu = g1*(phi_in-fc_in)**g2 
      else where
        myu = dmiss 
      end where
    case default;print*,trim(message)//'opt not recognized';stop
  end select
  end associate
  return
end function

! *********************************************************************
! Monthly LAI 
! *********************************************************************
function lai( lai_in, gammaPar, opt )
  implicit none
  ! input
  real(dp), intent(in)  :: lai_in(:,:)                        ! input: monthly LAI 
  real(dp), intent(in)  :: gammaPar(:)                        ! input: gamma parameter array 
  integer(i2b)          :: opt                                ! input: option for transfer function form
  ! output
  real(dp)              :: lai(size(lai_in,1),size(lai_in,2)) ! output: adjusted LAI
  ! local 
  real(dp),parameter    :: lai_min=0.0_dp                     ! minimum plausible LAI (0 m2/m2 for bare ground)
  real(dp),parameter    :: lai_max=10.0_dp                    ! maximum plausible LAI (10 m2/m2 for dense conifer forest)
  real(dp),allocatable  :: lai_temp(:,:)
  real(dp),allocatable  :: laislope(:,:)
  integer(i4b)          :: n1                                 ! number of 1st dimension 
  integer(i4b)          :: n2                                 ! number of 2nd dimension 
  character(len=strLen) :: message                           

  message="lai/"
  n1=size(lai_in,1)
  n2=size(lai_in,2)
  allocate(laislope(n1,n2))
  allocate(lai_temp(n1,n2))
  laislope=0.0_dp
  lai_temp=0.0_dp
  associate(g1=>gammaPar(ixGamma%lai1gamma1))
  select case(opt)
    case(1);  ! 
      where ( lai_in /= dmiss ) 
        lai_temp = g1*lai_in
        laislope=(lai_temp-lai_min)/(lai_max-lai_min)
        where ( laislope > 1.0_dp) laislope=1.0_dp
        where ( laislope < 0.0_dp) laislope=0.0_dp
        lai = laislope*(lai_max-lai_min)+lai_min
      else where
        lai = dmiss 
      end where
    case default;print*,trim(message)//'opt not recognized';stop
  end select
  end associate
  return
end function 

end module tf 
