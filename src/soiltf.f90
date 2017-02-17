module soiltf
! Compute model soil parameter transfer function 
use nrtype                                        ! variable types, etc.
use data_type                                     ! Including custum data structure definition
use public_var                                     ! Including common constant (physical constant, other e.g., missingVal, etc.)
use var_lookup,   only:ixVarSoilData, ixPar, nPar

implicit none

private

public::comp_soil_model_param

contains

! ********************************************************************************************
! Public subroutine: Execute computation of model-dependant soil parameter transfer function  
! *********************************************************************************************
subroutine comp_soil_model_param(parSxySz,          &  ! in/output: soil parameter values for all soil polygons that are input
                                 sdata,             &  ! input: soil data
                                 gammaParMasterMeta,&  ! input: gamma parameter meta file - val of calibrating parammeter is adjusted via calibration 
                                 nSLyr,             &  ! input: number of soil layers
                                 nSPoly,            &  ! input: number of soil polygons
                                 ierr,message) 

  use globalData, only:betaInGamma, betaNeeded
  use get_ixname, only:get_ixPar
  implicit none
  ! in/out
  type(namedvar2),      intent(inout) :: parSxySz(:)            ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  ! input
  type(namevar),        intent(in)    :: sdata(:)               ! storage of soil data strucuture
  type(par_meta)                      :: gammaParMasterMeta(:)
  integer(i4b),         intent(in)    :: nSLyr                  ! number of soil layer
  integer(i4b),         intent(in)    :: nSPoly                 ! number of soil polygons
  ! output
  integer(i4b),         intent(out)   :: ierr                   ! error code
  character(len=strLen),intent(out)   :: message                ! error message for current routine
  ! Local 
  type(namedvar2)                     :: ParTemp(nPar)          ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  integer(i4b)                        :: ix                     ! index of gamma parameter 
  integer(i4b)                        :: idBeta                 ! id of beta parameter array 
  integer(i4b)                        :: iParm                  ! Loop index of model parameters (e.g., VIC)
  logical(lgc)                        :: checkDone(nPar)        ! used to check if the VIC parameter is processed

  ierr=0; message="comp_soil_model_param/"
  first: associate(sclass  => sdata(ixVarSoilData%soilclass)%ivar2,   &
                   hslyrs  => sdata(ixVarSoilData%hslyrs)%dvar2,      &
                   sand    => sdata(ixVarSoilData%sand_frc)%dvar2,    & 
                   silt    => sdata(ixVarSoilData%silt_frc)%dvar2,    & 
                   clay    => sdata(ixVarSoilData%clay_frc)%dvar2,    & 
                   bulkd   => sdata(ixVarSoilData%bulk_density)%dvar2,& 
                   elestd  => sdata(ixVarSoilData%ele_std)%dvar1,     & 
                   elemean => sdata(ixVarSoilData%ele_mean)%dvar1,    & 
                   slpmean => sdata(ixVarSoilData%slp_mean)%dvar1,    &
                   gammaPar=> gammaParMasterMeta(:)%val)
  do iParm = 1,size(betaNeeded)
    ix = get_ixPar(betaNeeded(iParm)) 
    allocate(ParTemp(ix)%varData(nSLyr,nSPoly) ,stat=ierr); if(ierr/=0)then;message=trim(message)//'error allocating ParTemp';stop;endif
    associate (xPar => ParTemp(ix)%varData )
    select case(ix)
      case(ixPar%ks)
        checkDone(ix)=.true.
        xPar = ks( sand, clay, gammaPar, 1_i2b)
      case(ixPar%bd)
        checkDone(ix)=.true.
        xPar  = bd( bulkd, gammaPar )
      case(ixPar%phi)
        if(.not.checkDone(ixPar%bd)) then;ierr=10;message=trim(message)//'need to process bd before phi';return;endif
        checkDone(ix)=.true.
        xPar = phi( sand, clay, ParTemp(ixPar%bd)%varData,gammaPar, 1_i2b)
      case(ixPar%b)
        checkDone(ix)=.true.
        xPar = ret_curve( sand, clay, gammaPar, 1_i2b)
      case(ixPar%psis)
        checkDone(ix)=.true.
        xPar = psis( sand, silt, gammaPar, 1_i2b)
      case(ixPar%fc)
        if(.not.checkDone(ixPar%psis)) then;ierr=10;message=trim(message)//'need to process psis before fc';return;endif
        if(.not.checkDone(ixPar%phi))  then;ierr=10;message=trim(message)//'need to process phi before fc';return;endif
        if(.not.checkDone(ixPar%b))    then;ierr=10;message=trim(message)//'need to process b before fc';return;endif
        checkDone(ix)=.true.
        xPar = fc(sand, ParTemp(ixPar%phi)%varData, ParTemp(ixPar%psis)%varData, ParTemp(ixPar%b)%varData,gammaPar, 1_i2b)
      case(ixPar%wp)
        if(.not.checkDone(ixPar%psis)) then;ierr=10;message=trim(message)//'need to process psis before wp';return;endif
        if(.not.checkDone(ixPar%phi))  then;ierr=10;message=trim(message)//'need to process phi before wp';return;endif
        if(.not.checkDone(ixPar%b))    then;ierr=10;message=trim(message)//'need to process b before wp';return;endif
        checkDone(ix)=.true.
        xPar = wp( ParTemp(ixPar%phi)%varData, ParTemp(ixPar%psis)%varData, ParTemp(ixPar%b)%varData, gammaPar, 1_i2b) 
      case(ixPar%myu)
        if(.not.checkDone(ixPar%phi))  then;ierr=10;message=trim(message)//'need to process phi before myu';return;endif
        if(.not.checkDone(ixPar%fc))   then;ierr=10;message=trim(message)//'need to process fc before myu';return;endif
        checkDone(ix)=.true.
        xPar= myu( ParTemp(ixPar%phi)%varData, ParTemp(ixPar%fc)%varData, gammaPar, 1_i2b)
      case(ixPar%binfilt)
        checkDone(ix)=.true.
        xPar=spread( infilt( elestd, gammaPar ), 1, nSLyr)
      case(ixPar%D1)
        if(.not.checkDone(ixPar%ks))  then;ierr=10;message=trim(message)//'need to process "ks" before "D1"';return;endif
        if(.not.checkDone(ixPar%phi)) then;ierr=10;message=trim(message)//'need to process "phi" before "D1"';return;endif
        checkDone(ix)=.true. 
        xPar=D1( slpmean,                    &
                 ParTemp(ixPar%ks)%varData,  &
                 ParTemp(ixPar%phi)%varData, & 
                 hslyrs,                     &
                 gammaPar)
      case(ixPar%Ds)
        if(.not.checkDone(ixPar%D1))    then;ierr=10;message=trim(message)//'need to process "D1" before "Ds"';return;endif
        if(.not.checkDone(ixPar%D3))    then;ierr=10;message=trim(message)//'need to process "D3" before "Ds"';return;endif
        if(.not.checkDone(ixPar%Dsmax)) then;ierr=10;message=trim(message)//'need to process "Dsmax" before "Ds"';return;endif
        checkDone(ix)=.true. 
        xPar=Ds( ParTemp(ixPar%D1)%varData, ParTemp(ixPar%D3)%varData, ParTemp(ixPar%Dsmax)%varData )
      case(ixPar%D4)
        checkDone(ix)=.true.
        xPar=D4(gammaPar)
      case(ixPar%c)
        if(.not.checkDone(ixPar%D4)) then;ierr=10;message=trim(message)//'need to process "D4" before "c"';return;endif
        checkDone(ix)=.true.
        xPar=cexpt(ParTemp(ixPar%D4)%varData)
      case(ixPar%SD)
        checkDone(ix)=.true.
        xPar=soilDensity(ParTemp(ixPar%sd)%varData, gammaPar)
      case(ixPar%expt)
        checkDone(ix)=.true.
        xPar=expt( ParTemp(ixPar%b)%varData, gammaPar )
      case(ixPar%D2)
        if(.not.checkDone(ixPar%ks)) then;ierr=10;message=trim(message)//'need to process "ksat" before "D2"';return;endif
        if(.not.checkDone(ixPar%D4)) then;ierr=10;message=trim(message)//'need to process "D4" before "D2"';return;endif
        checkDone(ix)=.true.
        xPar=D2( slpmean, ParTemp(ixPar%ks)%varData, ParTemp(ixPar%D4)%varData, gammaPar )
      case(ixPar%Dsmax)
        if(.not.checkDone(ixPar%D1)) then;ierr=10;message=trim(message)//'need to process "D1" before "Dsmax"';return;endif
        if(.not.checkDone(ixPar%D2)) then;ierr=10;message=trim(message)//'need to process "D2" before "Dsmax"';return;endif
        if(.not.checkDone(ixPar%D3)) then;ierr=10;message=trim(message)//'need to process "D3" before "Dsmax"';return;endif
        if(.not.checkDone(ixPar%c))  then;ierr=10;message=trim(message)//'need to process "c" before "Dsmax"';return;endif
        checkDone(ix)=.true. 
        xPar=Dsmax( ParTemp(ixPar%D1)%varData,   &
                    ParTemp(ixPar%D2)%varData,   & 
                    ParTemp(ixPar%D3)%varData,   & 
                    ParTemp(ixPar%c)%varData,    & 
                    ParTemp(ixPar%phi)%varData,  &
                    hslyrs )
      case(ixPar%bbl)
        if(.not.checkDone(ixPar%expt)) then;ierr=10;message=trim(message)//'need to process "expt" before "bubble"';return;endif
        checkDone(ix)=.true.
        xPar=bubble( ParTemp(ixPar%expt)%varData, gammaPar)
      case(ixPar%WcrFrac)
        checkDone(ix)=.true.
        xPar=WcrFrac( ParTemp(ixPar%fc)%varData,ParTemp(ixPar%phi)%varData,gammaPar )
      case(ixPar%WpwpFrac)
        checkDone(ix)=.true.
        xPar=WpwpFrac( ParTemp(ixPar%wp)%varData, ParTemp(ixPar%phi)%varData, gammaPar)
      case(ixPar%D3)
        if(.not.checkDone(ixPar%fc))   then;ierr=10;message=trim(message)//'need to process fc before D3';return;endif
        checkDone(ix)=.true.
        xPar=D3( ParTemp(ixPar%fc)%varData, hslyrs, gammaPar )
      case(ixPar%Ws)
        if(.not.checkDone(ixPar%D3))then;ierr=10;message=trim(message)//'need to process "D3" before "Dsmax"';return;endif 
        checkDone(ix)=.true. 
        xPar=Ws( ParTemp(ixPar%D3)%varData, ParTemp(ixPar%phi)%varData,hslyrs )
      case(ixPar%twm)
        checkDone(ix)=.true. 
        xPar= twm( ParTemp(ixPar%fc)%varData,ParTemp(ixPar%wp)%varData,hslyrs )
      case(ixPar%fwm)
        checkDone(ix)=.true. 
        xPar= fwm( ParTemp(ixPar%phi)%varData,ParTemp(ixPar%fc)%varData,hslyrs )
      case(ixPar%fsm)
        if(.not.checkDone(ixPar%fwm))then;ierr=10;message=trim(message)//'need to process "fwm" before "fsm"';return;endif 
        checkDone(ix)=.true. 
        xPar= fsm( ParTemp(ixPar%fwm)%varData,ParTemp(ixPar%phi)%varData,ParTemp(ixPar%wp)%varData,gammaPar )
      case(ixPar%fpm)
        if(.not.checkDone(ixPar%fwm))then;ierr=10;message=trim(message)//'need to process "fwm" before "fpm"';return;endif 
        if(.not.checkDone(ixPar%fsm))then;ierr=10;message=trim(message)//'need to process "fsm" before "fpm"';return;endif 
        checkDone(ix)=.true. 
        xPar= fpm( ParTemp(ixPar%fwm)%varData, ParTemp(ixPar%fsm)%varData )
      case(ixPar%zk)
        checkDone(ix)=.true. 
        xPar= zk( ParTemp(ixPar%phi)%varData,ParTemp(ixPar%fc)%varData,gammaPar )
      case(ixPar%zsk)
        checkDone(ix)=.true. 
        xPar= zsk( ParTemp(ixPar%phi)%varData,ParTemp(ixPar%fc)%varData,ParTemp(ixPar%wp)%varData,gammaPar )
      case(ixPar%zpk)
        checkDone(ix)=.true. 
        xPar= zpk( ParTemp(ixPar%ks)%varData,ParTemp(ixPar%myu)%varData,hslyrs,gammaPar )
      case(ixPar%pfree)
        checkDone(ix)=.true. 
        xPar= pfree( ParTemp(ixPar%phi)%varData,ParTemp(ixPar%wp)%varData,gammaPar )
      case(ixPar%zperc)
        if(.not.checkDone(ixPar%twm))then;ierr=10;message=trim(message)//'need to process "twm" before "pfree"';return;endif 
        if(.not.checkDone(ixPar%fsm))then;ierr=10;message=trim(message)//'need to process "fsm" before "pfree"';return;endif 
        if(.not.checkDone(ixPar%zsk))then;ierr=10;message=trim(message)//'need to process "zsk" before "pfree"';return;endif 
        if(.not.checkDone(ixPar%fpm))then;ierr=10;message=trim(message)//'need to process "fpm" before "pfree"';return;endif 
        if(.not.checkDone(ixPar%zpk))then;ierr=10;message=trim(message)//'need to process "zpk" before "pfree"';return;endif 
        checkDone(ix)=.true. 
        xPar= zperc(ParTemp(ixPar%twm)%varData,ParTemp(ixPar%fsm)%varData,ParTemp(ixPar%zsk)%varData,ParTemp(ixPar%fpm)%varData,ParTemp(ixPar%zsk)%varData)
      case(ixPar%rexp)
        checkDone(ix)=.true. 
        xPar= rexp( ParTemp(ixPar%wp)%varData,gammaPar )
    end select ! end of parameter case
    end associate
  end do ! end of parameter loop
  end associate
  ! extract beta parameters in 'CalPar' list
  do iParm=1,size(betaInGamma)
    idBeta=get_ixPar(trim(betaInGamma(iParm)))
    parSxySz(iParm)%varData=parTemp(idBeta)%varData 
  enddo
  return
end subroutine 

! *********************************************************************
! Library of hydrologic model transfer functions
! *********************************************************************
! infilt parameter 
! *********************************************************************
function infilt(elestd_in, gammaPar)
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
  real(dp), intent(in)  :: elestd_in(:)
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: infilt(size(elestd_in)) 
  real(dp),parameter    :: infilt_min=0.03_dp
  real(dp),parameter    :: infilt_max=0.50_dp
  
  associate(g1=>gammaPar(ixPar%binfilt1gamma1), &
            g2=>gammaPar(ixPar%binfilt1gamma2))
  where ( elestd_in /= dmiss ) 
    infilt = (log(elestd_in+verySmall)-g1)/(log(elestd_in+verySmall)+g2*10) !transfer function 
  else where
    infilt = dmiss 
  end where
   !cap value with upper and lower bounds 
  where ( infilt > infilt_max ) infilt=infilt_max 
  where ( infilt /= dmiss .and. infilt < infilt_min ) infilt=infilt_min 
  end associate
  return
end function

! *********************************************************************
! residual_moist parameter 
! *********************************************************************
function residMoist()
 implicit none
 ! input
 ! output 
 ! local 
 real(dp) :: residMoist
 
 residMoist = 0._dp
 return
end function
         
! ***********
!  Nijssen basiflow D1 parameter
! *********************************************************************
function D1(slope_in, ks_in, phi_in, h_in, gammaPar)
  implicit none
  ! input
  real(dp), intent(in)  :: slope_in(:)   ! mean slope [%]
  real(dp), intent(in)  :: ks_in(:,:)    ! Ksat [mm/s]
  real(dp), intent(in)  :: phi_in(:,:)   ! porosity [-]
  real(dp), intent(in)  :: h_in(:,:)     ! layer thickness [m]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
! output 
! local 
  real(dp)              :: D1(size(ks_in,1),size(ks_in,2)) !output: [day-1]
  real(dp),allocatable  :: slope2d(:,:)
  integer(i4b)          :: n1            ! number of 1st dimension 
  integer(i4b)          :: n2            ! number of 1st dimension 
  real(dp), allocatable :: S(:,:)              ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp), parameter   :: D1_min=0.0001_dp
  real(dp), parameter   :: D1_max=1.0_dp
  !integer(i4b)          :: err                       ! error code
   
  ! local variable allocation
  n1=size(D1,1)
  n2=size(D1,2)
  allocate(S(n1,n2)) 
  allocate(slope2d(n1,n2)) 
  S=phi_in*h_in
  S=1.0_dp
  slope2d=spread(slope_in,1,n1)
  associate(g1=>gammaPar(ixPar%D11gamma1))
  ! compute parameters 
  where ( slope2d /= dmiss .and. Ks_in /= dmiss )
    D1 = S**(-1)*10**(-1*g1)*Ks_in*(60*60*24)*(slope2d*0.01)
  else where
    D1 = dmiss
  end where
  ! cap value with upper and lower bounds 
  where ( D1 > D1_max ) D1=D1_max 
  where ( D1 > 0._dp .and. D1 < D1_min ) D1=D1_min 
  end associate 
  return
end function

! ***********
!  Arno basiflow Ds parameter
! *********************************************************************
function Ds( D1, D3, Dsmax)
 implicit none
! input
 real(dp), intent(in)  :: D1(:,:)       ! nijssen baseflow D1 parameter [day-1]
 real(dp), intent(in)  :: D3(:,:)       ! nijssen baseflow D3 parameter [mm]
 real(dp), intent(in)  :: Dsmax(:,:)    ! ARNO Dsmax parameter [mm/day]
! output 
! local 
 real(dp)              :: Ds(size(D1,1),size(D1,2)) !output: [-]
 real(dp), parameter   :: Ds_min=0.0001_dp
 real(dp), parameter   :: Ds_max=1.0_dp

 where ( D1 /= dmiss .and. D3 /= dmiss .and. Dsmax /= dmiss )
   Ds = D1 * D3 / Dsmax
 else where
   Ds = dmiss
 end where
 ! cap value with upper and lower bounds 
 where ( Ds > Ds_max ) Ds=Ds_max 
 where ( Ds > 0._dp .and. Ds < Ds_min ) Ds=Ds_min 
 return
end function

! ***********
! Nijssen baseflow D2 parameter
! *********************************************************************
function D2(slope_in, ks_in, D4_in, gammaPar) 
  implicit none
  ! input
  real(dp), intent(in)  :: slope_in(:)   ! slope percent
  real(dp), intent(in)  :: Ks_in(:,:)    ! ksat [mm/s]
  real(dp), intent(in)  :: D4_in(:,:)    ! VIC D4 parameter [-]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: D2(size(Ks_in,1),size(Ks_in,2)) ! output [day^-D4]
  real(dp),allocatable  :: slope2d(:,:)
  integer(i4b)          :: n1           ! number of element for 1st dimension
  integer(i4b)          :: n2           ! number of element for 2nd dimension
  real(dp), allocatable :: S(:,:)             ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp), parameter   :: D2_min=0.0001_dp
  real(dp), parameter   :: D2_max=1.0_dp
  !integer(i4b)          :: err                       ! error code
  
 ! local variable allocation
  n1=size(D2,1)
  n2=size(D2,2)
  allocate(S(n1,n2))
  allocate(slope2d(n1,n2)) 
  S=1.0_dp
  slope2d=spread(slope_in,1,n1)
  associate(g1=>gammaPar(ixPar%D21gamma1))
 ! compute parameters 
  where ( slope2d /= dmiss .and. Ks_in /= dmiss )
    D2 = S**(-1*D4_in)*10**(-1*g1)*Ks_in*(60*60*24)*(slope2d*0.01)
  else where
    D2 = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( D2 > D2_max ) D2=D2_max
  where ( D2 > 0._dp .and. D2 < D2_min ) D2=D2_min
  end associate
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
                h_in)          ! input:  Soil layer thickness [m]
  implicit none
  ! input
  real(dp), intent(in)  :: D1(:,:)        ! Nijssen baseflow D1 parameter [day^-1]
  real(dp), intent(in)  :: D2(:,:)        ! Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
  real(dp), intent(in)  :: D3(:,:)        ! Nijssen baseflow D3 parameter [mm]
  real(dp), intent(in)  :: c(:,:)         ! c parameter [mm]
  real(dp), intent(in)  :: phi_in(:,:)    ! porosity [cm^3/cm^-3]
  real(dp), intent(in)  :: h_in(:,:)      ! Soil layer thickness [m]
  ! output 
  ! local 
  real(dp)              :: Dsmax(size(D1,1),size(D1,2))     ! Dsmax parameter for Arno baseflow [mm day-1]
  real(dp), parameter   :: Dsmax_min=0.1_dp
  real(dp), parameter   :: Dsmax_max=30.0_dp

  where ( phi_in /= dmiss .and. h_in /= dmiss )
    Dsmax = D2*(phi_in*h_in*1000-D3)**c+D1*(phi_in*h_in*1000)
  else where
    Dsmax = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( Dsmax > Dsmax_max ) Dsmax=Dsmax_max
  where ( Dsmax > 0._dp .and. Dsmax < Dsmax_min ) Dsmax=Dsmax_min
  return
end function

! ***********
!  Nijssen baseflow D3 parameter
! *********************************************************************
function D3( fc_in, h_in, gammaPar ) 
  implicit none
  ! input
  real(dp), intent(in)  :: fc_in(:,:)    ! input: field capacity [frac]
  real(dp), intent(in)  :: h_in(:,:)     ! input: layer thickness [m]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: D3(size(fc_in,1),size(fc_in,2)) !output: [mm]
  real(dp), parameter   :: D3_min=0.0001_dp
  real(dp), parameter   :: D3_max=1000.0_dp

  associate(g1=>gammaPar(ixPar%D31gamma1))
  where ( fc_in /= dmiss .and. h_in /= dmiss ) 
    D3 = g1* fc_in * (h_in*1000)
  else where
    D3 = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( D3 > D3_max ) D3=D3_max
  where ( D3 > 0._dp .and. D3 < D3_min ) D3=D3_min 
  end associate
  return
end function

! ***********
!  Arno baseflow Ws parameter (conversion equation)
! *********************************************************************
function Ws( D3,         & ! input:  D3 parameter [mm]
             phi_in,     & ! input:  porosity [cm^3/cm^-3]
             h_in)          ! input:  Soil layer thickness [m]
  implicit none
  ! input
  real(dp), intent(in)  :: D3(:,:)
  real(dp), intent(in)  :: phi_in(:,:)
  real(dp), intent(in)  :: h_in(:,:)
  ! output 
  ! local 
  real(dp)              :: Ws(size(D3,1),size(D3,2))   ! output [-]
  real(dp), parameter   :: Ws_min=0.05_dp
  real(dp), parameter   :: Ws_max=1.0_dp

  where ( phi_in /= dmiss .and. h_in /= dmiss .and. D3 /= dmiss ) 
    Ws = D3 / phi_in / (h_in*1000)
  else where
    Ws = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( Ws > Ws_max ) Ws=Ws_max
  where ( Ws > 0._dp .and. Ws < Ws_min ) Ws=Ws_min 
  return
end function

! ***********
!  Nijssen baseflow D4 parameter
! *********************************************************************
function D4( gammaPar )
  implicit none
  ! input
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)    :: D4
  
  associate(g1=>gammaPar(ixPar%D41gamma1))
  D4 = g1 
  end associate
  return
end function D4 

! ***********
!  c parameter
! *********************************************************************
function cexpt( D4 )
  implicit none
  ! input
  real(dp), intent(in)  :: D4(:,:)
  ! output 
  ! local 
  real(dp)              :: cexpt(size(D4,1),size(D4,2))
  
  cexpt = D4 
  return
end function cexpt 

! ***********
! computing expt parameter 
! *********************************************************************
function expt( b_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in)   :: b_in(:,:)     ! inpuy [-]
  real(dp), intent(in)   :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)               :: expt(size(b_in,1),size(b_in,2)) ! exponent in campbel equation [-]

  associate(g1=>gammaPar(ixPar%exp1gamma1), &
            g2=>gammaPar(ixPar%exp1gamma2))
  where ( b_in /= dmiss ) 
    expt = g1+g2*b_in
  else where
    expt = dmiss 
  end where
  end associate
  return
end function expt 

! ************
! computing init_moist parameter  
! *********************************************************************
function initMoist( phi_in, h_in)
  implicit none
  ! input
  real(dp), intent(in)  :: phi_in(:,:)       ! porosity [-]
  real(dp), intent(in)  :: h_in(:,:)         ! thickness [m]
  ! output 
  ! local  
  real(dp)              :: initMoist(size(phi_in,1),size(phi_in,2))
  
  where ( phi_in /= dmiss ) 
    initMoist = phi_in*(h_in*1000.0_dp)
  else where
    initMoist = dmiss 
  end where
  return
end function

! ***********
! bubble parameter 
! *********************************************************************
function bubble( expt_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in)  :: expt_in(:,:) 
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: bubble(size(expt_in,1),size(expt_in,2))
  
  associate(g1=>gammaPar(ixPar%bbl1gamma1), &
            g2=>gammaPar(ixPar%bbl1gamma2))
  where ( expt_in /= dmiss ) 
    bubble = g1*expt_in+g2
  else where
    bubble = dmiss 
  end where
  end associate
  return
end function

! ***********
! soil_density parameter 
! *********************************************************************
function soilDensity( srho_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in)  :: srho_in(:,:)  ! input: kg/m^3]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: soilDensity(size(srho_in,1),size(srho_in,2))
  
  associate(g1=>gammaPar(ixPar%sd1gamma1))
  where ( srho_in /= dmiss ) 
    soilDensity = g1*srho_in
  else where
    soilDensity = dmiss 
  end where
  end associate
  return
end function

! ***********
! WcrFrac parameter  
! *********************************************************************
function WcrFrac(fc_in, phi_in, gammaPar)
  implicit none
  ! input
  real(dp), intent(in)  :: fc_in(:,:)    ! input: field capacity [frac]
  real(dp), intent(in)  :: phi_in(:,:)   ! input: porosity [frac]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: WcrFrac(size(fc_in,1),size(fc_in,2)) !output: [frac]
  
  associate(g1=>gammaPar(ixPar%WcrFrac1gamma1))
  where ( fc_in /= dmiss .and. phi_in /= dmiss ) 
    WcrFrac = g1*fc_in/phi_in
  else where
    wcrFrac = dmiss
  end where
  end associate
  return
end function

! ************
! WpwpFrac parameter  
! *********************************************************************
function WpwpFrac( wp_in, phi_in, gammaPar)
  implicit none
  ! input
  real(dp), intent(in)  :: wp_in(:,:)    ! Wilting point [frac]
  real(dp), intent(in)  :: phi_in(:,:)   ! Porosity [frac]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! local 
  real(dp)              :: WpwpFrac(size(wp_in,1),size(wp_in,2)) !output: [frac] 
  
  associate(g1=>gammaPar(ixPar%WpwpFrac1gamma1))
  where ( wp_in /= dmiss .and. phi_in /= dmiss ) 
    WpwpFrac = g1*wp_in/phi_in
  else where
    wpwpFrac = dmiss
  end where
  end associate
  return
end function

! ************
! TWM parameter (UZTWM for top layer and LZTWM for bottom layer)
! *********************************************************************
function twm( fc_in, wp_in, h_in )
  implicit none
  ! input
  real(dp), intent(in) :: fc_in(:,:)                       ! input: field capacity [-]           
  real(dp), intent(in) :: wp_in(:,:)                       ! input: wilting point [-]
  real(dp), intent(in) :: h_in(:,:)                        ! input: thickness of layer [mm]
  ! local 
  real(dp)             :: twm(size(wp_in,1),size(wp_in,2)) ! output: tension water maximum [mm]

  where ( wp_in/=dmiss .and. fc_in/=dmiss ) 
    twm=(fc_in-wp_in)*h_in*1000.0_dp  ! convert m to mm 
  else where
    twm=dmiss
  end where
  return
end function 

! *********************************************************************
! FWM parameter (UZFWM for top layer and LZFWM for bottom layer)
! *********************************************************************
function fwm( phi_in, fc_in, h_in )
  implicit none
  ! input
  real(dp), intent(in) :: phi_in(:,:)                      ! input: porosity [-]
  real(dp), intent(in) :: fc_in(:,:)                       ! input: field capacity [-]           
  real(dp), intent(in) :: h_in(:,:)                        ! input: thickness of layer [mm]
  ! local 
  real(dp)             :: fwm(size(fc_in,1),size(fc_in,2)) ! output: free water maximum [mm]

  where ( phi_in/=dmiss .and. fc_in/=dmiss .and. h_in/=dmiss ) 
    fwm=(phi_in-fc_in)*h_in*1000.0_dp ! convert m to mm
  else where
    fwm=dmiss
  end where
  return
end function 

! *********************************************************************
! FSM parameter (LZFSM for bottom layer)
! *********************************************************************
function fsm( fwm_in, phi_in, wp_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in)  :: fwm_in(:,:)   ! Free water maximum 
  real(dp), intent(in)  :: phi_in(:,:)   ! Porosity
  real(dp), intent(in)  :: wp_in(:,:)    ! Wilting point
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! local 
  real(dp)             :: fsm(size(wp_in,1),size(wp_in,2)) ! output: supplementary free water maximum [mm]

  associate(g1=>gammaPar(ixPar%fsm1gamma1))
  where ( phi_in/=dmiss .and. wp_in/=dmiss .and. fwm_in/=dmiss ) 
    fsm=fwm_in*(wp_in/phi_in)**g1
  else where
    fsm=dmiss
  end where
  end associate
  return
end function 

! *********************************************************************
! FPM parameter (LZFPM for bottom layer)
! *********************************************************************
function fpm( fwm_in, fsm_in )
  implicit none
  ! input
  real(dp), intent(in)  :: fwm_in(:,:)   ! Free water maximum 
  real(dp), intent(in)  :: fsm_in(:,:)   ! supplementary Free water maximum 
  ! local 
  real(dp)             :: fpm(size(fwm_in,1),size(fwm_in,2)) ! output: primary free water maximum [mm]

  where ( fwm_in/=dmiss .and. fsm_in/=dmiss )
    fpm=fwm_in-fsm_in
  else where
    fpm=dmiss
  end where
  return
end function 

! *********************************************************************
! ZK parameter (UZK for top layer)
! *********************************************************************
function zk( phi_in, fc_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in) :: phi_in(:,:)   ! Porosity [frac]
  real(dp), intent(in) :: fc_in(:,:)    ! field capacity [frac]
  real(dp), intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  ! local 
  real(dp)             :: zk(size(fc_in,1),size(fc_in,2)) ! output: draw coefficient from free water content [/day] 

  associate(g1=>gammaPar(ixPar%zk1gamma1))
  where ( phi_in/=dmiss .and. fc_in/=dmiss ) 
    zk=1.0_dp-(fc_in/phi_in)**g1
  else where
    zk=dmiss
  end where
  end associate
  return
end function 

! *********************************************************************
! ZSK parameter (LZSK for top layer)
! *********************************************************************
function zsk( phi_in, fc_in, wp_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in) :: phi_in(:,:)   ! Porosity [frac]
  real(dp), intent(in) :: fc_in(:,:)    ! field capacity [frac]
  real(dp), intent(in) :: wp_in(:,:)    ! Wilting point [frac]
  real(dp), intent(in) :: gammaPar(:)   ! input: gamma parameter array 
  ! local 
  real(dp)             :: zsk(size(fc_in,1),size(fc_in,2)) ! output: draw coefficient from supplementary free water content [/day] 

  associate(g1=>gammaPar(ixPar%zsk1gamma1),&
            g2=>gammaPar(ixPar%zsk1gamma2))
  where ( phi_in/=dmiss .and. fc_in/=dmiss .and. wp_in/=dmiss ) 
    zsk=(1.0_dp-(fc_in/phi_in)**g1)/(1.0_dp+g2*(1.0_dp-wp_in)) 
  else where
    zsk=dmiss
  end where
  end associate
  return
end function 

! *********************************************************************
! ZPK parameter (LZPK for top layer)
! *********************************************************************
function zpk( ks_in, h_in, myu_in, gammaPar )
  implicit none
  ! input
  real(dp), intent(in) :: ks_in(:,:)                       ! Ksat [mm/s]
  real(dp), intent(in) :: myu_in(:,:)                      ! specific yield [-]
  real(dp), intent(in) :: h_in(:,:)                        ! input: thickness of layer [m]
  real(dp), intent(in) :: gammaPar(:)                      ! input: gamma parameter array 
  ! local 
  real(dp)             :: zpk(size(ks_in,1),size(ks_in,2)) ! output: draw coefficient from primary free water content [/day] 
  real(dp)             :: dt                               ! time step of the simulation [hr]

  dt=24.0_dp ! unit: hr
  associate(g1=>gammaPar(ixPar%zpk1gamma1))
  where ( ks_in/=dmiss .and. myu_in/=dmiss .and. h_in/= dmiss ) 
    zpk=1.0_dp-exp(-1.0_dp*g1**2.0_dp*pi*ks_in*60.0_dp*h_in*1000.0_dp*dt/myu_in)
  else where
    zpk=dmiss
  end where
  end associate
  return
end function 

! *********************************************************************
! PFREE parameter 
! *********************************************************************
function pfree( phi_in, wp_in, gammaPar)
  implicit none
  ! input
  real(dp), intent(in) :: phi_in(:,:)                        ! input: porosity [frac]
  real(dp), intent(in) :: wp_in(:,:)                         ! input: wilting point [frac]
  real(dp), intent(in) :: gammaPar(:)                        ! input: gamma parameter array 
  ! local 
  real(dp)             :: pfree(size(wp_in,1),size(wp_in,2)) ! output: tension water maximum [mm]

  associate(g1=>gammaPar(ixPar%pfree1gamma1))
  where ( phi_in/=dmiss .and. wp_in/=dmiss ) 
    pfree=(wp_in/phi_in)**g1
  else where 
    pfree=dmiss
  end where
  end associate
  return
end function

! *********************************************************************
! ZPERC parameter 
! *********************************************************************
function zperc( twm_in, fsm_in, zsk_in, fpm_in, zpk_in )
  implicit none
  ! input
  real(dp), intent(in) :: twm_in(:,:)   ! input: tension water maximum [mm] 
  real(dp), intent(in) :: fsm_in(:,:)   ! input: supplemental free water maximum [mm] 
  real(dp), intent(in) :: zsk_in(:,:)   ! input: flow rate from supplemental free water maximum [day-1] 
  real(dp), intent(in) :: fpm_in(:,:)   ! input: primary free water maximum [mm] 
  real(dp), intent(in) :: zpk_in(:,:)   ! input: flow rate from primary free water maximum [day-1] 
  ! output 
  real(dp)             :: zperc(size(twm_in,1),size(twm_in,2))    ! output: ratio of max and min percolation rates [day-1]
  ! local 
                       
  where ( twm_in/=dmiss .and. fsm_in/=dmiss ) 
    zperc=(twm_in+fsm_in*(1.0_dp-zsk_in)+fsm_in*(1.0_dp-zpk_in))/ &
          (fsm_in*zsk_in+fpm_in*zpk_in)
  else where
    zperc=dmiss
  end where
  return
end function

! *********************************************************************
! REXP parameter 
! *********************************************************************
function rexp( wp_in, gammaPar)
  implicit none
  ! input
  real(dp), intent(in) :: wp_in(:,:)                          ! input: wilting point [frac]
  real(dp), intent(in) :: gammaPar(:)                         ! input: gamma parameter array 
  ! local 
  real(dp)             :: rexp(size(wp_in,1),size(wp_in,2)) ! output: tension water maximum [mm]
  
  associate(g1=>gammaPar(ixPar%rexp1gamma1))
  where ( wp_in/=dmiss ) 
    rexp=sqrt(wp_in/(g1-0.001))
  else where
    rexp=dmiss
  end where
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

  ! opt 1: Cosby et al. WRR 1984
  ! opt 2: campbell & shiozawa 1994 
  associate(g1=>gammaPar(ixPar%ks1gamma1), &
            g2=>gammaPar(ixPar%ks1gamma2), &
            g3=>gammaPar(ixPar%ks1gamma3), &
            g4=>gammaPar(ixPar%ks2gamma1), &
            g5=>gammaPar(ixPar%ks2gamma2), &
            g6=>gammaPar(ixPar%ks2gamma3) )
  message="ks/"
  select case(opt)
    case(1); 
      where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
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
function bd( bd_in, gammaPar )
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: bd_in(:,:)    ! input: bd from dataset [kg/m3]
  real(dp), intent(in)  :: gammaPar(:)   ! input: gamma parameter array 
  ! output 
  ! local 
  real(dp)              :: bd(size(bd_in,1),size(bd_in,2))
  real(dp),parameter    :: bd_min=805.0_dp
  real(dp),parameter    :: bd_max=1880.0_dp
  real(dp),allocatable  :: bdslope(:,:)
  real(dp),allocatable  :: bd_temp(:,:)
  integer(i4b)          :: n1                  ! number of 1st dimension 
  integer(i4b)          :: n2                  ! number of 2nd dimension 

  n1=size(bd_in,1)
  n2=size(bd_in,2)
  allocate(bdslope(n1,n2))
  allocate(bd_temp(n1,n2))
  bdslope=0.0_dp
  bd_temp=0.0_dp
  associate(g1=>gammaPar(ixPar%bd1gamma1))
  where ( bd_in /= dmiss ) 
    bd_temp = g1*bd_in
    bdslope=(bd_temp-bd_min)/(bd_max-bd_min)
    where ( bdslope > 1.0_dp) bdslope=1.0_dp
    where ( bdslope < 0.0_dp) bdslope=0.0_dp
    bd = bdslope*(bd_max-bd_min)+bd_min
  else where
    bd = dmiss 
  end where
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

  associate(g1=>gammaPar(ixPar%phi1gamma1), &
            g2=>gammaPar(ixPar%phi1gamma2), &
            g3=>gammaPar(ixPar%phi1gamma3), &
            g4=>gammaPar(ixPar%phi2gamma1), &
            g5=>gammaPar(ixPar%phi2gamma2), &
            g6=>gammaPar(ixPar%phi2gamma3), &
            g7=>gammaPar(ixPar%phi2gamma4), &
            g8=>gammaPar(ixPar%phi2gamma5), &
            g9=>gammaPar(ixPar%phi2gamma6))
    select case(opt)
      case(1);  ! Cosby
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
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
  associate(g1=>gammaPar(ixPar%fc1gamma1))
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
  associate(g1=>gammaPar(ixPar%wp1gamma1))
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
  associate(g1=>gammaPar(ixPar%b1gamma1), &
            g2=>gammaPar(ixPar%b1gamma2), &
            g3=>gammaPar(ixPar%b1gamma2))
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
  associate(g1=>gammaPar(ixPar%psis1gamma1), &
            g2=>gammaPar(ixPar%psis1gamma2), &
            g3=>gammaPar(ixPar%psis1gamma2))
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
  associate(g1=>gammaPar(ixPar%myu1gamma1), &
            g2=>gammaPar(ixPar%myu1gamma2))
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

end module soiltf 
