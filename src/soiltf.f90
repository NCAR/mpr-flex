module soiltf

! Compute model soil parameter transfer function 

USE nrtype                                        ! variable types, etc.
USE data_type                                     ! Including custum data structure definition
USE public_var                                     ! Including common constant (physical constant, other e.g., missingVal, etc.)

implicit none

!Following accessible outside this module
public::comp_soil_model_param
!anything else
private

contains

! ********************************************************************************************
! Public subroutine: Execute computation of model-dependant soil parameter transfer function  
! *********************************************************************************************
subroutine comp_soil_model_param(ParSxySz,          &  ! in/output: soil parameter values for all soil polygons that are input
                                 sdata,             &  ! input: soil data
                                 gammaParMasterMeta,&  ! input: gamma parameter meta file - val of calibrating parammeter  is adjusted via calibration
                                 nSLyr,             &  ! input: number of soil layers
                                 nSPoly)               ! input: number of soil polygons

  use globalData, only:betaInGamma
  use get_ixname, only:get_ixPar
  use var_lookup, only:nPar

  implicit none

  ! in/out
  type(dat_d2d),intent(inout)         :: ParSxySz(:)            ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  ! input
  type(namevar), intent(in)           :: sdata(:)               ! storage of soil data strucuture
  type(par_meta)                      :: gammaParMasterMeta(:)
  integer(i4b),  intent(in)           :: nSLyr                  ! number of soil layer
  integer(i4b),  intent(in)           :: nSPoly          ! number of soil polygons
  ! Local 
  integer(i4b)                        :: ierr            ! error code
  character(len=strLen)               :: message         ! error message for current routine
  type(dat_d2d),intent(inout)         :: ParTemp(:)      ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  integer(i4b)                        :: iSLyr           ! Loop index of soil layer
  integer(i4b)                        :: iPoly           ! Loop index of polygon 
  integer(i4b)                        :: iparm           ! Loop index of model parameters (e.g., VIC)
  real(dp),pointer                    :: xPar(:)         ! pointer to soil parameters in model space 
  logical(lgt),allocatable            :: checkDone(:)    ! used to check if the VIC parameter is processed

  message="comp_soil_model_param/"
  associate(sclass  => sdata(ixVarSoilData%soilclass)%ivar2, &
            sand    => sdata(ixVarSoilData%sand)%dvar2, & 
            silt    => sdata(ixVarSoilData%silt)%dvar2, & 
            clay    => sdata(ixVarSoilData%clay)%dvar2, & 
            bulkd   => sdata(ixVarSoilData%bd)%dvar2, & 
            elestd  => sdata(ixVarSoilData%ele_std)%dvar1, & 
            elemean => sdata(ixVarSoilData%ele_mean)%dvar1, & 
            slpmean => sdata(ixVarSoilData%slp_mean)%dvar1)
  do iParm = 1,nPar
    do iSLyr = 1,nSlyr
    xPar => ParTemp(iparm)%dat(iSLyr,:) 
     select case(iParm)
       case(ixPar%ks);
         checkDone(iparm)=.true.
         call xPar = ks( sand, clay, ptfopt(1))
       case(ixParl%bd);
         checkDone(iparm)=.true.
         call xPar  = bd( bulkd, ptfopt(2))
       case(ixParl%phi);
         if(.not.checkDone(ixPar%bd)) stop trim(message)//'need to process bd before phi'
         checkDone(iparm)=.true.
         call xPar = phi( sand, clay, ParTemp(ixPar%bd)%dat(iSlyr,:), ptfopt(3))
       case(ixParl%b);
         checkDone(iparm)=.true.
         call xPar = ret_curve( sand, clay, ptfopt(6))
       case(ixPar%psis);
         checkDone(iparm)=.true.
         call xPar = psis( sand, silt, ptfopt(7))
       case(ixPar%fc);
         if(.not.checkDone(ixPar%psis)) stop trim(message)//'need to process psis before fc'
         if(.not.checkDone(ixPar%phi))  stop trim(message)//'need to process phi before fc'
         if(.not.checkDone(ixPar%b))    stop trim(message)//'need to process b before fc'
         checkDone(iparm)=.true.
         call  xPar = fc(sand, ParTemp(ixPar%phi)%dat, ParTemp(ixPar%psis)%dat, ParTemp(ixPar%b)%dat,  ptfopt(4))
       case(ixPar%wp);
         if(.not.checkDone(ixPar%psis)) stop trim(message)//'need to process psis before wp'
         if(.not.checkDone(ixPar%phi))  stop trim(message)//'need to process phi before wp'
         if(.not.checkDone(ixPar%b))    stop trim(message)//'need to process b before wp'
         checkDone(iparm)=.true.
         call  xPar = wp( ParTemp(ixPar%phi)%dat, ParTemp(ixPar%psis)%dat, ParTemp(ixPar%b)%dat, ptfopt(5))
       case(ixPar%myu);
         checkDone(iparm)=.true.
         call xpar= myu( ParTemp(ixPar%phi)%dat, ParTemp(ixPar%fc)%dat, ptfopt(8))
       case(ixPar%binfilt);
         checkDone(iparm)=.true.
         call xpar=infilt( elestd )
       case(ixPar%D1);
         if(.not.checkDone(ixPar%ks)) stop trim(message)//'need to process "ksat" before "D1"'
         checkDone(iparm)=.true. 
         call xPar=D1( slpmean,            &
                       ParTemp(ixPar%ks)%dat(iSLyr,:),       &
                       ParTemp(ixPar%phi)%dat(iSLyr,:), & 
                       hslyr(iSLyr,:) )
                       !ParTemp(ixPar%h)%dat(iSLyr,:) )
       case(ixPar%Ds);
         if(.not.checkDone(ixPar%D1))    stop trim(message)//'need to process "D1" before "Ds"'
         if(.not.checkDone(ixPar%D3))    stop trim(message)//'need to process "D3" before "Ds"'
         if(.not.checkDone(ixPar%Dsmax)) stop trim(message)//'need to process "Dsmax" before "Ds"'
         checkDone(iparm)=.true. 
         call xPar=Ds( ParTemp(ixPar%D1)%dat(iSLyr,:),       &
                       ParTemp(ixPar%D3)%dat(iSLyr,:),       & 
                       ParTemp(ixPar%Dsmax)%dat(iSLyr,:) )
       case(ixPar%D4);
         checkDone(iparm)=.true.
         call xPar=D4()
       case(ixPar%c);
         if(.not.checkDone(ixPar%D4)) stop trim(message)//'need to process "D4" before "c"'
         checkDone(iparm)=.true.
         call xPar=cexpt(ParTemp(ixPar%D4)%dat(iSLyr,:))
       case(ixPar%SD);
         checkDone(iparm)=.true.
         call xPar=soilDensity(ParTemp(ixPar%soil_density)%dat(iSLyr,:))
       case(ixPar%expt);
         checkDone(iparm)=.true.
         call xPar=expt(ParTemp(ixPar%slope_ret_curve)%dat(iSLyr,:))
       case(ixPar%init_moist);
         checkDone(iparm)=.true.
         call xPar=initMoist(xPar,ParTemp(ixPar%porosity)%dat(iSLyr,:), &
                                      hslyr(iSLyr,:) )
       case(ixPar%D2);
         if(.not.checkDone(ixPar%ks)) stop trim(message)//'need to process "ksat" before "D2"'
         if(.not.checkDone(ixPar%D4)) stop trim(message)//'need to process "D4" before "D2"'
         checkDone(iparm)=.true.
         call xPar=D2( topoSxy(ixVarTopo%slp_mean)%dat,  &
                       ParTemp(ixPar%ks)%dat(iSLyr,:),   &
                       ParTemp(ixPar%D4)%dat(iSLyr,:) )
                              
       case(ixPar%Dsmax);
         if(.not.checkDone(ixPar%D1)) stop trim(message)//'need to process "D1" before "Dsmax"'
         if(.not.checkDone(ixPar%D2)) stop trim(message)//'need to process "D2" before "Dsmax"'
         if(.not.checkDone(ixPar%D3)) stop trim(message)//'need to process "D3" before "Dsmax"'
         if(.not.checkDone(ixPar%c))  stop trim(message)//'need to process "c" before "Dsmax"'
         checkDone(iparm)=.true. 
         call xPar=Dsmax( ParTemp(ixPar%D1)%dat(iSLyr,:),         &
                                   ParTemp(ixPar%D2)%dat(iSLyr,:),         & 
                                   ParTemp(ixPar%D3)%dat(iSLyr,:),         & 
                                   ParTemp(ixPar%c)%dat(iSLyr,:),          & 
                                   ParTemp(ixPar%porosity)%dat(iSLyr,:), &
                                   hslyr(iSLyr,:) )
                                   
       case(ixPar%bbl);
         if(.not.checkDone(ixPar%expt)) stop trim(message)//'need to process "expt" before "bubble"'
         checkDone(iparm)=.true.
         call xPar=bubble( ParTemp(ixPar%expt)%dat(iSLyr,:))
       case(ixPar%WcrFrac);
         checkDone(iparm)=.true.
         call xPar=WcrFrac( ParTemp(ixPar%fc)%dat(iSLyr,:),&
                            ParTemp(ixPar%phi)%dat(iSLyr,:) )
       case(ixPar%WpwpFrac);
         checkDone(iparm)=.true.
         call xPar=WpwpFrac( ParTemp(ixPar%wp)%dat(iSLyr,:),&
                             ParTemp(ixPar%phi)%dat(iSLyr,:))
       case(ixPar%D3);
         checkDone(iparm)=.true.
         call xPar=D3( xPar,ParTemp(ixPar%fc)%dat(iSLyr,:), &
                        hslyr(iSLyr,:))
       case(ixPar%Ws);
         if(.not.checkDone(ixSoilParVic%D3)) stop trim(message)//'need to process "D3" before "Dsmax"'
         checkDone(iparm)=.true. 
         call xPar=Ws( ParTemp(ixSoilParVic%D3)%dat(iSLyr,:), &
                       ParTemp(ixPar%phi)%dat(iSLyr,:),       &
                       hslyr(iSLyr,:))
     end select ! end of parameter case
   end do ! end of layer loop
  end do ! end of parameter loop

  do ipar=1,nSoilParModel
    idBeta=get_ixPar(betaInGamma(iParm))
    ParSxySz(ipar)=ParTemp(idBeta) 
  enddo
  return
end subroutine comp_soil_model_param


! *********************************************************************
! Library of hydrologic model transfer functions
! *********************************************************************
! infilt parameter 
! *********************************************************************
function infilt(elestd_in)
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

  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: elestd_in(:)
  ! output 
  ! local 
  real(dp)              :: infilt(:) 
  real(dp),parameter    :: infilt_min=0.03_dp
  real(dp),parameter    :: infilt_max=0.50_dp
  
  associate(g1=>gammaParMasterMeta(ixPar%binfilt1gamma1)%val, &
            g2=>gammaParMasterMeta(ixPar%binfilt1gamma2)%val)
    where ( elestd_in /= dmiss ) 
      infilt = (log(elestd_in+verySmall)-g1)/(log(elestd_in+verySmall)+g2*10) !transfer function 
    else where
      infilt = dmiss 
    end where
     !cap value with upper and lower bounds 
    where ( infilt > infilt_max ) infilt=infilt_max 
    where ( infilt > 0._dp .and. infilt < infilt_min ) infilt=infilt_min 
  end associate
  return
end function infilt 

! *********************************************************************
! residual_moist parameter 
! *********************************************************************
function residMoist()
 ! Define variables
 implicit none
 ! input
 ! output 
 ! local 
 real(dp) :: resid_moist(:)
 
 resid_moist = 0._dp
  
end function residMoist
         
! ***********
!  Nijssen basiflow D1 parameter
! *********************************************************************
function D1(slope_in, ks_in, phi_in, h_in)
! Define variables
 implicit none

! input
 real(dp), intent(in)  :: slope_in(:)   ! mean slope [%]
 real(dp), intent(in)  :: Ks_in(:)      ! Ksat [mm/s]
 real(dp), intent(in)  :: phi_in(:)     ! porosity [-]
 real(dp), intent(in)  :: h_in(:)       ! layer thickness [m]
! output 
! local 
 real(dp),             :: D1(:)
 integer(i4b)          :: nSpoly            ! number of element
 real(dp), allocatable :: S(:)              ! length scaling term [mm]: 1, Max. soil storage etc
 real(dp), parameter   :: D1_min=0.0001_dp
 real(dp), parameter   :: D1_max=1.0_dp
 integer(i4b)          :: err                       ! error code
  
 ! local variable allocation
  nSpoly=size(D1)
  allocate(S(nSpoly),stat=err) 
  S(:)=phi_in*h_in*1000
  S(:)=1.0_dp
 ! compute parameters 
 where ( slope_in /= dmiss .and. Ks_in /= dmiss )
   D1 = S**(-1)*10**(-1*a_D1)*Ks_in*(slope_in*0.01)
 else where
   D1 = dmiss
 end where
 ! cap value with upper and lower bounds 
 where ( D1 > D1_max ) D1=D1_max 
 where ( D1 > 0._dp .and. D1 < D1_min ) D1=D1_min 

end function D1 

! ***********
!  Arno basiflow Ds parameter
! *********************************************************************
function Ds( D1, D3, Dsmax)
! Define variables
 implicit none

! input
 real(dp), intent(in)  :: D1(:)       ! nijssen baseflow D1 parameter [day-1]
 real(dp), intent(in)  :: D3(:)       ! nijssen baseflow D3 parameter [mm]
 real(dp), intent(in)  :: Dsmax(:)    ! ARNO Dsmax parameter [mm/day]
! output 
! local 
 real(dp)              :: Ds(:)
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

end function Ds 

! ***********
! Nijssen baseflow D2 parameter
! *********************************************************************
function D2(slope_in, ks_in, D4_in) 
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: slope_in(:) ! slope percent
  real(dp), intent(in)  :: Ks_in(:)    ! ksat [mm/s]
  real(dp), intent(in)  :: D4_in(:)    ! VIC D4 paramet
  ! output 
  ! local 
  real(dp),             :: D2(:)
  integer(i4b)          :: nSpoly           ! number of element
  real(dp), allocatable :: S(:)             ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp), parameter   :: D2_min=0.0001_dp
  real(dp), parameter   :: D2_max=1.0_dp
  integer(i4b)          :: err                       ! error code
  
 ! local variable allocation
  nSpoly=size(D2)
  allocate(S(nSpoly),stat=err)
  S(:)=1.0_dp
 ! compute parameters 
  where ( slope_in /= dmiss .and. Ks_in /= dmiss )
    D2 = S**(-1*D4_in)*10**(-1*a_D2)*Ks_in*(slope_in*0.01)
  else where
    D2 = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( D2 > D2_max ) D2=D2_max
  where ( D2 > 0._dp .and. D2 < D2_min ) D2=D2_min

end function D2 

! ***********
! Arno baseflow Dsmax parameter
! *********************************************************************
function Dsmax( D1,           & ! input:  Nijssen baseflow D1 parameter [day^-1]
                D2,           & ! input:  Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
                D3,           & ! input:  Nijssen baseflow D3 parameter [mm]
                c,            & ! input:  c parameter [mm]
                phi_in,       & ! input:  porosity [cm^3/cm^-3]
                h_in)          ! input:  Soil layer thickness [m]
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: D1(:)        ! Nijssen baseflow D1 parameter [day^-1]
  real(dp), intent(in)  :: D2(:)        ! Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
  real(dp), intent(in)  :: D3(:)        ! Nijssen baseflow D3 parameter [mm]
  real(dp), intent(in)  :: c(:)         ! c parameter [mm]
  real(dp), intent(in)  :: phi_in(:)    ! porosity [cm^3/cm^-3]
  real(dp), intent(in)  :: h_in(:)      ! Soil layer thickness [m]
  ! output 
  ! local 
  real(dp)              :: Dsmax(:)     ! Dsmax parameter for Arno baseflow [mm day-1]
  real(dp), parameter   :: Dsmax_min=0.1_dp
  real(dp), parameter   :: Dsmax_max=30.0_dp

  where ( phi_in /= dmiss .and. h_in /= dmiss )
    Dsmax = D2*(phi_in*h_in-D3)**c+D1*(phi_in*h_in)
  else where
    Dsmax = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( Dsmax > Dsmax_max ) Dsmax=Dsmax_max
  where ( Dsmax > 0._dp .and. Dsmax < Dsmax_min ) Dsmax=Dsmax_min

end function Dsmax 

! ***********
!  Nijssen baseflow D3 parameter
! *********************************************************************
function D3( fc_in,        & ! input:  Field capacity [cm^3/cm^-3]
             h_in)          ! input:  Soil layer thickness [m]
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: fc_in(:)
  real(dp), intent(in)  :: h_in(:)
  real(dp), intent(in)  :: a_D3 
  ! output 
  ! local 
  real(dp)              :: D3(:)
  real(dp), parameter   :: D3_min=0.0001_dp
  real(dp), parameter   :: D3_max=1000.0_dp

  where ( fc_in /= dmiss .and. h_in /= dmiss ) 
    D3 = a_D3* fc_in * (h_in*1000)
  else where
    D3 = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( D3 > D3_max ) D3=D3_max
  where ( D3 > 0._dp .and. D3 < D3_min ) D3=D3_min 

end function D3 

! ***********
!  Arno baseflow Ws parameter (conversion equation)
! *********************************************************************
function Ws( D3,         & ! input:  D3 parameter [mm]
             phi_in,     & ! input:  porosity [cm^3/cm^-3]
             h_in)          ! input:  Soil layer thickness [m]
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: D3(:)
  real(dp), intent(in)  :: phi_in(:)
  real(dp), intent(in)  :: h_in(:)
  ! output 
  ! local 
  real(dp)              :: Ws(:)
  real(dp), parameter   :: Ws_min=0.05_dp
  real(dp), parameter   :: Ws_max=1.0_dp

  where ( phi_in /= dmiss .and. h_in /= dmiss .and. D3 /= dmiss ) 
    Ws = D3 / phi_in / h_in
  else where
    Ws = dmiss
  end where
 ! cap value with upper and lower bounds 
  where ( Ws > Ws_max ) Ws=Ws_max
  where ( Ws > 0._dp .and. Ws < Ws_min ) Ws=Ws_min 

end function Ws 

! ***********
!  Nijssen baseflow D4 parameter
! *********************************************************************
function D4()
  ! Define variables
  implicit none
  ! input
  ! output 
  ! local 
  real(dp)    :: D4(:)
  
  associate(g1=>gammaParMasterMeta(ixPar%D41gamma1)%val)
  D4 = g1 
  end associate

end function D4 

! ***********
!  c parameter
! *********************************************************************
function cexpt( D4 )
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: D4(:)
  ! output 
  ! local 
  real(dp)              :: cexpt(:)
  
  cexpt = D4 

end function cexpt 

! ***********
! computing expt parameter 
! *********************************************************************
function expt( b_in )
! Define variables
implicit none
! input
real(dp), intent(in)   :: b_in(:)
! output 
! local 
real(dp)               :: expt(:)

where ( b_in /= dmiss ) 
  expt = a_expt+b_expt*b_in
else where
  expt = dmiss 
end where

end function expt 

! ************
! computing init_moist parameter  
! *********************************************************************
function initMoist( phi_in, h_in)

! Define variables
implicit none
! input
real(dp), intent(in)  :: phi_in(:)       ! porosity [-]
real(dp), intent(in)  :: h_in(:)         ! thickness [m]
! output 
! local  
real(dp)              :: initMoist(:)

where ( phi_in /= dmiss ) 
  initMoist = a_initmoist*phi_in*(h_in*1000.0_dp)
else where
  initMoist = dmiss 
end where

end function initMoist 

! ***********
! bubble parameter 
! *********************************************************************
function bubble( expt_in )
! Requre expt computation first
! Define variables
implicit none
! input
real(dp), intent(in)  :: expt_in(:) 
! output 
! local 
real(dp)              :: bubble(:)

where ( expt_VIC /= dmiss ) 
  bubble = a_bubble*expt_VIC+b_bubble
else where
  bubble = dmiss 
end where

end function bubble 

! ***********
! soil_density parameter 
! *********************************************************************
function soilDensity( srho_in )
! Define variables
implicit none
! input
real(dp), intent(in)  :: srho_in(:)
! output 
! local 
real(dp), intent(out) :: soilDensity(:)

where ( srho_in /= dmiss ) 
  soilDensity = a_soilDensity*srho_in
else where
  soilDensity = dmiss 
end where

end function soilDensity 

! ***********
! WcrFrac parameter  
! *********************************************************************
function WcrFrac(fc_in, phi_in)
! Define variables
implicit none
! input
real(dp), intent(in)  :: fc_in(:)
real(dp), intent(in)  :: phi_in(:)
! output 
! local 
real(dp)              :: WcrFrac(:) 

where ( fc_in /= dmiss .and. phi_in /= dmiss ) 
  WcrFrac = a_WcrFrac*fc_in/phi_in
else where
  wcrFrac = dmiss
end where

end function WcrFrac 

! ************
! WpwpFrac parameter  
! *********************************************************************
function WpwpFrac( wp_in, phi_in, a_WpwpFrac)
! Define variables
implicit none
! input
real(dp), intent(in)  :: wp_in(:)   ! Wilting point
real(dp), intent(in)  :: phi_in(:)   ! Porosity
real(dp), intent(in)  :: a_WpwpFrac 
! local 
real(dp), intent(out) :: WpwpFrac(:) 

where ( wp_in /= dmiss .and. phi_in /= dmiss ) 
  WpwpFrac = a_WpwpFrac*wp_in/phi_in
else where
  wpwpFrac = dmiss
end where

end function WpwpFrac 

! *********************************************************************
! pedo-transfer function for saturated hydraulic conductivity (ks)
! *********************************************************************
function ks( sand_in, clay_in, opt)
  implicit none
  ! input
  real(dp), intent(in)       :: sand_in(:,:) ! input: sand [percent] 
  real(dp), intent(in)       :: clay_in(:,:) ! input: clay [percent] 
  integer(i2b)               :: opt          ! input: option for transfer function form
  ! local 
  character(*)               :: message      ! error message
  integer(i4b)               :: ierr         ! error code
  real(dp)                   :: ks(:,:)      ! ks 

  ! opt 1: Cosby et al. WRR 1984
  ! opt 2: campbell & shiozawa 1994 
  message="ks/"
  select case(opt)
    case(1); 
      where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
        ks = gamma1ks + gamma2ks*sand_in + gamma3ks*clay_in
        ks = 10**ks*2.54   ! 2.54 cm/inch. Cosby give Ksat in inch/hr 
      else where
        ks = dmiss 
      end where
    case(2); 
      where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
        ks = gamma1ks*exp(gamma2ks*sand_in+gamma3ks*clay_in)
      else where
        ks = dmiss 
      end where
    case default; stop trim(message)//'opt not recognized'
  end select
end function ks

! *********************************************************************
! pedo-transfer function for bulk density 
! *********************************************************************
function bd( bd_in )
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: bd_in(:) 
  ! output 
  ! local 
  real(dp), intent(out) :: bd(:)
  real(dp),parameter    :: bd_min=805.0_dp
  real(dp),parameter    :: bd_max=1880.0_dp
  real(dp),allocatable  :: bdslope(:)
  real(dp),allocatable  :: bd_temp(:)
  integer(i4b)          :: nSpoly                    ! number of element
  integer(i4b)          :: err                       ! error code

  nSpoly=size(rho_in)
  allocate(bdslope(nSpoly),stat=err)
  allocate(bd_temp(nSpoly),stat=err)
  bdslope(:)=0.0_dp
  bd_temp(:)=0.0_dp
  associate(g1=>gammaParMasterMeta(ixPar%bd1gamma1)%val)
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
end function bd

! *********************************************************************
! pedo-transfer function for porosity 
! *********************************************************************
function phi(sand_in, clay_in, db_in, opt)      
  implicit none
  ! input
  real(dp), intent(in)       :: sand_in(:,:)   ! input: sand  
  real(dp), intent(in)       :: clay_in(:,:)   ! input: clay 
  real(dp), intent(in)       :: db_in(:,:)     ! input: bulk density 
  integer(i2b)               :: opt            ! option for transfer function form
  ! local 
  character(*)               :: message        ! error message
  integer(i4b)               :: ierr           ! error code
  real(dp)                   :: phi(:,:)  ! estimated porosity [fraction]

  ! opt 1: Cosby et al. WRR 1984
  ! opt 2: Zacharias & Wessolek 2007
  message="comp_porosity/"

  associate(g1=>gammaParMasterMeta(ixPar%phi1gamma1)%val, &
            g2=>gammaParMasterMeta(ixPar%phi1gamma2)%val, &
            g3=>gammaParMasterMeta(ixPar%phi1gamma3)%val, &
            g4=>gammaParMasterMeta(ixPar%phi2gamma1)%val, &
            g5=>gammaParMasterMeta(ixPar%phi2gamma2)%val, &
            g6=>gammaParMasterMeta(ixPar%phi2gamma3)%val, &
            g7=>gammaParMasterMeta(ixPar%phi2gamma4)%val, &
            g8=>gammaParMasterMeta(ixPar%phi2gamma5)%val, &
            g9=>gammaParMasterMeta(ixPar%phi2gamma6)%val)
    select case(opt)
      case(1);  ! Cosby
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          phi = g1+g2*sand_in+g3*clay_in
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
      case default; stop trim(message)//'opt not recognized'
    end select
  end associate
  end function phi 

! *********************************************************************
! pedo-transfer function for field capacity 
! *********************************************************************
  function fc(sand_in, phi_in, psis_in, b_in, opt)
    implicit none
    ! input
    real(dp), intent(in)       :: sand_in(:,:)   ! input: sand  
    real(dp), intent(in)       :: phi_in(:,:)    ! input: porosity [fraction]  
    real(dp), intent(in)       :: psis_in(:,:)   ! input: saturation matric potential [ 
    real(dp), intent(in)       :: b_in(:,:)      ! input: slope of cambell retention curve 
    integer(i2b)               :: opt            ! id for transfer function form
    ! output 
    ! inout
    character(*)               :: message        ! error message
    ! local 
    real(dp)                   :: fc(:,:)        ! estimated field capacity [fraction]
    real(dp),allocatable       :: psi_fc(:,:)    ! matric potential at field capacity [kPa]  
    integer(i4b)               :: nSpoly         ! number of soil polygon 
    integer(i4b)               :: nSlyr          ! number of soil layer 

    ! opt 1: Campbell 1974 
    message="fc/"

    nSpoly=size(porosity,1)
    nSlyr=size(porosity,2)
    allocate(psi_fc(nSpoly,nSlyr),stat=ierr)
    psi_fc(:,:)=-20
    where (sand_in > 69) psi_fc=-10

    select case(opt)
      case(1);  !campbell
        where ( porosity /= dmiss .and. sand_in /= dmiss .and. b /= dmiss .and. psi_sat /= dmiss ) 
          fc = gamma1fc*porosity*(psi_fc/psi_sat)**(-1/b)
        else where
          fc = dmiss 
        end where
      case default; stop trim(message)//'opt not recognized' 
    end select
  end function fc

! *********************************************************************
! pedo-transfer function for wilting point 
! *********************************************************************
  function wp( phi_in, psis_in, b_in, opt) 
    implicit none
    ! input
    real(dp), intent(in)       :: phi_in(:,:)    ! input: porosity [fraction]  
    real(dp), intent(in)       :: psis_in(:,:)   ! input: saturation matric potential [kPa]  
    real(dp), intent(in)       :: b_in(:,:)      ! input: slope of cambell retention curve [-]
    integer(i2b)               :: opt            ! input: option for transfer function form
    ! local 
    integer(i4b)               :: ierr           ! error code
    character(*)               :: message        ! error message
    real(dp)                   :: wp(:,:)        ! estimated field capacity [frac]
    real(dp),allocatable       :: psi_wp(:,:)    ! matric potential at wilting point [kPa]  
    integer(i4b)               :: nSpoly         ! number of soil polygon 
    integer(i4b)               :: nSlyr          ! number of soil layer 

    ! opt 1: Campbell 1974 
    message="wp/"

    nSpoly=size(porosity,1)
    nSlyr=size(porosity,2)
    allocate(psi_wp(nSpoly,nSlyr),stat=ierr)
    psi_wp(:,:)=-1500
    select case(opt)
      case(1);  !Cosby et al. 
        where ( porosity /= dmiss .and. b /= dmiss .and. psi_sat /= dmiss ) 
          wp = gamma1wp*porosity*(psi_wp/psi_sat)**(-1/b)
        else where
          wp = dmiss 
        end where
      case default; stop trim(message)//'opt not recognized'
  end function wp

! *********************************************************************
! pedo-transfer function for b (slope of retention curve in log space)
! *********************************************************************
function ret_curve(sand_in, clay_in, opt)
  implicit none
  ! input
  real(dp), intent(in)       :: sand_in(:,:)  ! input: sand  
  real(dp), intent(in)       :: clay_in(:,:)  ! input: clay 
  integer(i2b)               :: opt           ! input: option for transfer function form
  ! local 
  character(*)               :: message        ! error message
  real(dp)                   :: ret_curv(:,:)        ! computed [-] 

  ! opt 1: Cosby et al. WRR 1984
  message="ret_curv/"
  associate(g1=>gammaParMasterMeta(ixPar%b1gamma1)%val, &
            g2=>gammaParMasterMeta(ixPar%b1gamma2)%val, &
            g3=>gammaParMasterMeta(ixPar%b1gamma2)%val)
    select case(opt)
      case(1); 
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          ret_curv = g1+g2*sand_in+g3*clay_in
        else where
          ret_curv = dmiss 
        end where
      case default; stop trim(message)//'opt not recognized'
    end select
  end associate
end function ret_curve

! *********************************************************************
! pedo-transfer function for saturation matric potential 
! *********************************************************************
function psis( sand_in, silt_in, opt)
  implicit none
  ! input
  real(dp), intent(in)       :: sand_in(:,:)   ! input: sand percentage [percentage]
  real(dp), intent(in)       :: silt_in(:,:)   ! input: silt percentage [percentage]
  integer(i2b), intent(in)   :: opt            ! input: option for transfer function form
  ! local 
  character(*)               :: message        ! error message
  real(dp)                   :: psis(:,:)      ! output: saturation matric potential [kPa]  

  ! opt 1: Cosby et al. WRR 1984
  message="psis/"
  associate(g1=>gammaParMasterMeta(ixPar%psis1gamma1)%val, &
            g2=>gammaParMasterMeta(ixPar%psis1gamma2)%val, &
            g3=>gammaParMasterMeta(ixPar%psis1gamma2)%val)
    select case(opt)
      case(1);  !Cosby et al. 
        where ( sand_in /= dmiss .and. silt_in /= dmiss ) 
          psis = g1 + g2*sand_in + g3*silt_in
          psis = -1.0_dp*10**psi_sat*0.0980665_dp        ! 0.0980665 kPa/cm-H2O. Cosby give psi_sat in cm of water (cm-H2O)
        else where
          psis = dmiss 
        end where
      case default; stop message=trim(message)//'opt not recognized'
    end select
  end associate
end function psis

! *********************************************************************
! pedo-transfer function for specific yield  
! *********************************************************************
function myu(phi_in, fc_in, opt)
  implicit none
  ! input
  real(dp),    intent(in)    :: phi_in(:,:)     ! input: porosity [frac]
  real(dp),    intent(in)    :: fc_in(:,:)      ! input: field capacity [frac]
  integer(i2b),intent(in)    :: opt             ! input: option for transfer function form
  ! local 
  character(*)               :: message        ! error message
  real(dp)                   :: myu(:,:)        ! output: specific yield [-]  

  ! opt 1: Koren et al. 2003
  message="myu/"
  associate(g1=>gammaParMasterMeta(ixPar%myu1gamma1)%val, &
            g2=>gammaParMasterMeta(ixPar%myu1gamma2)%val)
  select case(opt)
    case(1);  ! koren
      where ( phi_in /= dmiss .and. fc_in /= dmiss ) 
        myu = g1*(phi_in-fc_in)**g2 
      else where
        myu = dmiss 
      end where
    case default; stop message=trim(message)//'opt not recognized'
  end select
  end associate
end function myu

end module soiltf 
