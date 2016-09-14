module soiltf

! Compute model soil parameter transfer function 

USE nrtype                                        ! variable types, etc.
USE data_type                                     ! Including custum data structure definition
USE multconst                                     ! Including common constant (physical constant, other e.g., missingVal, etc.)
USE def_soiltfParm                                ! define/read transfer function parameters
USE var_lookup,only:ixSoilParVic,nSoilParVic      ! index of Model soil parameter variables and number of variables
USE var_lookup,only:ixSoilParSumma,nSoilParSumma  ! index of Model soil parameter variables and number of variables
USE var_lookup,only:ixVarTopo,nVarTopo            ! index of soil polygon variables and number of variables
USE var_lookup,only:ixPar,nPrpSoil                ! index of soil properties and number of properties

implicit none

!Following accessible outside this module
public::comp_soil_model_param
!anything else
private

contains

! ********************************************************************************************
! Public subroutine: Execute computation of model-dependant soil parameter transfer function  
! *********************************************************************************************
subroutine comp_soil_model_param(ParSxySz,   &  ! in/output: soil parameter values for all soil polygons that are input
                                 sdata,      &
                                 nSLyr,      &  ! input: number of soil layers
                                 nSPoly,     &  ! input: number of soil polygons
                                 ierr, message ) ! output: error handling  

  implicit none

  ! input
  type(namevar), intent(in)           :: sdata(:)        ! storage of soil data strucuture
  integer(i4b),intent(in)             :: nSLyr           ! number of soil layer
  integer(i4b),intent(in)             :: nSPoly          ! number of soil polygons
  ! Output
  type(dat_d2d),intent(inout)         :: ParSxySz(:)     ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  integer(i4b),intent(out)            :: ierr            ! error code
  character(len=strLen),intent(out)   :: message         ! error message for current routine
  ! Local 
  integer(i4b)                        :: iSLyr           ! Loop index of soil layer
  integer(i4b)                        :: iPoly           ! Loop index of polygon 
  integer(i4b)                        :: iparm           ! Loop index of model parameters (e.g., VIC)
  real(dp),pointer                    :: xPar(:)         ! pointer to soil parameters in model space 
  logical(lgt),allocatable            :: checkDone(:)    ! used to check if the VIC parameter is processed

  ! Initialize error control
  ierr=0; message="comp_soil_model_param/"

  allocate(checkDone(nSoilParModel),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'error allocating for checkDone'; return; endif
  checkDone(:) = .false.
  sclass  => sdata(ixVarSoilData%soilclass)%ivar2 
  sand    => sdata(ixVarSoilData%sand)%dvar2 
  silt    => sdata(ixVarSoilData%silt)%dvar2 
  clay    => sdata(ixVarSoilData%clay)%dvar2 
  bd      => sdata(ixVarSoilData%bd)%dvar2 
  elestd  => sdata(ixVarSoilData%ele_std)%dvar1 
  elemean => sdata(ixVarSoilData%ele_mean)%dvar1 
  slpmean => sdata(ixVarSoilData%slp_mean)%dvar1 

  do iparm = 1,nSoilParModel
    do iSLyr = 1,nSlyr
     xPar => ParSxySz(iparm)%dat(iSLyr,:) 
     select case(iparm) !!!  neeed to change this
       case(ixPar%ks);
         checkDone(iparm)=.true.
         call xPar = ks( sand, clay, ptfopt(1), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixParl%bd);
         checkDone(iparm)=.true.
         call xPar  = bd( bd, ptfopt(2), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixParl%phi);
         if(.not.checkDone(ixPar%bd))then; ierr=20; message=trim(message)//'need to process bd before phi'; return; endif
         checkDone(iparm)=.true.
         call  xPar = phi( sand, clay, bd, ptfopt(3), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixParl%b);
         checkDone(iparm)=.true.
         call  xPar = ret_curve( sand, clay, ptfopt(6), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixPar%psis);
         checkDone(iparm)=.true.
         call xPar = psis( sand, silt, ptfopt(7), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixPar%fc);
         if(.not.checkDone(ixPar%psis))then; ierr=20; message=trim(message)//'need to process psis before fc'; return; endif
         if(.not.checkDone(ixPar%phi))then;  ierr=20; message=trim(message)//'need to process phi before fc'; return; endif
         if(.not.checkDone(ixPar%b))then;    ierr=20; message=trim(message)//'need to process b before fc'; return; endif
         checkDone(iparm)=.true.
         call  xPar = fc(sand, sprp(ixPar%phi)%varData, sprp(ixPar%psis)%varData, sprp(ixPar%b)%varData,  ptfopt(4), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixPar%wp);
         if(.not.checkDone(ixPar%psis))then; ierr=20; message=trim(message)//'need to process psis before wp'; return; endif
         if(.not.checkDone(ixPar%phi))then;  ierr=20; message=trim(message)//'need to process phi before wp'; return; endif
         if(.not.checkDone(ixPar%b))then;    ierr=20; message=trim(message)//'need to process b before wp'; return; endif
         checkDone(iparm)=.true.
         call  xPar = wp( sprp(ixPar%phi)%varData, sprp(ixPar%psis)%varData, sprp(ixPar%b)%varData, ptfopt(5), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixPar%myu);
         checkDone(iparm)=.true.
         call xpar= myu( sprp(ixPar%phi)%varData, sprp(ixPar%fc)%varData, ptfopt(8), ierr, message)
         if(ierr/=0)then; return;endif
       case(ixPar%binfilt);
         checkDone(iparm)=.true.
         call xpar=infilt( topoSxy(ixvarTopo%ele_std)%varData )
       case(ixPar%D1);
         if(.not.checkDone(ixPar%ks))then; ierr=20; message=trim(message)//'need to process "ksat" before "D1"';return; endif
         checkDone(iparm)=.true. 
         call xPar=D1( topoSxy(ixVarTopo%slp_mean)%varData,            &
                       ParSxySz(ixPar%ks)%dat(iSLyr,:),       &
                       sprpSxySz(ixPar%porosity)%varData(iSLyr,:), & 
                       sprpSxySz(ixPar%h)%varData(iSLyr,:) )
       case(ixPar%Ds);
         if(.not.checkDone(ixPar%D1))then;    ierr=20; message=trim(message)//'need to process "D1" before "Ds"'; return; endif
         if(.not.checkDone(ixPar%D3))then;    ierr=20; message=trim(message)//'need to process "D3" before "Ds"'; return; endif
         if(.not.checkDone(ixPar%Dsmax))then; ierr=20; message=trim(message)//'need to process "Dsmax" before "Ds"'; return; endif
         checkDone(iparm)=.true. 
         call xPar=Ds( ParSxySz(ixPar%D1)%dat(iSLyr,:),       &
                       ParSxySz(ixPar%D3)%dat(iSLyr,:),       & 
                       ParSxySz(ixPar%Dsmax)%dat(iSLyr,:) )
       case(ixPar%D4);
         checkDone(iparm)=.true.
         call xPar=D4(xPar,(/(a_D4_vic,iPoly=1,nSpoly)/))
       case(ixPar%c);
         if(.not.checkDone(ixPar%D4))then; ierr=20; message=trim(message)//'need to process "D4" before "c"'; return; endif
         checkDone(iparm)=.true.
         call xPar=cexpt(xPar,ParSxySz(ixPar%D4)%dat(iSLyr,:))
       case(ixPar%SD);
         checkDone(iparm)=.true.
         call xPar=soilDensity(xPar,sprpSxySz(ixPar%soil_density)%varData(iSLyr,:), a_soilDensity_vic, dmiss)
       case(ixPar%expt);
         checkDone(iparm)=.true.
         call xPar=expt(xPar,sprpSxySz(ixPar%slope_ret_curve)%varData(iSLyr,:), a_expt_vic, b_expt_vic, dmiss)
       case(ixPar%ks);
         checkDone(iparm)=.true.
         call xPar=ks( sprpSxySz(ixPar%ks)%varData(iSLyr,:) )
       case(ixPar%init_moist);
         checkDone(iparm)=.true.
         call xPar=initMoist(xPar,sprpSxySz(ixPar%porosity)%varData(iSLyr,:), &
                                      sprpSxySz(ixPar%h)%varData(iSLyr,:),        &
                                      a_initmoist_vic, dmiss)
       case(ixPar%D2);
         if(.not.checkDone(ixPar%ks))then; ierr=20; message=trim(message)//'need to process "ksat" before "D2"';return;endif
         if(.not.checkDone(ixPar%D4))then;   ierr=20; message=trim(message)//'need to process "D4" before "D2"'; return;endif
         checkDone(iparm)=.true.
         call xPar=D2( topoSxy(ixVarTopo%slp_mean)%varData,         &
                       ParSxySz(ixPar%ks)%dat(iSLyr,:), &
                       ParSxySz(ixPar%D4)%dat(iSLyr,:) )
                              
       case(ixPar%Dsmax);
         if(.not.checkDone(ixPar%D1))then; ierr=20; message=trim(message)//'need to process "D1" before "Dsmax"'; return;endif
         if(.not.checkDone(ixPar%D2))then; ierr=20; message=trim(message)//'need to process "D2" before "Dsmax"'; return;endif
         if(.not.checkDone(ixPar%D3))then; ierr=20; message=trim(message)//'need to process "D3" before "Dsmax"'; return;endif
         if(.not.checkDone(ixPar%c))then;  ierr=20; message=trim(message)//'need to process "c" before "Dsmax"'; return;endif
         checkDone(iparm)=.true. 
         call xPar=Dsmax( ParSxySz(ixPar%D1)%dat(iSLyr,:),         &
                                   ParSxySz(ixPar%D2)%dat(iSLyr,:),         & 
                                   ParSxySz(ixPar%D3)%dat(iSLyr,:),         & 
                                   ParSxySz(ixPar%c)%dat(iSLyr,:),          & 
                                   sprpSxySz(ixPar%porosity)%varData(iSLyr,:), &
                                   sprpSxySz(ixPar%h)%varData(iSLyr,:) )
                                   
       case(ixPar%bbl);
         if(.not.checkDone(ixPar%expt))then; ierr=20; message=trim(message)//'need to process "expt" before "bubble"'; return;endif
         checkDone(iparm)=.true.
         call xPar=bubble( ParSxySz(ixPar%expt)%dat(iSLyr,:))
       case(ixPar%BD);
         checkDone(iparm)=.true.
         call xPar=bulkDensity( sprpSxySz(ixPar%BD)%varData(iSLyr,:))
       case(ixPar%WcrFrac);
         checkDone(iparm)=.true.
         call xPar=WcrFrac( sprpSxySz(ixPar%field_capacity)%varData(iSLyr,:),&
                                    sprpSxySz(ixPar%porosity)%varData(iSLyr,:),      &
                                    a_WcrFrac_vic,dmiss)
       case(ixPar%WpwpFrac);
         checkDone(iparm)=.true.
         call xPar=WpwpFrac( sprpSxySz(ixPar%wilting_point)%varData(iSLyr,:),&
                                     sprpSxySz(ixPar%porosity)%varData(iSLyr,:),     &
                                     a_WpwpFrac_vic,dmiss)
       case(ixPar%D3);
         checkDone(iparm)=.true.
         call xPar=D3( xPar,sprpSxySz(ixPar%field_capacity)%varData(iSLyr,:), &
                               sprpSxySz(ixPar%h)%varData(iSLyr,:),              &
                               a_D3_vic, dmiss)
       case(ixPar%Ws);
         if(.not.checkDone(ixSoilParVic%D3))then; ierr=20; message=trim(message)//'need to process "D3" before "Dsmax"'; return;endif
         checkDone(iparm)=.true. 
         call xPar=Ws( ParSxySz(ixSoilParVic%D3)%dat(iSLyr,:),         &
                                sprpSxySz(ixPar%porosity)%varData(iSLyr,:), &
                                sprpSxySz(ixPar%h)%varData(iSLyr,:),        &
                                dmiss)
     end select ! end of parameter case
   end do ! end of layer loop
  end do ! end of parameter loop
  deallocate(checkDone,stat=ierr); 
  if(ierr/=0)then; message=trim(message)//'problem deallocating space for checkDone(nSoilParVic)';return;endif

end subroutine comp_soil_model_param

! *********************************************************************
! VIC model transfer functions
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
  real(dp), intent(in)  :: a_infilt
  real(dp), intent(in)  :: b_infilt
  real(dp), intent(in)  :: elestd_in(:)
  ! output 
  ! local 
  real(dp)              :: infilt(:) 
  real(dp),parameter    :: infilt_min=0.03_dp
  real(dp),parameter    :: infilt_max=0.50_dp

  where ( elestd_in /= dmiss ) 
    infilt = (log(elestd_in+verySmall)-a_infilt)/(log(elestd_in+verySmall)+b_infilt*10) !transfer function 
  else where
    infilt = dmiss 
  end where
 ! cap value with upper and lower bounds 
  where ( infilt > infilt_max ) infilt=infilt_max 
  where ( infilt > 0._dp .and. infilt < infilt_min ) infilt=infilt_min 

end function infilt 

! *********************************************************************
! residual_moist parameter 
! *********************************************************************
function residMoist(resid_moist)
 ! Define variables
 implicit none
 ! input
 ! output 
 real(dp), intent(out) :: resid_moist(:)
 ! local 
 
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
  allocate(S(nSpoly),stat=err); if(err/=0) call handle_err(err,'problem allocating space for S')
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
function Dsmax(Dsmax,        & ! output: Dsmax parameter for Arno baseflow [mm day^-1] 
                          D1,           & ! input:  Nijssen baseflow D1 parameter [day^-1]
                          D2,           & ! input:  Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
                          D3,           & ! input:  Nijssen baseflow D3 parameter [mm]
                          c,            & ! input:  c parameter [mm]
                          phi_in,       & ! input:  porosity [cm^3/cm^-3]
                          h_in,         & ! input:  Soil layer thickness [m]
                          dmiss)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: D1(:)        ! Nijssen baseflow D1 parameter [day^-1]
  real(dp), intent(in)  :: D2(:)        ! Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
  real(dp), intent(in)  :: D3(:)        ! Nijssen baseflow D3 parameter [mm]
  real(dp), intent(in)  :: c(:)         ! c parameter [mm]
  real(dp), intent(in)  :: phi_in(:)    ! porosity [cm^3/cm^-3]
  real(dp), intent(in)  :: h_in(:)      ! Soil layer thickness [m]
  real(dp), intent(in)  :: dmiss 
  ! output 
  real(dp), intent(out) :: Dsmax(:)     ! Dsmax parameter for Arno baseflow [mm day-1]
  ! local 
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
function D3(D3,           & ! output: D3 parameter for Nijssen baseflow [mm]
                       fc_in,        & ! input:  Field capacity [cm^3/cm^-3]
                       h_in,         & ! input:  Soil layer thickness [m]
                       a_D3,         & ! input:  TF parameter
                       dmiss)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: fc_in(:)
  real(dp), intent(in)  :: h_in(:)
  real(dp), intent(in)  :: a_D3 
  real(dp), intent(in)  :: dmiss 
  ! output 
  real(dp), intent(out) :: D3(:)
  ! local 
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
function Ws(Ws,           & ! output: Ws parameter for Arno baseflow [frac]
                       D3,           & ! input:  D3 parameter [mm]
                       phi_in,     & ! input:  porosity [cm^3/cm^-3]
                       h_in,         & ! input:  Soil layer thickness [m]
                       dmiss)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: D3(:)
  real(dp), intent(in)  :: phi_in(:)
  real(dp), intent(in)  :: h_in(:)
  real(dp), intent(in)  :: dmiss 
  ! output 
  real(dp), intent(out) :: Ws(:)
  ! local 
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
function D4(D4, a_D4)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: a_D4(:) 
  ! output 
  real(dp), intent(out) :: D4(:)
  ! local 
  
  D4 = a_D4

end function D4 

! ***********
!  c parameter
! *********************************************************************
function cexpt(c,D4)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)  :: D4(:)
  ! output 
  real(dp), intent(out) :: c(:)
  ! local 
  
  c = D4 

end function cexpt 

! ***********
! computing expt parameter 
! *********************************************************************
function expt(expt, b_in, a_expt, b_expt, dmiss)
! Define variables
implicit none
! input
real(dp), intent(in)   :: b_in(:)
real(dp), intent(in)   :: a_expt 
real(dp), intent(in)   :: b_expt 
real(dp), intent(in)   :: dmiss 
! output 
real(dp), intent(out)  :: expt(:)
! local 

where ( b_in /= dmiss ) 
  expt = a_expt+b_expt*b_in
else where
  expt = dmiss 
end where

end function expt 

! ************
! omputing Ksat parameter 
! *********************************************************************
function Ksat(Ksat, Ks_in, a_Ks, dmiss)
! Define variables
implicit none
! input
real(dp), intent(in)   :: Ks_in(:) 
real(dp), intent(in)   :: a_Ks
real(dp), intent(in)   :: dmiss 
! output 
real(dp), intent(out)  :: Ksat(:)
! local 

where ( Ks_in /= dmiss ) 
  ksat = a_Ks*(Ks_in*10*24)
else where
  ksat = dmiss 
end where

end function Ksat 

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
! bulk_density parameter -VIC 
! *********************************************************************
function bd( bd_in, a_bulkDensity, dmiss)
! Define variables
implicit none
! input
real(dp), intent(in)  :: bd_in(:) 
! output 
! local 
real(dp), intent(out) :: bd(:)
real(dp),parameter    :: rho_bulk_min=805.0_dp
real(dp),parameter    :: rho_bulk_max=1880.0_dp
real(dp),allocatable  :: beta(:)
real(dp),allocatable  :: bulkDensity_temp(:)
integer(i4b)          :: nSpoly                    ! number of element
integer(i4b)          :: err                       ! error code

nSpoly=size(rho_in)
allocate(beta(nSpoly),stat=err); if(err/=0) call handle_err(err,'problem allocating space for PolyID')
allocate(bulkDensity_temp(nSpoly),stat=err); if(err/=0) call handle_err(err,'problem allocating space for PolyID')

beta(:)=0.0_dp
bulkDensity_temp(:)=0.0_dp

where ( rho_in /= dmiss ) 
  bulkDensity_temp = a_bulkDensity*rho_in
  beta=(bulkDensity_temp-rho_bulk_min)/(rho_bulk_max-rho_bulk_min)
  where ( beta > 1.0_dp) beta=1.0_dp
  where ( beta < 0.0_dp) beta=0.0_dp
  bulkDensity = beta*(rho_bulk_max-rho_bulk_min)+rho_bulk_min
else where
  bulkDensity = dmiss 
end where

end function bulkDensity 

! ***********
! soil_density parameter 
! *********************************************************************
function soilDensity( srho_in)
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

! ***********
! soil porosity
! *********************************************************************
function phi(phi_in)
! Define variables
implicit none
! input
real(dp), intent(in)  :: phi_in(:)      ! Porosity
! local 
real(dp)              :: phi(:)

where ( phi_in /= dmiss ) 
  phi = a_porosity* phi_in
else where
  phi = dmiss 
end where

end function phi 

end module soiltf 
