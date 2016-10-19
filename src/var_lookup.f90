MODULE var_lookup

 ! Define index arrays for named variables
 ! list of all the parameters including gamma and beta parameters

 USE nrtype
 USE public_var

 implicit none

 private

! ***********************************************************************************************************
! 1.Define indices for gamma (global) parameters
! ***********************************************************************************************************
 type, public  ::  iLook_Par
   !Gamma parameter
   integer(i4b)     :: ks1gamma1       = imiss  ! 
   integer(i4b)     :: ks1gamma2       = imiss  ! 
   integer(i4b)     :: ks1gamma3       = imiss  ! 
   integer(i4b)     :: ks2gamma1       = imiss  ! 
   integer(i4b)     :: ks2gamma2       = imiss  ! 
   integer(i4b)     :: ks2gamma3       = imiss  ! 
   integer(i4b)     :: phi1gamma1      = imiss  ! 
   integer(i4b)     :: phi1gamma2      = imiss  ! 
   integer(i4b)     :: phi1gamma3      = imiss  ! 
   integer(i4b)     :: phi2gamma1      = imiss  !
   integer(i4b)     :: phi2gamma2      = imiss  ! 
   integer(i4b)     :: phi2gamma3      = imiss  ! 
   integer(i4b)     :: phi2gamma4      = imiss  ! 
   integer(i4b)     :: phi2gamma5      = imiss  ! 
   integer(i4b)     :: phi2gamma6      = imiss  ! 
   integer(i4b)     :: fc1gamma1       = imiss  ! 
   integer(i4b)     :: wp1gamma1       = imiss  ! 
   integer(i4b)     :: b1gamma1        = imiss  ! 
   integer(i4b)     :: b1gamma2        = imiss  ! 
   integer(i4b)     :: b1gamma3        = imiss  ! 
   integer(i4b)     :: psis1gamma1     = imiss  ! 
   integer(i4b)     :: psis1gamma2     = imiss  ! 
   integer(i4b)     :: psis1gamma3     = imiss  ! 
   integer(i4b)     :: myu1gamma1      = imiss  ! 
   integer(i4b)     :: myu1gamma2      = imiss  ! 
   integer(i4b)     :: z1gamma1        = imiss  ! total depth mulitplier 
   integer(i4b)     :: h1gamma1        = imiss  ! fraction of top layer to total depth  
   integer(i4b)     :: h1gamma2        = imiss  ! fraction of 2nd layer to total depth
   integer(i4b)     :: binfilt1gamma1  = imiss  ! 
   integer(i4b)     :: binfilt1gamma2  = imiss  ! 
   integer(i4b)     :: D11gamma1       = imiss  ! 
   integer(i4b)     :: D21gamma1       = imiss  ! 
   integer(i4b)     :: D31gamma1       = imiss  ! 
   integer(i4b)     :: D41gamma1       = imiss  ! 
   integer(i4b)     :: exp1gamma1      = imiss  ! 
   integer(i4b)     :: exp1gamma2      = imiss  ! 
   integer(i4b)     :: bbl1gamma1      = imiss  ! 
   integer(i4b)     :: bbl1gamma2      = imiss  ! 
   integer(i4b)     :: BD1gamma1       = imiss  ! 
   integer(i4b)     :: SD1gamma1       = imiss  ! 
   integer(i4b)     :: WcrFrac1gamma1  = imiss  ! 
   integer(i4b)     :: WpwpFrac1gamma1 = imiss  ! 
   integer(i4b)     :: fsm1gamma1      = imiss  ! 
   integer(i4b)     :: zk1gamma1       = imiss  ! 
   integer(i4b)     :: zsk1gamma1      = imiss  ! 
   integer(i4b)     :: zsk1gamma2      = imiss  ! 
   integer(i4b)     :: zpk1gamma1      = imiss  ! 
   integer(i4b)     :: pfree1gamma1    = imiss  ! 
   integer(i4b)     :: rexp1gamma1     = imiss  ! 
   ! Beta parameter
   integer(i4b)     :: uhshape         = imiss  ! uh gamma pdf shape parameter
   integer(i4b)     :: uhscale         = imiss  ! uh gamma pdf scale parameter
   integer(i4b)     :: ks              = imiss  ! Saturated conductivity 
   integer(i4b)     :: bd              = imiss  ! bulk density
   integer(i4b)     :: sd              = imiss  ! soil density
   integer(i4b)     :: psis            = imiss  ! saturation matric potential 
   integer(i4b)     :: b               = imiss  ! retension courve slope in log space
   integer(i4b)     :: phi             = imiss  ! porosity 
   integer(i4b)     :: fc              = imiss  ! field capacity
   integer(i4b)     :: wp              = imiss  ! wilting point 
   integer(i4b)     :: myu             = imiss  ! 
   integer(i4b)     :: binfilt         = imiss  ! 
   integer(i4b)     :: D1              = imiss  ! 
   integer(i4b)     :: D4              = imiss  ! 
   integer(i4b)     :: D2              = imiss  ! 
   integer(i4b)     :: D3              = imiss  ! 
   integer(i4b)     :: c               = imiss  ! 
   integer(i4b)     :: Dsmax           = imiss  ! 
   integer(i4b)     :: Ds              = imiss  ! 
   integer(i4b)     :: Ws              = imiss  ! 
   integer(i4b)     :: expt            = imiss  ! 
   integer(i4b)     :: bbl             = imiss  ! 
   integer(i4b)     :: h1              = imiss  ! 
   integer(i4b)     :: h2              = imiss  ! 
   integer(i4b)     :: h3              = imiss  !
   integer(i4b)     :: h4              = imiss  ! 
   integer(i4b)     :: h5              = imiss  ! 
   integer(i4b)     :: WcrFrac         = imiss  ! 
   integer(i4b)     :: WpwpFrac        = imiss  ! 
   integer(i4b)     :: twm             = imiss  ! 
   integer(i4b)     :: fwm             = imiss  ! 
   integer(i4b)     :: fsm             = imiss  ! 
   integer(i4b)     :: fpm             = imiss  ! 
   integer(i4b)     :: zk              = imiss  ! 
   integer(i4b)     :: zsk             = imiss  ! 
   integer(i4b)     :: zpk             = imiss  ! 
   integer(i4b)     :: pfree           = imiss  ! 
   integer(i4b)     :: zperc           = imiss  ! 
   integer(i4b)     :: rexp            = imiss  ! 
   integer(i4b)     :: rmin            = imiss  ! minimum stomatal resistance
   integer(i4b)     :: lai             = imiss  ! Lai
 endtype iLook_par

! ***********************************************************************************************************
!  Define indices for variables for mapping data
! ***********************************************************************************************************
 type, public  ::  iLook_VarMapData
  integer(i4b)     :: hru_id         = imiss        ! hru id 
  integer(i4b)     :: weight         = imiss        ! Areal weight of overlapping polygon for hru 
  integer(i4b)     :: overlapPolyId  = imiss        ! Id of overlapping polygon for hru 
 endtype iLook_varMapData

! ***********************************************************************************************************
!  Define index for variables in Veg data 
! ***********************************************************************************************************
 type, public  ::  iLook_VarVegData
  integer(i4b)     :: polyid          = imiss       ! veg polygon id
  integer(i4b)     :: grnfrc          = imiss       ! greeness fraction 
  integer(i4b)     :: lai             = imiss       ! lai 
  integer(i4b)     :: vegclass        = imiss       ! veg class 
 endtype iLook_VarVegData

! ***********************************************************************************************************
!  Define index for variables in soil data 
! ***********************************************************************************************************
 type, public  ::  iLook_VarSoilData
  integer(i4b)     :: polyid       = imiss      ! soil polygon id
  integer(i4b)     :: hslyrs       = imiss      ! soil layer thickness 
  integer(i4b)     :: soilclass    = imiss      ! soil polygon id
  integer(i4b)     :: sand_frc     = imiss      ! sand fraction in soil polygon and layer
  integer(i4b)     :: silt_frc     = imiss      ! silt fraction in soil polygon and layer
  integer(i4b)     :: clay_frc     = imiss      ! clay fraction in soil polygon and layer
  integer(i4b)     :: bulk_density = imiss      ! bulk density in soil polygon and layer 
  integer(i4b)     :: ele_mean     = imiss      ! mean elevation over a soil polygon 
  integer(i4b)     :: ele_std      = imiss      ! standard deviation of elevation over a soil polygon 
  integer(i4b)     :: slp_mean     = imiss      ! mean slope over a soil polygon 
 endtype iLook_VarSoilData

! ***********************************************************************************************************
!  Define indices for Veg properties
! ***********************************************************************************************************
 type, public  ::  iLook_PrpVeg
  integer(i4b)     :: lai01          = imiss       ! Jan LAI  m^2/m^2
  integer(i4b)     :: lai02          = imiss       ! Feb LAI
  integer(i4b)     :: lai03          = imiss       ! Mar LAI
  integer(i4b)     :: lai04          = imiss       ! Apr LAI
  integer(i4b)     :: lai05          = imiss       ! May LAI
  integer(i4b)     :: lai06          = imiss       ! Jun LAI
  integer(i4b)     :: lai07          = imiss       ! Jul LAI
  integer(i4b)     :: lai08          = imiss       ! Aug LAI
  integer(i4b)     :: lai09          = imiss       ! Sep LAI
  integer(i4b)     :: lai10          = imiss       ! Oct LAI
  integer(i4b)     :: lai11          = imiss       ! Nov LAI
  integer(i4b)     :: lai12          = imiss       ! Dec LAI
  integer(i4b)     :: vegtype        = imiss       ! vege type
  integer(i4b)     :: nroot          = imiss       ! rooitng depth m 
  integer(i4b)     :: snup           = imiss       ! threshold SWE depth that implies 100% snow cover 
  integer(i4b)     :: rs             = imiss       ! stomatal resistance
  integer(i4b)     :: mrs            = imiss       ! minimum stomatal resistance
  integer(i4b)     :: leafDim        = imiss       ! characteristic leaf dimension
  integer(i4b)     :: can_top_h      = imiss       ! height of top of vegetation canopy above ground
  integer(i4b)     :: can_bot_h      = imiss       ! height of bottom of vegetation canopy above ground
  integer(i4b)     :: c_veg          = imiss       ! specific heat of vegetation J kg-1 K-1 
  integer(i4b)     :: maxMassVeg     = imiss       ! max. mass of vegetation with full foliage km m-2 
 endtype iLook_PrpVeg

! ***********************************************************************************************************
! define data vectors
! ***********************************************************************************************************
 type(iLook_par),        public,parameter  :: ixPar          = iLook_Par         (1,2,3,4,5,6,7,8,9,10,&
                                                                                 11,12,13,14,15,16,17,18,19,20,&
                                                                                 21,22,23,24,25,26,27,28,29,30,&
                                                                                 31,32,33,34,35,36,37,38,39,40,&
                                                                                 41,42,43,44,45,46,47,48,49,50,&
                                                                                 51,52,53,54,55,56,57,58,59,60,&
                                                                                 61,62,63,64,65,66,67,68,69,70,&
                                                                                 71,72,73,74,75,76,77,78,79,80.&
                                                                                 81,82,83,84,85,86,87,88,89,90)
 type(iLook_VarMapData),  public,parameter :: ixVarMapData   = iLook_VarMapData  (1,2,3)
 type(iLook_VarSoilData), public,parameter :: ixVarSoilData  = iLook_VarSoilData (1,2,3,4,5,6,7,8,9,10)
 type(iLook_VarVegData),  public,parameter :: ixVarVegData   = iLook_VarVegData  (1,2,3,4)
 type(iLook_PrpVeg),      public,parameter :: ixPrpVeg       = iLook_PrpVeg      (1,2,3,4,5,6,7,8,9,10,&
                                                                                  11,12,13,14,15,16,17,18,19,20,&
                                                                                  21,22)

! ***********************************************************************************************************
! define size of data vectors
! ***********************************************************************************************************
! Number of vairables defined
 integer(i4b),parameter,public    :: nPar = 90 
 integer(i4b),parameter,public    :: nVarMapData=3
 integer(i4b),parameter,public    :: nVarSoilData=10
 integer(i4b),parameter,public    :: nVarVegData=4
 integer(i4b),parameter,public    :: nPrpVeg=22

END MODULE var_lookup
