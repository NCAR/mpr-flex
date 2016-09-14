module pedotf

! Compute model soil parameter transfer function 

USE nrtype                                        ! variable types, etc.
USE data_type                                     ! Including custum data structure definition
USE public_var                                    ! Including common constant (physical constant, other e.g., dmiss, etc.)
USE def_pedotfParm                                ! define/read pedo transfer function parameters
USE var_lookup,only:ixPrpSoil,nPrpSoil            ! index of soil properties and number of properties
USE var_lookup,only:ixVarSoilData,nVarSoilData    ! index of soil data variables and number of properties

implicit none

!Following accessible outside this module
public::comp_pedotf
!anything else
private

contains

! ********************************************************************************************
! Public subroutine: Execute computation of soil properties with pedotransfer function  
! *********************************************************************************************
! Soil properties to be computed
!  --------------------------------------------------------------------------------------------
!  Property Name                          Notation    opt-1      opt-2        
!  --------------------------------------------------------------------------------------------
!  1. saturated hydraulic conductivity    ks          cosby      campbell & shiozawa 
!  2. bulk density                        db          Frm data   *** 
!  3. porosity                            phi         cosby      Aacharias & Wessolek 
!  4. Field capacity                      fc          campbel    ***
!  5. Wilting point                       wp          campbel    ***
!  6. Rtention curve slope in log space   b           cosby      ***
!  7. saturation matric potentioal        psi_sat     cosby      ***
!  8. specific yield                      myu         koren      ***
! ---------------------------------------------------------------------------------------------
!  *** indicates there is no equations avaialble 

subroutine comp_pedotf(sdata,     &      ! input:     soil data for all soil polygons  and layers
                       sprp,      &      ! in/output: soil properties for all soil polygons and layers
                       ierr, message)    ! output: err id and message 

  implicit none

  ! input
  type(namevar), intent(in)    :: sdata(:)             ! storage of soil data strucuture
  ! input/output
  type(namedvar2),intent(inout):: sprp(:)              ! soil property values for polygon and layers 
  ! Output
  integer(i4b),intent(out)     :: ierr                 ! error code
  character(*),intent(out)     :: message              ! error message
  ! Local 
  integer(i4b)                 :: iPrpSoil             ! Loop index of soil property
  real(dp),pointer             :: sand(:,:)            ! pointer to sand fraction array 
  real(dp),pointer             :: silt(:,:)            ! pointer to silt fraction array 
  real(dp),pointer             :: clay(:,:)            ! pointer to clay fraction array 
  real(dp),pointer             :: bd(:,:)              ! pointer to bulk density array 
  logical(lgt),allocatable     :: checkDone(:)         ! used to check if the soil property is processed
  
  ierr=0; message="comp_pedotf/"

  allocate(checkDone(nPrpSoil),stat=ierr); if(ierr/=0) then; message=trim(message)//'error allocating for checkDone'; return; endif
  checkDone(:) = .false.
  sand => sdata(ixVarSoilData%sand)%dvar2 
  silt => sdata(ixVarSoilData%silt)%dvar2 
  clay => sdata(ixVarSoilData%clay)%dvar2 
  bd   => sdata(ixVarSoilData%bd)%dvar2 
  do iPrpSoil = 1,nPrpSoil
    select case(iPrpSoil)
      case(ixPrpSoil%ks);
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = ks( sand, clay, ptfopt(1), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%bd);
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = bd( bd, ptfopt(2), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%phi);
        if(.not.checkDone(ixPrpSoil%bd))then; ierr=20; message=trim(message)//'need to process bd before phi'; return; endif
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = phi( sand, clay, bd, ptfopt(3), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%b);
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = ret_curve( sand, clay, ptfopt(6), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%psis);
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = psis( sand, silt, ptfopt(7), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%fc);
        if(.not.checkDone(ixPrpSoil%psis))then; ierr=20; message=trim(message)//'need to process psis before fc'; return; endif
        if(.not.checkDone(ixPrpSoil%phi))then;  ierr=20; message=trim(message)//'need to process phi before fc'; return; endif
        if(.not.checkDone(ixPrpSoil%b))then;    ierr=20; message=trim(message)//'need to process b before fc'; return; endif
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = fc(sand, sprp(ixPrpSoil%phi)%varData, sprp(ixPrpSoil%psis)%varData, sprp(ixPrpSoil%b)%varData,  ptfopt(4), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%wp);
        if(.not.checkDone(ixPrpSoil%psis))then; ierr=20; message=trim(message)//'need to process psis before wp'; return; endif
        if(.not.checkDone(ixPrpSoil%phi))then;  ierr=20; message=trim(message)//'need to process phi before wp'; return; endif
        if(.not.checkDone(ixPrpSoil%b))then;    ierr=20; message=trim(message)//'need to process b before wp'; return; endif
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = wp( sprp(ixPrpSoil%phi)%varData, sprp(ixPrpSoil%psis)%varData, sprp(ixPrpSoil%b)%varData, ptfopt(5), ierr, message)
        if(ierr/=0)then; return;endif
      case(ixPrpSoil%myu);
        checkDone(iPrpSoil)=.true.
        call sprp(iPrpSoil)%varData = myu( sprp(ixPrpSoil%phi)%varData, sprp(ixPrpSoil%fc)%varData, ptfopt(8), ierr, message)
        if(ierr/=0)then; return;endif
     end select ! end of parameter case
  end do ! end of parameter loop

  deallocate(checkDone,stat=ierr); if(ierr/=0) then; message=trim(message)//'error deallocating for checkDone(nVegParVic)'; return; endif
  return
end subroutine comp_pedotf

! Private subroutine
! *********************************************************************
! pedo-transfer function for saturated hydraulic conductivity (ks)
! *********************************************************************
  function ks( sand_in, clay_in, opt, ierr, message)
    implicit none
    ! input
    real(dp), intent(in)       :: sand_in(:,:) ! input: sand [percent] 
    real(dp), intent(in)       :: clay_in(:,:) ! input: clay [percent] 
    integer(i2b)               :: opt          ! input: option for transfer function form
    ! inout
    character(*),intent(inout) :: message      ! error message
    integer(i4b),intent(inout) :: ierr         ! error code
    ! local 
    real(dp)                   :: ks(:,:)      ! ks 

    ! opt 1: Cosby et al. WRR 1984
    ! opt 2: campbell & shiozawa 1994 
    message=trim(message)//"ks/"

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
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function ks

! *********************************************************************
! pedo-transfer function for bulk density 
! *********************************************************************
! There is bulk density from look up table [g/cm3] and soil data [kg/m3]
  function bd(Db_in, opt, ierr, message)
    implicit none
    ! input
    real(dp), intent(in)       :: Db_in(:,:)  ! input: bulk density [kg/m3]
    integer(i2b)               :: opt         ! input: option for transfer function form
    ! inout
    character(*),intent(inout) :: message     ! error message
    integer(i4b),intent(inout) :: ierr        ! error code
    ! local 
    real(dp)                   :: bd(:,:)     ! computed bulk density  

    ! opt 1: use soil data
    message=trim(message)//"bd/"

    select case(opt)
      case(1); 
        where ( Db_in /= dmiss ) 
          bd = gamma1db*Db_in 
        else where
          bd = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function bd 

! *********************************************************************
! pedo-transfer function for porosity 
! *********************************************************************
  function phi(sand_in, clay_in, db_in, opt, ierr, message)      
    implicit none
    ! input
    real(dp), intent(in)       :: sand_in(:,:)   ! input: sand  
    real(dp), intent(in)       :: clay_in(:,:)   ! input: clay 
    real(dp), intent(in)       :: db_in(:,:)     ! input: bulk density 
    integer(i2b)               :: opt            ! option for transfer function form
    ! inout
    character(*),intent(inout) :: message        ! error message
    integer(i4b),intent(inout) :: ierr           ! error code
    ! local 
    real(dp)                   :: phi(:,:)  ! estimated porosity [fraction]

    ! opt 1: Cosby et al. WRR 1984
    ! opt 2: Zacharias & Wessolek 2007
    message=trim(message)//"comp_porosity/"

    select case(opt)
      case(1);  ! Cosby
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          phi = gamma1phi+gamma2phi*sand_in+gamma3phi*clay_in
        else where
          phi = dmiss 
        end where
      case(2);  ! Zacharias & Wessolek 2007  
        where ( sand_in /= dmiss .and. clay_in /= dmiss .and. Db_in /= dmiss ) 
          where ( sand_in < 66.5_dp) 
            phi = gamma1phi+gamma2phi*clay_in+gamma3phi*Db_in/1000._dp
          else where
            phi = gamma4phi+gamma5phi*clay_in+gamma6phi*Db_in/1000._dp
          end where 
        else where
          phi = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function phi 

! *********************************************************************
! pedo-transfer function for field capacity 
! *********************************************************************
  function fc(sand_in, phi_in, psis_in, b_in, opt, ierr, message)
    implicit none
    ! input
    real(dp), intent(in)       :: sand_in(:,:)   ! input: sand  
    real(dp), intent(in)       :: phi_in(:,:)    ! input: porosity [fraction]  
    real(dp), intent(in)       :: psis_in(:,:)   ! input: saturation matric potential [ 
    real(dp), intent(in)       :: b_in(:,:)      ! input: slope of cambell retention curve 
    integer(i2b)               :: opt            ! id for transfer function form
    ! output 
    ! inout
    character(*),intent(inout) :: message        ! error message
    integer(i4b),intent(inout) :: ierr           ! error code
    ! local 
    real(dp)                   :: fc(:,:)        ! estimated field capacity [fraction]
    real(dp),allocatable       :: psi_fc(:,:)    ! matric potential at field capacity [kPa]  
    integer(i4b)               :: nSpoly         ! number of soil polygon 
    integer(i4b)               :: nSlyr          ! number of soil layer 

    ! opt 1: Campbell 1974 
    message=trim(message)//"fc/"

    nSpoly=size(porosity,1)
    nSlyr=size(porosity,2)
    allocate(psi_fc(nSpoly,nSlyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating space for psi_fc'; return; endif
    psi_fc(:,:)=-20
    where (sand_in > 69) psi_fc=-10

    select case(opt)
      case(1);  !campbell
        where ( porosity /= dmiss .and. sand_in /= dmiss .and. b /= dmiss .and. psi_sat /= dmiss ) 
          fc = gamma1fc*porosity*(psi_fc/psi_sat)**(-1/b)
        else where
          fc = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function fc

! *********************************************************************
! pedo-transfer function for wilting point 
! *********************************************************************
  function wp( phi_in, psis_in, b_in, opt, ierr, message) 
    implicit none
    ! input
    real(dp), intent(in)       :: phi_in(:,:)    ! input: porosity [fraction]  
    real(dp), intent(in)       :: psis_in(:,:)   ! input: saturation matric potential [kPa]  
    real(dp), intent(in)       :: b_in(:,:)      ! input: slope of cambell retention curve [-]
    integer(i2b)               :: opt            ! input: option for transfer function form
    ! inout
    character(*),intent(inout) :: message        ! error message
    integer(i4b),intent(inout) :: ierr           ! error code
    ! local 
    real(dp)                   :: wp(:,:)        ! estimated field capacity [frac]
    real(dp),allocatable       :: psi_wp(:,:)    ! matric potential at wilting point [kPa]  
    integer(i4b)               :: nSpoly         ! number of soil polygon 
    integer(i4b)               :: nSlyr          ! number of soil layer 

    ! opt 1: Campbell 1974 
    message=trim(message)//"wp/"

    nSpoly=size(porosity,1)
    nSlyr=size(porosity,2)
    allocate(psi_wp(nSpoly,nSlyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'problem allocating space for psi_wp';return;endif
    psi_wp(:,:)=-1500

    select case(opt)
      case(1);  !Cosby et al. 
        where ( porosity /= dmiss .and. b /= dmiss .and. psi_sat /= dmiss ) 
          wp = gamma1wp*porosity*(psi_wp/psi_sat)**(-1/b)
        else where
          wp = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
  end function wp

! *********************************************************************
! pedo-transfer function for b (slope of retention curve in log space)
! *********************************************************************
  subroutine ret_curve(sand_in, clay_in, opt, ierr, message)
    implicit none
    ! input
    real(dp), intent(in)       :: sand_in(:,:)  ! input: sand  
    real(dp), intent(in)       :: clay_in(:,:)  ! input: clay 
    integer(i2b)               :: opt           ! input: option for transfer function form
    ! inout
    character(*),intent(inout) :: message       ! error message
    integer(i4b),intent(inout) :: ierr          ! error code
    ! local 
    real(dp)                   :: ret_curv(:,:)        ! computed [-] 

    ! opt 1: Cosby et al. WRR 1984
    message=trim(message)//"ret_curv/"

    select case(opt)
      case(1); 
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          ret_curv = gamma1b+gamma2b*sand_in+gamma3b*clay_in
        else where
          ret_curv = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function ret_curve

! *********************************************************************
! pedo-transfer function for saturation matric potential 
! *********************************************************************
  function psis( sand_in, silt_in, opt, ierr, message)
    implicit none
    ! input
    real(dp), intent(in)       :: sand_in(:,:)   ! input: sand percentage [percentage]
    real(dp), intent(in)       :: silt_in(:,:)   ! input: silt percentage [percentage]
    integer(i2b)               :: opt            ! input: option for transfer function form
    ! inout
    character(*),intent(inout) :: message        ! error message
    integer(i4b),intent(inout) :: ierr           ! error code
    ! local 
    real(dp)                   :: psis(:,:)      ! output: saturation matric potential [kPa]  

    ! opt 1: Cosby et al. WRR 1984
    message=trim(message)//"psis/"

    select case(opt)
      case(1);  !Cosby et al. 
        where ( sand_in /= dmiss .and. silt_in /= dmiss ) 
          psis = gamma1psis + gamma2psis*sand_in + gamma3psis*silt_in
          psis = -1.0_dp*10**psi_sat*0.0980665_dp        ! 0.0980665 kPa/cm-H2O. Cosby give psi_sat in cm of water (cm-H2O)
        else where
          psis = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function psis

! *********************************************************************
! pedo-transfer function for specific yield  
! *********************************************************************
  function myu(phi_in, fc_in, opt, ierr, message)
    implicit none
    ! input
    real(dp), intent(in)       :: phi_in(:,:)       ! input: porosity [frac]
    real(dp), intent(in)       :: fc_in(:,:) ! input: field capacity [frac]
    integer(i2b)               :: opt                 ! input: option for transfer function form
    ! inout
    character(*),intent(inout) :: message             ! error message
    integer(i4b),intent(inout) :: ierr                ! error code
    ! local 
    real(dp)                   :: myu(:,:)            ! output: specific yield [-]  

    ! opt 1: Koren et al. 2003
    message=trim(message)//"myu/"

    select case(opt)
      case(1);  ! koren
        where ( porosity /= dmiss .and. field_capacity /= dmiss ) 
          myu = gamma1myu*(porosity-field_capacity)**gamma2myu 
        else where
          myu = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized'; ierr=35; return
    end select
  end function myu

end module pedotf 
