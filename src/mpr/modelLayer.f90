module modelLayer

! Compute model soil parameter transfer function 

USE nrtype                                        ! variable types, etc.
USE data_type                                     ! Including custum data structure definition
USE public_var                                    ! Including common constant (physical constant, other e.g., missingVal, etc.)
USE def_soiltfParm                                ! define/read transfer function parameters
USE var_lookup,only:ixSoilParVic,nSoilParVic      ! index of Model soil parameter variables and number of variables
USE var_lookup,only:ixSoilParSumma,nSoilParSumma  ! index of Model soil parameter variables and number of variables
USE var_lookup,only:ixPrpSoil,nPrpSoil            ! index of soil properties and number of properties

implicit none

!Following accessible outside this module
public::comp_model_depth
!anything else
private

contains

! *****************************************************************************
! Public subroutine: Execute computation of model soil layer depth/thickness 
! ******************************************************************************
subroutine comp_model_depth(ParSxyMz,    &        ! in/output: model soil depth for all soil polygons that are input
                            nMLyr,       &        ! input: number of model layers
                            nSPoly,      &        ! input: number of soil polygons
                            PrpSxySz,    &        ! input: soil property storage containing soil data layer depth and thickness
                            sclass,      &        ! input: soil class
                            ierr, message)        ! output: error control

  ! input
  integer(i4b),intent(in)     :: nSPoly           ! number of soil polygons
  integer(i4b),intent(in)     :: nMLyr            ! number of model layer
  type(namedvar2),intent(in)  :: PrpSxySz(:)      ! soil property storage containing soil layer depth and thickness 
  integer(i4b), intent(in)    :: sclass(:,:)      ! soil class
  ! Output
  type(dat_d2d),intent(inout) :: ParSxyMz(:)      ! storage of model soil parameter in model vertical domain ParSxyMz(:)%dat(mlyr,poly) 
  integer(i4b), intent(out)   :: ierr             ! error code
  character(*), intent(out)   :: message          ! error message
  ! Local 
  integer(i4b)                :: iPoly            ! Loop index of polygon 

  ! Initialize error control
  ierr=0; message='comp_model_depth/'

   do iPoly=1,nSPoly
     call comp_model_lyr(ParSxyMz(ixSoilParVic%h)%dat(:,iPoly), &
                         ParSxyMz(ixSoilParVic%z)%dat(:,iPoly), &
                         nMLyr,                                 &
                         PrpSxySz(ixPrpSoil%h)%varData(:,iPoly),&
                         a_d_vic,                               &
                         sclass(:,iPoly),                       &
                         dmiss,                                 &
                         ierr,message)
   enddo

end subroutine comp_model_depth 

! ***********
! subroutine to compute model layer thickness and bottom depth
! *********************************************************************
 subroutine comp_model_lyr(hLyrModel,   &      ! output: model layer thickness [mm]
                           zLyrModel,   &      ! output: model layer bottom depth [mm]
                           nMLyr,       &      ! input : number of model layer
                           hLyrSoil,    &      ! input : soil layer thickness [mm]
                           fracLyrCoef, &      ! input : fraction of total soil layer for each model layer
                           sclass,      &      ! input : usda soil class
                           missingVal,  &
                           ierr,message)       ! output: error control
 
 implicit none
 
 ! Define variables
 ! input
 integer(i4b),intent(in)     :: nMLyr            ! number of model layer
 real(dp),    intent(in)     :: hLyrSoil(:)      ! soil layer thickness [mm]
 real(dp),    intent(in)     :: fracLyrCoef(:)   ! fraction of total soil layer for each model layer
 integer(i4b),intent(in)     :: sclass(:)        ! usda soil class
 real(dp),    intent(in)     :: missingVal       
 ! output 
 real(dp), intent(inout)     :: hLyrModel(:)     ! model layer thickness [mm]
 real(dp), intent(inout)     :: zLyrModel(:)     ! depth of model layer bottom [mm]
 integer(i4b)                :: ierr             ! error code
 character(*),intent(inout)  :: message          ! error message
 ! local 
 integer(i4b)                :: iLyr             ! loop index of model layer 
 integer(i4b)                :: nSLyr            ! number of Soil layer (to be computed based on input soil data)
 integer(i4b)                :: nFrac            ! number of fraction coefficients 
 real(dp)                    :: topZ             ! depth of top of model layers
 real(dp)                    :: Ztot_in          ! total depth of soil layers
 logical(lgt),allocatable    :: mask(:)
 real(dp),allocatable        :: lyr_packed(:)
 integer(i4b)                :: nElm

 ! Initialize error control
 ierr=0; message=trim(message)//'comp_model_lyr/'

 ! Get number of soil layer 
 nSLyr=size(hLyrSoil) 

 ! Make mask to exclude layers with water/bedrock/other  
 allocate(mask(nSLyr),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for mask(nElm)';return;endif
 mask = (sclass < 13) !Exclude 14=water, 15=bedrock, 16=other(??)
 
 allocate(lyr_packed(count(mask)),stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for lyr_packed(nElm)';return;endif
 
 lyr_packed = pack(hLyrSoil,mask)
 nElm = size(lyr_packed) !number of valid soil layer
 
 if (nElm > 0) then ! if actually soil layers exist
   Ztot_in = sum(lyr_packed) ! Total depth for soil layer
   nFrac=size(fracLyrCoef)   ! number of layer fraction that is input
 
   if (nMLyr == 1) then
     if (nFrac /= 1) then
       ierr=10; message=trim(message)//'if model has single soil layer, nFrac has to be one'
       return
     endif
     hLyrModel(1)=Ztot_in
     zLyrModel(1)=hLyrModel(1)
   else  
     if (nFrac+1 /= nMLyr)then
       ierr=15; message=trim(message)//'number of nFrac does not match with number of model layer'
       return
     endif
     topZ=0.0_dp ! depth of model layer top (1st model layer = 0 m)
     do iLyr=1,nMLyr-1
       hLyrModel(iLyr)=fracLyrCoef(iLyr)*(Ztot_in-topZ)
       zLyrModel(iLyr)=topZ+hLyrModel(iLyr)
       topZ=zLyrModel(iLyr)   !update depth of layr top
     enddo
     hLyrModel(nMLyr)=Ztot_in-topZ
     zLyrModel(nMLyr)=topZ+hLyrModel(nMLyr)
   endif
 else ! if there are no real soil layers
   hLyrModel = missingVal
   zLyrModel = missingVal
 endif 
 
 end subroutine comp_model_lyr 

end module modelLayer 
