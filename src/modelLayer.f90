module modelLayer
! Compute model soil parameter transfer function 
USE nrtype                                        ! variable types, etc.
USE data_type                                     ! Including custum data structure definition
USE public_var                                    ! Including common constant (physical constant, other e.g., missingVal, etc.)

implicit none

private

public::comp_model_depth
public::map_slyr2mlyr

contains
! ***********
! public subroutine to compute model layer thickness and bottom depth
! *********************************************************************
  subroutine comp_model_depth(hModel,   &          ! output: model layer thickness [mm]
                              zModel,   &          ! output: model layer bottom depth [mm]
                              hfrac,    &          ! input : fraction of total soil layer for each model layer
                              soildata, &          ! input : soil data 
                              ierr, message)       ! output: error id and message 

  use var_lookup,  only:ixVarSoilData,nVarSoilData ! index of soil data variables and number of variables 
  implicit none
  ! input
  real(dp),     intent(in)    :: hfrac(:)          ! fraction of total soil layer for each model layer
  type(namevar),intent(in)    :: soilData(:)       ! soil data container for local soil polygon 
  ! output 
  real(dp),     intent(out)   :: hModel(:,:)       ! model layer thickness [mm]
  real(dp),     intent(out)   :: zModel(:,:)       ! depth of model layer bottom [mm]
  integer(i4b), intent(out)   :: ierr             ! error code
  character(*), intent(out)   :: message          ! error message
  ! local 
  integer(i4b)                :: iPoly             ! loop index of polygon 
  integer(i4b)                :: iLyr              ! loop index of model layer 
  integer(i4b)                :: nSLyr             ! number of Soil layer 
  integer(i4b)                :: nPoly             ! number of soil polygon 
  integer(i4b)                :: nFrac             ! number of fraction coefficients 
  real(dp)                    :: topZ              ! depth of top of model layers
  real(dp)                    :: Ztot_in           ! total depth of soil layers
  logical(lgc),allocatable    :: mask(:,:)
  real(dp),allocatable        :: lyr_packed(:)
  integer(i4b)                :: nElm
 
  ! Initialize error control
  ierr=0; message=trim(message)//'comp_model_depth/'
  ! Get number of soil layer 
  associate( hslyrs => soilData(ixVarSoilData%hslyrs)%dvar2 )
  nSLyr=size(hslyrs,1) 
  nPoly=size(hslyrs,2) 
  ! Make mask to exclude layers with water/bedrock/other  
  allocate(mask(nSLyr,nPoly),stat=ierr); if(ierr/=0)then;message=trim(message)//'error allocating mask';return;endif
  mask = (soilData(ixVarSoilData%soilclass)%ivar2 < 13) !Exclude 14=water, 15=bedrock, 16=other(??)
  do iPoly=1,nPoly
    allocate(lyr_packed(count(mask(:,iPoly))),stat=ierr); if(ierr/=0)then;message=trim(message)//'error allocating lyr_packed';return;endif
    lyr_packed = pack( hslyrs(:,iPoly), mask(:,iPoly) )
    nElm = size(lyr_packed)      ! number of valid soil layer
    if (nElm > 0) then           ! if actually soil layers exist
      Ztot_in = sum(lyr_packed)  ! Total depth for soil layer
      nFrac=size(hfrac)          ! number of layer fraction that is input
      if (nLyr == 1) then
        if (nFrac /= 1) then
          ierr=10;message=trim(message)//'if model has single soil layer, nFrac has to be one';return
        endif
        hModel(1,iPoly)=Ztot_in
        zModel(1,iPoly)=hModel(1,iPoly)
      else  
        if (nFrac+1 /= nLyr)then
          ierr=15;message=trim(message)//'number of nFrac does not match with number of model layer';return
        endif
        topZ=0.0_dp ! depth of model layer top (1st model layer = 0 m)
        do iLyr=1,nLyr-1
          hModel(iLyr,iPoly)=hfrac(iLyr)*(Ztot_in-topZ)
          zModel(iLyr,iPoly)=topZ+hModel(iLyr,iPoly)
          topZ=zModel(iLyr,iPoly)          !update depth of layr top
        enddo
        hModel(nLyr,iPoly)=Ztot_in-topZ
        zModel(nLyr,iPoly)=topZ+hModel(nLyr,iPoly)
      endif
    else ! if there are no real soil layers
      hModel(:,iPoly) = dmiss
      zModel(:,iPoly) = dmiss
    endif 
    deallocate(lyr_packed)
  enddo
  end associate
  return 
end subroutine 

! *****************************************************************************
! Public subroutine: weight of soil layers within each mode layer 
! ******************************************************************************
 subroutine map_slyr2mlyr( hSoil, zModel, lyrmap, ierr, message)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)        :: hSoil(:,:)      ! thickness of soil layer [m]
  real(dp), intent(in)        :: zModel(:,:)     ! depth of model layer bottom [m]
  ! output 
  type(poly),intent(inout)    :: lyrmap(:)        ! data type storing weight and intersecting soil layer index for each model layer
  integer(i4b),intent(out)    :: ierr             ! error code
  character(*),intent(out)    :: message          ! error message
  ! local 
  real(dp),    allocatable    :: zSoil(:,:)      ! thickness of soil layer [m]
  integer(i4b),parameter      :: nSub=11          ! max. number of Soil layer within Model layer
  integer(i4b)                :: ctr              ! counter 
  integer(i4b)                :: iPoly            ! loop index of polygon 
  integer(i4b)                :: iSlyr            ! loop index of soil layer 
  integer(i4b)                :: iMlyr            ! loop index of model layer 
  integer(i4b)                :: nPoly            ! number of polygon 
  integer(i4b)                :: nSLyr            ! number of Soil layer (to be computed based on input soil data)
  real(dp),allocatable        :: Zs_top(:)        ! depth to top of ith soil layer (i=1..nSlyr)
  real(dp),allocatable        :: Zs_bot(:)        ! depth to bottom of ith soil layer (i=1..nSlyr)
  real(dp),allocatable        :: Zm_top(:)        ! depth to top of ith model layer (i=1..nLyr)
  real(dp),allocatable        :: Zm_bot(:)        ! depth to bottom of ith model layer (i=1..nLyr)
  integer(i4b),allocatable    :: idxTop(:)        ! index of soil layer of which top is within ith model layer (i=1..nLyr)
  integer(i4b),allocatable    :: idxBot(:)        ! index of the lowest soil layer of which bottom is within ith model layer (i=1..nLyr)
  
  ! initialize error control
  ierr=0; message='map_slyr2mlyr/'
  ! dimensions
  nSLyr=size(hSoil,1)  !get soil layer number
  nPoly=size(hSoil,2)  !get polygon 
  if (nPoly /= size(zModel,2))then;ierr=30;message=trim(message)//'number of polygon mismatch'; return; endif
  allocate(zSoil,source=hSoil)
  do iSLyr=2,nSlyr
    zSoil(iSlyr,:)=hSoil(iSlyr,:)+zSoil(iSlyr-1,:)
  enddo
  do iPoly=1,nPoly
    do iMLyr=1,nLyr
      lyrmap(iPoly)%layer(iMLyr)%weight = dmiss
      lyrmap(iPoly)%layer(iMLyr)%ixSubLyr = imiss
    enddo
    if ( zModel(1,iPoly) >= dmiss .and. zModel(1,iPoly) <= dmiss) then !if model layer depth has missing value
      do iMLyr=1,nLyr
        lyrmap(iPoly)%layer(iMLyr)%weight  = dmiss
        lyrmap(iPoly)%layer(iMLyr)%ixSubLyr= imiss
      enddo
    else
      allocate(Zs_bot(nSlyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating Zs_bot'; return; endif
      allocate(Zs_top(nSlyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating Zs_top'; return; endif
      allocate(Zm_top(nLyr),stat=ierr);  if(ierr/=0)then; message=trim(message)//'error allocating Zm_top'; return; endif
      allocate(Zm_bot(nLyr),stat=ierr);  if(ierr/=0)then; message=trim(message)//'error allocating Zm_bot'; return; endif
      allocate(idxTop(nLyr),stat=ierr);  if(ierr/=0)then; message=trim(message)//'error allocating idxTop'; return; endif
      allocate(idxBot(nLyr),stat=ierr);  if(ierr/=0)then; message=trim(message)//'error allocating idxBot'; return; endif
      !-- Compute for depths to 1)top and 2) bottom of soil and model layer 
      Zm_top(1)=0.0_dp
      Zs_top(1)=0.0_dp
      do iMLyr=2,nLyr
        Zm_top(iMLyr) = zModel(iMLyr-1,iPoly)
      end do
      do iSLyr=2,nSlyr
        Zs_top(iSLyr) = zSoil(iSLyr-1,iPoly)
      end do
      Zm_bot = zModel(:,iPoly)
      Zs_bot = zSoil(:,iPoly)
      !-- Find index of upper-most soil layer which gets within model layer (for each model layer)
      ! condition: from top to bottom of soil layer, 1st soil layer whose bottom gets below top of i-th model layer 
      do iMLyr=1,nLyr   
        do iSLyr = 1,nSlyr  
          if ( Zm_top(iMlyr)-Zs_bot(iSLyr)<0.0_dp ) then 
            idxTop(iMlyr) = iSLyr; exit
          endif
        enddo 
      enddo 
      !-- Find index of lowest soil layer which get within model layer (for each layer)
      ! condition: from top to bottom of soil layer 1st soil layer whose top get deeper than bottom of i-th model layer 
      do iMLyr=1,nLyr
        do iSLyr = 1,nSlyr
          if ( Zm_bot(iMlyr)-Zs_bot(iSLyr)<=valMin ) then 
              idxBot(iMlyr) = iSLyr; exit
          endif
        enddo 
      enddo 
      ! Error check
      do iMLyr=1,nLyr
        if (idxTop(iMlyr)>11)then;             message=trim(message)//'index of idxTop not assinged'; return; endif 
        if (idxBot(iMlyr)>11)then;             message=trim(message)//'index of idxBot not assinged'; return; endif
        if (idxTop(iMlyr)-idxBot(iMlyr)>0)then; message=trim(message)//'index of idxTop lower than idxBot'; return; endif
      enddo
      !-- Compute weight of soil layer contributing to each model layer and populate lyrmap variable
      do iMLyr=1,nLyr
        ctr = 1
        ! loop frm the upper most soil layer to the lowest soil layer, but only soil layers that intersect current model layer 
        do iSLyr=idxTop(iMlyr),idxBot(iMLyr)         
          if ( idxTop(iMlyr) == idxBot(iMlyr) )then ! if model layer is completely within soil layer
              lyrmap(iPoly)%layer(iMLyr)%weight(ctr)   = 1.0 
              lyrmap(iPoly)%layer(iMLyr)%ixSubLyr(ctr) = iSLyr
          else                                      ! if model layer contains multiple soil layers  
            if ( iSLyr == idxTop(iMLyr) )then      ! for the upper most soil layer that intersect model layer 
              lyrmap(iPoly)%layer(iMLyr)%weight(ctr)   = (Zs_bot(iSLyr)-Zm_top(iMlyr))/zModel(iMlyr,iPoly)
              lyrmap(iPoly)%layer(iMLyr)%ixSubLyr(ctr) = iSLyr 
            elseif ( iSLyr == idxBot(iMLyr) ) then  ! for the lowest soil layer that intersect model layer 
              lyrmap(iPoly)%layer(iMLyr)%weight(ctr)   = (Zm_bot(iMlyr)-Zs_top(iSLyr))/zModel(iMLyr,iPoly)
              lyrmap(iPoly)%layer(iMLyr)%ixSubLyr(ctr) = iSLyr 
            else                                    ! for soil layers that completely in model layer 
              lyrmap(iPoly)%layer(iMLyr)%weight(ctr)   = hSoil(iSLyr,iPoly)/zModel(iMlyr,iPoly)
              lyrmap(iPoly)%layer(iMLyr)%ixSubLyr(ctr) = iSLyr 
            endif
          endif
          ctr = ctr+1
        enddo 
      enddo
      deallocate(Zs_bot,stat=ierr); if(ierr/=0)then; message=trim(message)//'error deallocating Zs_bot'; return; endif
      deallocate(Zs_top,stat=ierr); if(ierr/=0)then; message=trim(message)//'error deallocating Zs_top'; return; endif
      deallocate(Zm_top,stat=ierr); if(ierr/=0)then; message=trim(message)//'error deallocating Zm_top'; return; endif
      deallocate(Zm_bot,stat=ierr); if(ierr/=0)then; message=trim(message)//'error deallocating Zm_bot'; return; endif
      deallocate(idxTop,stat=ierr); if(ierr/=0)then; message=trim(message)//'error deallocating idxTop'; return; endif
      deallocate(idxBot,stat=ierr); if(ierr/=0)then; message=trim(message)//'error deallocating idxBot'; return; endif
    endif 
  enddo 
  return
 end subroutine

end module modelLayer 
