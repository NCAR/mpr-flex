module layerWeight 

USE nrtype                             ! variable types, etc.
USE data_type                          ! Including custum data structure definition
USE public_var                         ! Including public constant values 

implicit none

!Following accessible outside this module
public::map_slyr2mlyr
!anything else
private

contains

! *****************************************************************************
! Public subroutine: weight of soil layers within each mode layer 
! ******************************************************************************
 subroutine map_slyr2mlyr(zLyrSoil, hLyrSoil, zLyrModel, lyrmap, ierr, message)
  ! Define variables
  implicit none
  ! input
  real(dp), intent(in)        :: zLyrSoil(:)      ! depth of soil layer bottom [m]
  real(dp), intent(in)        :: hLyrSoil(:)      ! thickness of soil layer [m]
  real(dp), intent(in)        :: zLyrModel(:)     ! depth of model layer bottom [m]
  ! output 
  type(mapping),intent(inout) :: lyrmap(:)        ! data type storing weight and intersecting soil layer index for each model layer
  integer(i4b),intent(out)    :: ierr             ! error code
  character(*),intent(out)    :: message          ! error message
  ! local 
  integer(i4b)                :: ctr              ! counter 
  integer(i4b)                :: iSlyr            ! loop index of soil layer 
  integer(i4b)                :: iMlyr            ! loop index of model layer 
  integer(i4b)                :: nSLyr            ! number of Soil layer (to be computed based on input soil data)
  integer(i4b)                :: nMLyr            ! number of Soil layer (to be computed based on input soil data)
  real(dp),allocatable        :: Zs_top(:)        ! depth to top of ith soil layer (i=1..nSlyr)
  real(dp),allocatable        :: Zs_bot(:)        ! depth to bottom of ith soil layer (i=1..nSlyr)
  real(dp),allocatable        :: Zm_top(:)        ! depth to top of ith model layer (i=1..nMlyr)
  real(dp),allocatable        :: Zm_bot(:)        ! depth to bottom of ith model layer (i=1..nMlyr)
  integer(i4b),allocatable    :: idxTop(:)        ! index of soil layer of which top is within ith model layer (i=1..nMlyr)
  integer(i4b),allocatable    :: idxBot(:)        ! index of the lowest soil layer of which bottom is within ith model layer (i=1..nMlyr)
  
  ! initialize error control
  ierr=0; message='map_slyr2mlyr/'

  nSLyr=size(zLyrSoil)  !get soil layer number
  nMLyr=size(zLyrModel) !get model layer number
  
  !if model layer depth has missing value
  if ( zLyrModel(1) >= dmiss .and. zLyrModel(1) <= dmiss) then
    do iMLyr=1,nMLyr
      lyrmap(iMLyr)%weight  = dmiss
      lyrmap(iMLyr)%ixSubLyr= imiss
    enddo
  else
    allocate(Zs_bot(nSlyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating for Zs_bot'; return; endif
    allocate(Zs_top(nSlyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating for Zs_top'; return; endif
    allocate(Zm_top(nMLyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating for Zm_top'; return; endif
    allocate(Zm_bot(nMLyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating for Zm_bot'; return; endif
    allocate(idxTop(nMLyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating for idxTop'; return; endif
    allocate(idxBot(nMLyr),stat=ierr); if(ierr/=0)then; message=trim(message)//'error allocating for idxBot'; return; endif

    !-- Compute for depths to 1)top and 2) bottom of soil and model layer 
    Zm_top(1)=0.0_dp
    Zs_top(1)=0.0_dp
    do iMLyr=2,nMLyr
      Zm_top(iMLyr) = zLyrModel(iMLyr-1)
    end do
    do iSLyr=2,nSlyr
      Zs_top(iSLyr) = zLyrSoil(iSLyr-1)
    end do
    Zm_bot = zLyrModel(:)
    Zs_bot = zLyrSoil(:)

    !-- Find index of upper-most soil layer which gets within model layer (for each model layer)
    ! condition: from top to bottom of soil layer, 1st soil layer whose bottom gets below top of i-th model layer 
    do iMLyr=1,nMLyr   
      do iSLyr = 1,nSlyr  
        if ( Zm_top(iMLyr)-Zs_bot(iSLyr)<0.0_dp ) then 
          idxTop(iMLyr) = iSLyr; exit
        endif
      enddo 
    enddo 
    !-- Find index of lowest soil layer which get within model layer (for each layer)
    ! condition: from top to bottom of soil layer 1st soil layer whose top get deeper than bottom of i-th model layer 
    do iMLyr=1,nMLyr
      do iSLyr = 1,nSlyr
        if ( Zm_bot(iMLyr)-Zs_bot(iSLyr)<=valMin ) then 
            idxBot(iMLyr) = iSLyr; exit
        endif
      enddo 
    enddo 
    ! Error check
    do iMLyr=1,nMLyr
      if (idxTop(iMLyr)>11)then;              message=trim(message)//'index of idxTop not assinged'; return; endif 
      if (idxBot(iMLyr)>11)then;              message=trim(message)//'index of idxBot not assinged'; return; endif
      if (idxTop(iMLyr)-idxBot(iMLyr)>0)then; message=trim(message)//'index of idxTop lower than idxBot'; return; endif
    enddo

    !-- Compute weight of soil layer contributing to each model layer and populate lyrmap variable
    do iMLyr=1,nMLyr
      ctr = 1
      ! loop frm the upper most soil layer to the lowest soil layer, but only soil layers that intersect current model layer 
      do iSLyr=idxTop(iMlyr),idxBot(iMLyr)         
        if ( idxTop(iMlyr) == idxBot(iMlyr) ) then ! if model layer is completely within soil layer
            lyrmap(iMLyr)%weight(ctr)= 1.0 
            lyrmap(iMLyr)%ixSubLyr(ctr)= iSLyr
        else                                      ! if model layer contains multiple soil layers  
          if ( iSLyr == idxTop(iMLyr) ) then      ! for the upper most soil layer that intersect model layer 
            lyrmap(iMLyr)%weight(ctr)= (Zs_bot(iSLyr)-Zm_top(iMLyr))/zLyrModel(iMLyr)
            lyrmap(iMLyr)%ixSubLyr(ctr)= iSLyr 
          elseif ( iSLyr == idxBot(iMLyr) ) then  ! for the lowest soil layer that intersect model layer 
            lyrmap(iMLyr)%weight(ctr)= (Zm_bot(iMLyr)-Zs_top(iSLyr))/zLyrModel(iMLyr)
            lyrmap(iMLyr)%ixSubLyr(ctr)= iSLyr 
          else                                    ! for soil layers that completely in model layer 
            lyrmap(iMLyr)%weight(ctr)=hLyrSoil(iSLyr)/zLyrModel(iMLyr)
            lyrmap(iMLyr)%ixSubLyr(ctr)= iSLyr 
          endif
        endif
        ctr = ctr+1
      enddo 
    enddo
  endif 

 end subroutine map_slyr2mlyr 

end module layerWeight
