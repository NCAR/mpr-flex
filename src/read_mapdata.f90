module read_mapdata

USE nrtype
USE netcdf
USE data_type 
USE public_var 
USE var_lookup,only:ixVarMapData,nVarMapData    ! index of variables in mapping data netCDF and number of variables

implicit none

private

public::getMapData

contains

! *********************************************************************
! Subroutine: get mapping data between geophysical polygon and model hru  
! *********************************************************************
 subroutine getMapData(fname,          &   ! input: file name
                       vartype,        &   ! input: geophysical data type 
                       mapdata_meta,   &   ! input: soil data meta
                       dname_hru,      &   ! input: dimension name for soil polygon
                       dname_overPoly, &   ! input: dimension name for soil layer 
                       mapdata,        &   ! input-output: soil data structure
                       nHru,           &   ! output: number of soil polygon
                       nOverPoly,      &   ! output: number of soil layer  
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)           :: fname           ! filename
 character(*), intent(in)           :: vartype         ! dimension name for hru
 type(var_meta),intent(in)          :: mapdata_meta(:) ! data meta
 character(*), intent(in)           :: dname_hru       ! dimension name for hru
 character(*), intent(in)           :: dname_overPoly  ! dimension name for overlap polygon 
 ! input-output
 type(mapvar),intent(inout)         :: mapdata(:)      ! map data container
 ! output variables
 integer(i4b), intent(out)          :: nHru            ! number of soil polygons
 integer(i4b), intent(out)          :: nOverPoly       ! number of soil layers 
 integer(i4b), intent(out)          :: ierr            ! error code
 character(*), intent(out)          :: message         ! error message
 ! local variables
 integer(i4b)                       :: ix              ! geophysical data index
 integer(i4b)                       :: iVar            ! variable index
 integer(i4b)                       :: ncid            ! NetCDF file ID
 integer(i4b)                       :: idimID_hru      ! dimension ID for stream segments
 integer(i4b)                       :: idimID_overPoly ! dimension ID for HRUs
 integer(i4b)                       :: iVarID          ! variable ID
 ! initialize error control
 ierr=0; message='getMapData/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the ID of the hru dimension
 ierr = nf90_inq_dimid(ncid, dname_hru, idimID_hru)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_hru); return; endif
 ! get the length of the hru dimension
 ierr = nf90_inquire_dimension(ncid, idimID_hru, len=nHru)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_hru); return; endif
 ! get the ID of the overlap polygon dimension
 ierr = nf90_inq_dimid(ncid, dname_overPoly, idimID_overPoly)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_overPoly); return; endif
 ! get the length of the overlap polygon dimension
 ierr = nf90_inquire_dimension(ncid, idimID_overPoly, len=nOverPoly)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_overPoly); return; endif

 select case(trim(vartype))
   case('soil'); ix =1
   case('veg');  ix =2 
 end select
 allocate(mapdata(ix)%var(nVarMapData),stat=ierr);
 if(ierr/=0)then; message=trim(message)//'problem allocating mapdata(ix)%var '; return; endif

 ! go through variables specified in mapdata meta file to read them
 do iVar=1,size(mapdata_meta)
   ! get the variable ID
   ierr = nf90_inq_varid(ncid, trim(mapdata_meta(ivar)%varName), ivarID)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(mapdata_meta(ivar)%varName); return; endif
   select case(mapdata_meta(iVar)%varName)
     case('hru_id')
       ! allocate space for the 2D integer array 
       allocate(mapdata(ix)%var(ivar)%ivar1(nHru),stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating hru_id space for mapdata data structure'; return; endif
       ! get the data
       ierr = nf90_get_var(ncid, ivarID, mapdata(ix)%var(ivar)%ivar1)
     case('weight')
       ! allocate space for the 1D integer array 
       allocate(mapdata(ix)%var(ivar)%dvar2(nOverPoly,nHru),stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating weight space for mapdata data structure'; return; endif
       ! get the data
       ierr = nf90_get_var(ncid, ivarID, mapdata(ix)%var(ivar)%dvar2)
     case('overlapPolyId')
       ! allocate space for the 2D integer array 
       allocate(mapdata(ix)%var(ivar)%ivar2(nOverPoly,nHru),stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating overlapID space for mapdata data structure'; return; endif
       ! get the data
       ierr = nf90_get_var(ncid, ivarID, mapdata(ix)%var(ivar)%ivar2)
   end select 
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end do  ! (looping through variables)

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine getMapData

end module read_mapdata
