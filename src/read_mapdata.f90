module read_mapdata

use nrtype
use netcdf
use data_type 
use public_var 
use var_lookup,only:ixVarMapData,nVarMapData    ! index of variables in mapping data netCDF and number of variables

implicit none

private

public::getMapData

contains

! *********************************************************************
! public subroutine: get mapping data between geophysical polygon and model hru  
! *********************************************************************
subroutine getMapData(fname,          &   ! input: file name
                      mapdata_meta,   &   ! input: soil data meta
                      dname_hru,      &   ! input: dimension name for model hru polygon
                      dname_overPoly, &   ! input: dimension name for intersection geophysical data polygon 
                      mapdata,        &   ! input-output: soil data structure
                      nHruMap,        &   ! output: number of model hru polygons 
                      err, message)       ! output: error control
  implicit none
  ! input variables
  character(*), intent(in)           :: fname           ! filename
  type(var_meta),intent(in)          :: mapdata_meta(:) ! data meta
  character(*), intent(in)           :: dname_hru       ! dimension name for hru
  character(*), intent(in)           :: dname_overPoly  ! dimension name for overlap polygon 
  ! input-output
  type(mapvar),intent(inout)         :: mapdata(:)      ! map data container
  ! output variables
  integer(i4b), intent(out)          :: nHruMap         ! number of soil polygons
  integer(i4b), intent(out)          :: err             ! error code
  character(*), intent(out)          :: message         ! error message
  ! local variables
  integer(i4b)                       :: nOverPoly       ! dimension size for intersecting geophysical datat polygon 
  integer(i4b)                       :: iVar            ! variable index
  integer(i4b)                       :: ncid            ! NetCDF file ID
  integer(i4b)                       :: idimID_hru      ! dimension ID for stream segments
  integer(i4b)                       :: idimID_overPoly ! dimension ID for HRUs
  integer(i4b)                       :: iVarID          ! variable ID
 
  ! initialize error control
  err=0; message='getMapData/'
  ! open file for reading
  err = nf90_open(fname, nf90_nowrite, ncid)
  if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
  ! get the ID of the hru dimension
  err = nf90_inq_dimid(ncid, dname_hru, idimID_hru)
  if(err/=0)then; message=trim(message)//trim(nf90_strerror(err))//'; name='//trim(dname_hru); return; endif
  ! get the length of the hru dimension
  err = nf90_inquire_dimension(ncid, idimID_hru, len=nHruMap)
  if(err/=0)then; message=trim(message)//trim(nf90_strerror(err))//'; name='//trim(dname_hru); return; endif
  ! get the ID of the overlap polygon dimension
  err = nf90_inq_dimid(ncid, dname_overPoly, idimID_overPoly)
  if(err/=0)then; message=trim(message)//trim(nf90_strerror(err))//'; name='//trim(dname_overPoly); return; endif
  ! get the length of the overlap polygon dimension
  err = nf90_inquire_dimension(ncid, idimID_overPoly, len=nOverPoly)
  if(err/=0)then; message=trim(message)//trim(nf90_strerror(err))//'; name='//trim(dname_overPoly); return; endif
  allocate(mapdata(1)%var(nVarMapData),stat=err);
  if(err/=0)then; message=trim(message)//'problem allocating mapdata%var '; return; endif
  var:do iVar=1,size(mapdata_meta)
    ! get the variable ID
    err = nf90_inq_varid(ncid, trim(mapdata_meta(ivar)%varName), iVarID)
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err))//'; name='//trim(mapdata_meta(ivar)%varName); return; endif
    select case(mapdata_meta(iVar)%varName)
      case('hru_id') ! 1D integer array
        allocate(mapdata(1)%var(ivar)%ivar1(nHruMap),stat=err)
        if(err/=0)then; message=trim(message)//'problem allocating hru_id for mapdata data structure'; return; endif
        err = nf90_get_var(ncid, iVarID, mapdata(1)%var(ivar)%ivar1)
      case('weight') ! 1D double vector
        allocate(mapdata(1)%var(ivar)%dvar1(nOverPoly),stat=err)
        if(err/=0)then; message=trim(message)//'problem allocating weight for mapdata data structure'; return; endif
        err = nf90_get_var(ncid, iVarID, mapdata(1)%var(ivar)%dvar1)
      case('intersector') ! 1D integer vector 
        allocate(mapdata(1)%var(ivar)%ivar1(nOverPoly),stat=err)
        if(err/=0)then; message=trim(message)//'problem allocating overlapID for mapdata data structure'; return; endif
        err = nf90_get_var(ncid, iVarID, mapdata(1)%var(ivar)%ivar1)
      case('overlaps') ! 1D integer vector
        allocate(mapdata(1)%var(ivar)%ivar1(nHruMap),stat=err)
        if(err/=0)then; message=trim(message)//'problem allocating overlaps for mapdata data structure'; return; endif
        err = nf90_get_var(ncid, iVarID, mapdata(1)%var(ivar)%ivar1)
      case default;  err=10; message=trim(message)//'variableName_'//mapdata_meta(iVar)%varName//'_NotAvaiable'; return
    end select 
    if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
  enddo var
  err = nf90_close(ncid)
  if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
  return 
end subroutine

end module read_mapdata
