module read_topodata

USE nrtype
USE netcdf
USE data_type 
USE public_var                                   ! Including common constant (physical constant, other e.g., missingVal, etc.)
USE var_lookup,only:ixVarTopoData,nVarTopoData   ! index of topo polygon variables and number of variables

implicit none

private

public::getTopoData

contains

 ! *********************************************************************
 ! Subroutine: get topo data for gridbox and layers 
 ! *********************************************************************
  subroutine getTopoData(fname,       &   ! input: file name
                         tdata_meta,  &   ! input: topo data meta
                         dname_tpoly, &   ! input: dimension name for gridbox 
                         tdata,       &   ! input-output: topo data structure
                         nTpoly,      &   ! output: number of gridbox 
                         ierr, message)   ! output: error control
  implicit none
  ! input variables
  character(*),  intent(in)       :: fname          ! filename
  type(var_meta),intent(in)       :: tdata_meta(:)  ! topo data meta
  character(*),  intent(in)       :: dname_tpoly    ! dimension name for polygon
  ! input-output
  type(namevar), intent(inout)    :: tdata(:)       ! topo data container
  ! output variables
  integer(i4b),  intent(out)      :: nTpoly         ! number of gridbox 
  integer(i4b),  intent(out)      :: ierr           ! error code
  character(*),  intent(out)      :: message        ! error message
  ! local variables
  integer(i4b)                    :: iVar           ! variable index
  integer(i4b)                    :: ncid           ! NetCDF file ID
  integer(i4b)                    :: idimID_poly    ! dimension ID for HRUs
  integer(i4b)                    :: iVarID         ! variable ID
  ! initialize error control
  ierr=0; message='getTopoData/'
 
  ! open file for reading
  ierr = nf90_open(fname, nf90_nowrite, ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  ! get the ID of the poly dimension
  ierr = nf90_inq_dimid(ncid, dname_tpoly, idimID_poly)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_tpoly); return; endif
  ! get the length of the poly dimension
  ierr = nf90_inquire_dimension(ncid, idimID_poly, len=nTpoly)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_tpoly); return; endif
  ! ** read in Topo poly variables
  do iVar=1,size(tdata)
    ! get the variable ID
    ierr = nf90_inq_varid(ncid, trim(tdata_meta(ivar)%varName), iVarID)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(tdata_meta(ivar)%varName); return; endif
    select case(tdata_meta(iVar)%vartype)
      case('integer')
       select case(tdata_meta(iVar)%vardims)
         case('2D')
           ierr=20; message=trim(message)//'There is no 2D data for topographic data'; return
         case('1D')
           ! allocate space for the 1D integer array 
           allocate(tdata(ivar)%ivar1(nTpoly),stat=ierr)
           if(ierr/=0)then; message=trim(message)//'problem allocating 1D int space for tdata data structure'; return; endif
           ! get the data
           ierr = nf90_get_var(ncid, iVarID, tdata(ivar)%ivar1)
         end select
     case('double')
       select case(tdata_meta(iVar)%vardims)
         case('2D')
           ierr=20; message=trim(message)//'There is no 2D data for topographic data'; return
         case('1D')
           ! allocate space for the 1D integer array 
           allocate(tdata(ivar)%dvar1(nTpoly),stat=ierr)
           if(ierr/=0)then; message=trim(message)//'problem allocating 1D real space for tdata data structure'; return; endif
           ! get the data
           ierr = nf90_get_var(ncid, iVarID, tdata(ivar)%dvar1)
       end select
   end select 
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 
  end do  ! (looping through variables)
 
  ! close the NetCDF file
  ierr = nf90_close(ncid)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 
  end subroutine
 
end module read_topodata
