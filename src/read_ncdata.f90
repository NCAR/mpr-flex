module read_ncdata

! Collection of netCDF reading routines

USE netcdf
USE nrtype

implicit none

public::get_array2_dim
public::get_array2_dvar
public::get_array2_ivar
public::get_vec_dim
public::get_vec_ivar
public::get_scl_ivar
public::get_vec_dvar
public::get_scl_dvar

private

contains

! *********************************************************************
! subroutine: get vector dimension from netCDF 
! *********************************************************************
subroutine get_vec_dim(fname,           &  ! input: filename
                       dname,           &  ! input: variable name
                       nDim,            &  ! output: Size of dimension 
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: dname        ! dimension name
 ! output variables
 integer(i4b), intent(out)       :: nDim         ! size of dimension
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iDimID       ! NetCDF dimension ID

 ! initialize error control
 ierr=0; message='get_vec_dim/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the ID of the dimension
 ierr = nf90_inq_dimid(ncid, dname, iDimID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname); return; endif
 ! get the length of the dimension
 ierr = nf90_inquire_dimension(ncid, iDimID, len=nDim)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_vec_dim

! *********************************************************************
! subroutine: get 2D array dimension from netCDF 
! *********************************************************************
subroutine get_array2_dim(fname,         &  ! input: filename
                       dname1,           &  ! input: 1st dimension name
                       dname2,           &  ! input: 2nd dimension name
                       nDim1,            &  ! output: Size of 1st dimension 
                       nDim2,            &  ! output: Size of 2nd dimension 
                       ierr, message)       ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: dname1       ! 1st dimension name
 character(*), intent(in)        :: dname2       ! 2nd dimension name
 ! output variables
 integer(i4b), intent(out)       :: nDim1        ! size of 1st dimension
 integer(i4b), intent(out)       :: nDim2        ! size of 2nd dimension
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iDim1ID      ! NetCDF 1st dimension ID
 integer(i4b)                    :: iDim2ID      ! NetCDF 2nd dimension ID

 ! initialize error control
 ierr=0; message='get_array2_dim/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the ID of the 1st dimension
 ierr = nf90_inq_dimid(ncid, dname1, iDim1ID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname1); return; endif
 ! get the length of the 1st dimension
 ierr = nf90_inquire_dimension(ncid, iDim1ID, len=nDim1)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the ID of the 2nd dimension
 ierr = nf90_inq_dimid(ncid, dname2, iDim2ID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname1); return; endif
 ! get the length of the 2nd dimension
 ierr = nf90_inquire_dimension(ncid, iDim2ID, len=nDim2)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_array2_dim

! *********************************************************************
! subroutine: get integer scaler from netCDF
! *********************************************************************
subroutine get_scl_ivar(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        iStart,          &  ! input: start index
                        iScl,            &  ! output: outputvariable data
                        ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart    ! start index
 ! output variables
 integer(i4b), intent(out)                          :: iScl      ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b),dimension(1)                          :: iDummy    ! temporary vector of length 1
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_scl_ivar/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the data
 ierr = nf90_get_var(ncid, ivarID, iDummy, start=(/iStart/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! save the output
 iScl = iDummy(1)
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_scl_ivar

! *********************************************************************
! subroutine: get double precision scaler value from netCDF
! *********************************************************************
subroutine get_scl_dvar(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        iStart,          &  ! input: start index
                        dScl,            &  ! output: outputvariable data
                        ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)               :: fname     ! filename
 character(*), intent(in)               :: vname     ! variable name
 integer(i4b), intent(in)               :: iStart    ! start index
 ! output variables
 real(dp), intent(out)                  :: dScl      ! output variable data
 integer(i4b), intent(out)              :: ierr      ! error code
 character(*), intent(out)              :: message   ! error message
 ! local variables
 real(dp),dimension(1)                  :: dDummy    ! temporary vector of length 1
 integer(i4b)                           :: ncid      ! NetCDF file ID
 integer(i4b)                           :: iVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_scl_dvar/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the data
 ierr = nf90_get_var(ncid, ivarID, dDummy, start=(/iStart/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! save the output
 dScl = dDummy(1)
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_scl_dvar

! *********************************************************************
! subroutine: get integer vector value from netCDF
! *********************************************************************
subroutine get_vec_ivar(fname,           &  ! input: filename
                        vname,           &  ! input: variable name
                        iStart,          &  ! input: start index
                        iCount,          &  ! input: length of vector
                        iVec,            &  ! output: outputvariable data
                        ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart    ! start index
 integer(i4b), intent(in)                           :: iCount    ! length of vector to be read in
 ! output variables
 integer(i4b), intent(out),dimension(:)             :: iVec      ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_vec_ivar/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the data
 ierr = nf90_get_var(ncid, ivarID, iVec, start=(/iStart/), count=(/iCount/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

end subroutine get_vec_ivar

 ! *********************************************************************
 ! subroutine: read a double precision vector
 ! *********************************************************************
 subroutine get_vec_dvar(fname,           &  ! input: filename
                         vname,           &  ! input: variable name
                         dVec,            &  ! input: variable data
                         iStart,          &  ! input: start index
                         iCount,          &  ! input: length of vector
                         ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart    ! start index
 integer(i4b), intent(in)                           :: iCount    ! length of vector to be read in
 ! output variables
 !real(dp), intent(out), dimension(:), allocatable   :: dVec      ! output variable data
 real(dp), intent(out), dimension(:)                :: dVec      ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_vec_dVar/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! allocate space for the output 
! allocate(dVec(iCount),stat=ierr)
! if(ierr/=0)then; message=trim(message)//'problem allocating space for dVec'; return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the data
 ierr = nf90_get_var(ncid, ivarID, dVec, start=(/iStart/), count=(/iCount/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_vec_dvar

 ! *********************************************************************
 ! subroutine: read a double precision 2D array
 ! *********************************************************************
 subroutine get_array2_dvar(fname,           &  ! input: filename
                            vname,           &  ! input: variable name
                            iStart,          &  ! input: start row index (vertical direction) 
                            iCount,          &  ! input: length of row 
                            jStart,          &  ! input: start index of column (horizontal direction)
                            jCount,          &  ! input: length of column 
                            dArray,          &  ! output: variable data
                            ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart    ! start columb index (horizontal direction)
 integer(i4b), intent(in)                           :: iCount    ! length of columbs to be read in
 integer(i4b), intent(in)                           :: jStart    ! start row index (vertical direction)
 integer(i4b), intent(in)                           :: jCount    ! length of rows to be read in
 ! output variables
 real(dp), intent(out), dimension(:,:)              :: dArray    ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: dVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_array2_dVar/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),dVarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the data
 ierr = nf90_get_var(ncid, dVarID, dArray, start=(/jStart,iStart/), count=(/jCount,iCount/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_array2_dvar

 ! *********************************************************************
 ! subroutine: read a integer 2D array
 ! *********************************************************************
 subroutine get_array2_ivar(fname,           &  ! input: filename
                            vname,           &  ! input: variable name
                            iStart,          &  ! input: start row index (vertical direction) 
                            iCount,          &  ! input: length of row 
                            jStart,          &  ! input: start index of column (horizontal direction)
                            jCount,          &  ! input: length of column 
                            iArray,          &  ! output: variable data
                            ierr, message)      ! output: error control

 implicit none
 ! input variables
 character(*), intent(in)                           :: fname     ! filename
 character(*), intent(in)                           :: vname     ! variable name
 integer(i4b), intent(in)                           :: iStart    ! start columb index (horizontal direction)
 integer(i4b), intent(in)                           :: iCount    ! length of columbs to be read in
 integer(i4b), intent(in)                           :: jStart    ! start row index (vertical direction)
 integer(i4b), intent(in)                           :: jCount    ! length of rows to be read in
 ! output variables
 integer(i4b), intent(out), dimension(:,:)          :: iArray    ! output variable data
 integer(i4b), intent(out)                          :: ierr      ! error code
 character(*), intent(out)                          :: message   ! error message
 ! local variables
 integer(i4b)                                       :: ncid      ! NetCDF file ID
 integer(i4b)                                       :: iVarId    ! NetCDF variable ID

 ! initialize error control
 ierr=0; message='get_array2_iVar/'
 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the data
 ierr = nf90_get_var(ncid, iVarID, iArray, start=(/jStart,iStart/), count=(/jCount,iCount/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_array2_ivar

end module read_ncdata
