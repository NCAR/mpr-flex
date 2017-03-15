module read_vegdata

USE netcdf
USE nrtype                                   ! Including numerical type definition 
USE data_type                                ! Including custum data structure definition
USE public_var                               ! Including common constant (physical constant, other e.g., missingVal, etc.)
USE ascii_util,only:file_open                !
USE ascii_util,only:split_line               !
USE ascii_util,only:get_vlines               !
USE var_lookup,only:ixVarVegData,nVarVegData ! index of veg properties and number of properties
USE var_lookup,only:ixPrpVeg,nPrpVeg         ! index of veg properties and number of properties

private

public::getVegData
public::getvegClassLookup

contains

! *********************************************************************
! Public subroutine: get veg data for polygons (and months) 
! *********************************************************************
 subroutine getVegData(fname,       &   ! input: file name
                       vdata_meta,  &   ! input: veg data meta
                       dname_vpoly, &   ! input: dimension name for veg polygon
                       vdata,       &   ! input-output: veg data structure
                       nVpoly,      &   ! output: number of veg polygon
                       ierr, message)   ! output: error control
 
   implicit none
  
   character(*), intent(in)        :: fname          ! filename
   type(var_meta),intent(in)       :: vdata_meta(:)  ! veg data meta
   character(*), intent(in)        :: dname_vpoly    ! dimension name for polygon
   ! input-output
   type(namevar), intent(inout)    :: vdata(:)       ! veg data container
   ! output variables
   integer(i4b), intent(out)       :: nVpoly         ! number of veg polygons
   integer(i4b), intent(out)       :: ierr           ! error code
   character(*), intent(out)       :: message        ! error message
   ! local variables
   integer(i4b)                    :: iVar           ! variable index
   integer(i4b)                    :: ncid           ! NetCDF file ID
   integer(i4b)                    :: idimID_poly    ! dimension ID for HRUs
   integer(i4b)                    :: iVarID         ! variable ID
   integer(i4b)                    :: nMonth=12      ! number of month 
 
   ! initialize error control
   ierr=0; message='getVegData/'
   ! open file for reading
   ierr = nf90_open(fname, nf90_nowrite, ncid)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   ! get the ID of the poly dimension 
   ierr = nf90_inq_dimid(ncid, dname_vpoly, idimID_poly)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_vpoly); return; endif
   ! get the length of the poly dimension
   ierr = nf90_inquire_dimension(ncid, idimID_poly, len=nVpoly)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_vpoly); return; endif
   ! ** read in veg poly variables
   do iVar=1,size(vdata)
     ! get the variable ID
     ierr = nf90_inq_varid(ncid, trim(vdata_meta(ivar)%varName), ivarID)
     if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vdata_meta(ivar)%varName); return; endif
     select case(vdata_meta(iVar)%vartype)
     case('integer')
       select case(vdata_meta(iVar)%vardims)
       case('2D')
         ! allocate space for the 2D integer array
         allocate(vdata(ivar)%ivar2(nMonth,nVpoly),stat=ierr)
         if(ierr/=0)then; message=trim(message)//'err allocating 2D int space for vdata data structure'; return; endif
         ! get the data 
         ierr = nf90_get_var(ncid, ivarID, vdata(ivar)%ivar2)
       case('1D')
         ! allocate space for the 1D integer array
         allocate(vdata(ivar)%ivar1(nVpoly),stat=ierr)
         if(ierr/=0)then; message=trim(message)//'err allocating 2D dbl space for vdata data structure'; return; endif
         ! get the data 
         ierr = nf90_get_var(ncid, ivarID, vdata(ivar)%ivar1)
       end select
     case('double')
       select case(vdata_meta(iVar)%vardims)
       case('2D')
         ! allocate space for the 2D integer array
         allocate(vdata(ivar)%dvar2(nMonth,nVpoly),stat=ierr)
         if(ierr/=0)then; message=trim(message)//'err allocating 2D int space for vdata data structure'; return; endif
         ! get the data 
         ierr = nf90_get_var(ncid, ivarID, vdata(ivar)%dvar2)
       case('1D')
         ! allocate space for the 1D integer array
         allocate(vdata(ivar)%dvar1(nVpoly),stat=ierr)
         if(ierr/=0)then; message=trim(message)//'err allocating 2D int space for vdata data structure'; return; endif
         ! get the data 
         ierr = nf90_get_var(ncid, ivarID, vdata(ivar)%dvar1)
       end select
     end select
     if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   end do ! (looping through variables) 
   ! close the NetCDF file
   ierr = nf90_close(ncid)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end subroutine getVegData
 
 ! *****
 ! Public Subroutine: Read in MODIS land cover IGBP class text data...
 ! *********************************************
  subroutine getVegClassLookup(in_vegTable,   &  ! input: look up table to map veg class and veg properties
                               nVclass,       &  ! input:
                               vegClass,      &  ! output: veg class
                               vegPropt,      &  ! output: veg properties in 
                               err, message)    ! output: number of layer
 
 ! Purpose: read look-up table to map veg class and veg properties
   
   implicit none
 
   ! input
   character(*),intent(in)           :: in_vegTable       ! ascii table that containing veg class and properties
   integer(i4b),intent(in)           :: nVclass            ! number of veg classes (e.g., MODIS-IGBP=18)
   ! output
   integer(i4b),intent(out)          :: vegClass(:)       ! storage of veg class ID
   type(var_d),intent(out)           :: vegPropt(:)       ! storage of property value for each veg class
   integer(i4b),intent(out)          :: err                ! error code
   character(*),intent(out)          :: message            ! error message
   ! local
   character(len=strLen)             :: spName,spData      ! name and data from cLines(iLine) veg data, respectively
   character(len=strLen),allocatable :: sLines(:)          ! vector of character strings
   integer(i4b)                      :: ibeg_name          ! start index of variable name in string sLines(iLine)
   integer(i4b)                      :: iend_name          ! end index of variable name in string sLines(iLine)
   integer(i4b)                      :: iVclass            ! ID (= index) of veg class 
   integer(i4b)                      :: iLine              ! loop index of line in sLines
   integer(i4b)                      :: iunit              ! file unit
 
   err=0; message="getVegClassLookup/"
   !allocation
   do iVclass=1,nVclass; allocate(vegPropt(iVclass)%var(nPrpVeg)); end do
   ! open file (also returns un-used file unit used to open the file)
   call file_open(in_vegTable,iunit,err,message)
   ! get a list of character strings from non-comment lines
   call get_vlines(iunit,sLines,err,message)
   ! close the file unit
   close(iunit)
   do iLine=1,size(sLines) ! looping through lines in the veg data file 
     ! identify start and end of the name and the data
     ibeg_name = index(sLines(iLine),'<'); if(ibeg_name==0) err=20
     iend_name = index(sLines(iLine),'>'); if(iend_name==0) err=20
     if(err/=0)then; message='problem disentangling sLines(iLine) [string='//trim(sLines(iLine))//']'; return; endif
     ! extract name of the information, and the information itself
     spName = adjustl(sLines(iLine)(ibeg_name:iend_name))
     spData = adjustl(sLines(iLine)(iend_name+1:))
     select case(trim(spName))
       ! define directories 
       case('<vclass>')
         read(spData,*,iostat=err) vegClass(:)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<nroot>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%nroot),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<snup>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%snup),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<rs>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%rs),iVclass=1,nVclass) 
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<mrs>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%mrs),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<leafDIM>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%leafDIM),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<can_top_h>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%can_top_h),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<can_bot_h>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%can_bot_h),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<c_veg>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%c_veg),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case('<maxMassVeg>')
         read(spData,*,iostat=err) (vegPropt(iVclass)%var(ixPrpVeg%maxMassVeg),iVclass=1,nVclass)
         if(err/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
       case default; message=trim(message)//'spName not recognized'; err=35; return
     end select
   end do  ! end of loop thru veg data 
   return 
 end subroutine
 
end module read_vegdata
