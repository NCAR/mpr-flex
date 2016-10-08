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

!Following accessible outside this module
public::getVegData
public::getvegClassLookup
public::map_vcls2prp

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
       case('vector')
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
       case('vector')
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
                               nVclass,      &  ! input:
                               vegClass,     &  ! output: veg class
                               vegPropt,     &  ! output: veg properties in 
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
 
 end subroutine getVegClassLookup
 
 ! *****
 ! Public Subroutine: Map veg properties from veg class with lookup table  
 ! *********************************************
 subroutine map_vcls2prp(vdata,   &  ! input: veg data holder containing all the vege data variables 
                         vcls2prp,&  ! input: veg class-property mapping 
                         vegClass,&  ! input: list of veg class id
                         vprp,    &  ! output: veg property array for veg polygon
                         ierr,    &  ! error code
                         message)    ! error message   
   implicit none
   !input
   type(namevar),intent(in)      :: vdata(:)
   type(var_d),  intent(in)      :: vcls2prp(:) ! veg properties for polygon 
   integer(i4b), intent(in)      :: vegClass(:) ! veg properties for polygon 
   !output
   type(namedvar),intent(inout)  :: vprp(:)     ! veg properties for polygon
   integer(i4b), intent(out)     :: ierr      ! error code
   character(*), intent(out)     :: message   ! error message

   !local
   integer(i4b)                  :: iPrpVeg       ! Loop index of veg properties 
   integer(i4b)                  :: ivcls         ! Loop index of veg class 
   integer(i4b)                  :: iSelect       ! matching indix  
   integer(i4b)                  :: iVpoly        ! Loop index of veg polygon 
   integer(i4b)                  :: nVpoly        ! number of veg polygons
 
   nVpoly=size(vdata(ixVarVegData%vegclass)%ivar1)
   ! from straight from netCDF
   do iVpoly = 1,nVpoly
     vprp(ixPrpVeg%lai01)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(1, iVpoly)
     vprp(ixPrpVeg%lai02)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(2, iVpoly)
     vprp(ixPrpVeg%lai03)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(3, iVpoly)
     vprp(ixPrpVeg%lai04)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(4, iVpoly)
     vprp(ixPrpVeg%lai05)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(5, iVpoly)
     vprp(ixPrpVeg%lai06)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(6, iVpoly)
     vprp(ixPrpVeg%lai07)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(7, iVpoly)
     vprp(ixPrpVeg%lai08)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(8, iVpoly)
     vprp(ixPrpVeg%lai09)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(9, iVpoly)
     vprp(ixPrpVeg%lai10)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(10,iVpoly)
     vprp(ixPrpVeg%lai11)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(11,iVpoly)
     vprp(ixPrpVeg%lai12)%varData(iVpoly) = vdata(ixVarVegData%lai)%dvar2(12,iVpoly)
     vprp(ixPrpVeg%vegtype)%varData(iVpoly) = vdata(ixVarVegData%vegclass)%ivar1(iVpoly)
   enddo
   do iPrpVeg=14,nPrpVeg ! go through veg properties
     do iVpoly = 1,nVpoly
        ivcls = vdata(ixVarVegData%vegclass)%ivar1(iVpoly)
        call findix(ivcls,vegClass,iSelect,ierr,message)  !find index in vegClass array matching iVclass
        vprp(iPrpVeg)%varData(iVpoly) = vcls2prp(iSelect)%var(iPrpVeg) ! assign veg properties based on veg class
     enddo
   enddo

   return   

   contains

     subroutine findix(scl,vec,iSelect, ierr, message)
       ! Find vec index where the value match up with scl  
       implicit none
     
       integer(i4b),intent(in)              :: scl
       integer(i4b),intent(in)              :: vec(:)
       integer(i4b),intent(out)             :: iSelect
       integer(i4b)                         :: i(1)
       integer(i4b), intent(out)            :: ierr      ! error code
       character(*), intent(out)            :: message   ! error message

       ! initialize error control
       ierr=0; message='findix/' 

       i = minloc(abs(vec - scl))
       iSelect = i(1)  ! de-vectorize the desired stream segment
       if(vec(iSelect) /= scl)&
         ierr=60; message=trim(message)//'unable to find matched value'; return  
     end subroutine findix
   
 end subroutine map_vcls2prp

end module read_vegdata
