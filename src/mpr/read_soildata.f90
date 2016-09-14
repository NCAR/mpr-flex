module read_soildata

USE nrtype
USE netcdf
USE data_type 
USE multconst                                     ! Including common constant (physical constant, other e.g., missingVal, etc.)
USE ascii_util,only:file_open                     !
USE ascii_util,only:split_line                    !
USE ascii_util,only:get_vlines                    !
USE var_lookup,only:ixVarSoilData,nVarSoilData    ! index of soil polygon variables and number of variables
USE var_lookup,only:ixVarTopo,nVarTopo            ! index of topographic variables and number of variables
USE var_lookup,only:ixVarSoil,nVarSoil            ! index of topographic variables and number of variables
USE var_lookup,only:ixPrpSoil,nPrpSoil            ! index of soil variables and number of variables
USE def_slyrmod                                   ! Public variables- multipliers of soil layer thickness

implicit none

private

public::getData
public::read_scls
public::mod_hslyrs
public::map_scls2prp
public::get_topoinfo

contains

! *********************************************************************
! Subroutine: get soil data for polygons and layers 
! *********************************************************************
 subroutine getData(fname,       &   ! input: file name
                    sdata_meta,  &   ! input: soil data meta
                    dname_spoly, &   ! input: dimension name for soil polygon
                    dname_slyrs, &   ! input: dimension name for soil layer 
                    sdata,       &   ! input-output: soil data structure
                    nSpoly,      &   ! output: number of soil polygon
                    nSLyrs,      &   ! output: number of soil layer  
                    ierr, message)   ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname          ! filename
 type(var_info),intent(in)       :: sdata_meta(:)  ! soil data meta
 character(*), intent(in)        :: dname_spoly    ! dimension name for polygon
 character(*), intent(in)        :: dname_slyrs    ! dimension name for layer 
 ! input-output
 type(namevar), intent(inout)    :: sdata(:)       ! soil data container
 ! output variables
 integer(i4b), intent(out)       :: nSpoly         ! number of soil polygons
 integer(i4b), intent(out)       :: nSlyrs         ! number of soil layers 
 integer(i4b), intent(out)       :: ierr           ! error code
 character(*), intent(out)       :: message        ! error message
 ! local variables
 integer(i4b)                    :: iVar           ! variable index
 integer(i4b)                    :: ncid           ! NetCDF file ID
 integer(i4b)                    :: idimID_poly    ! dimension ID for HRUs
 integer(i4b)                    :: idimID_slyr    ! dimension ID for stream segments
 integer(i4b)                    :: iVarID         ! variable ID
 ! initialize error control
 ierr=0; message='getData/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ! get the ID of the poly dimension
 ierr = nf90_inq_dimid(ncid, dname_spoly, idimID_poly)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_spoly); return; endif
 ! get the length of the poly dimension
 ierr = nf90_inquire_dimension(ncid, idimID_poly, len=nSpoly)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_spoly); return; endif
 ! get the ID of the soil layer dimension
 ierr = nf90_inq_dimid(ncid, dname_slyrs, idimID_slyr)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_slyrs); return; endif
 ! get the length of the soil layer dimension
 ierr = nf90_inquire_dimension(ncid, idimID_slyr, len=nSlyrs)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(dname_slyrs); return; endif
 ! ** read in soil poly variables
 do iVar=1,size(sdata)
   ! get the variable ID
   ierr = nf90_inq_varid(ncid, trim(sdata_meta(ivar)%varName), ivarID)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(sdata_meta(ivar)%varName); return; endif

  select case(sdata_meta(iVar)%vartype)
    case('integer')
      select case(sdata_meta(iVar)%vardims)
        case('2D')
          ! allocate space for the 2D integer array 
          allocate(sdata(ivar)%ivar2(nSlyrs,nSpoly),stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating 2D int space for sdata data structure'; return; endif
          ! get the data
          ierr = nf90_get_var(ncid, ivarID, sdata(ivar)%ivar2)
        case('vector')
          ! allocate space for the 1D integer array 
          allocate(sdata(ivar)%ivar1(nSpoly),stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating 1D int space for sdata data structure'; return; endif
          ! get the data
          ierr = nf90_get_var(ncid, ivarID, sdata(ivar)%ivar1)
        end select
    case('double')
      select case(sdata_meta(iVar)%vardims)
        case('2D')
          ! allocate space for the 2D integer array 
          allocate(sdata(ivar)%dvar2(nSlyrs,nSpoly),stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating 2D real space for sdata data structure'; return; endif
          ! get the data
          ierr = nf90_get_var(ncid, ivarID, sdata(ivar)%dvar2)
        case('vector')
          ! allocate space for the 1D integer array 
          allocate(sdata(ivar)%dvar1(nSpoly),stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating 1D real space for sdata data structure'; return; endif
          ! get the data
          ierr = nf90_get_var(ncid, ivarID, sdata(ivar)%dvar1)
      end select
  end select 
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end do  ! (looping through variables)

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine getData

! *****
! Public Subroutine: Read in USDA soil class text data...
! *********************************************
 subroutine read_scls(in_soilTable,  &  ! input: look up table to map soil class and soil properties
                      nSclass,       &  ! input:
                      soilClass,     &  ! output: soil class
                      scls2prp,      &  ! output: soil properties in 
                      ierr, message)    ! output: error control 

! Purpose: read look-up table to map soil class and soil properties
  
  implicit none

  ! input
  character(*),intent(in)           :: in_soilTable       ! ascii table that containing soil class and properties
  integer(i4b),intent(in)           :: nSclass            ! number of soil classes (e.g., STATSGO=16)
  ! output
  integer(i4b),intent(out)          :: soilClass(:)       ! storage of soil class ID
  type(var_d),intent(out)           :: scls2prp(:)       ! storage of property value for each soil class
  integer(i4b),intent(out)          :: ierr               ! error code
  character(*),intent(out)          :: message            ! error message
  ! local
  character(len=strLen)             :: spName,spData      ! name and data from cLines(iLine) soil data, respectively
  character(len=strLen),allocatable :: sLines(:)          ! vector of character strings
  integer(i4b)                      :: ibeg_name          ! start index of variable name in string sLines(iLine)
  integer(i4b)                      :: iend_name          ! end index of variable name in string sLines(iLine)
  integer(i4b)                      :: iSclass            ! ID (= index) of soil class 
  integer(i4b)                      :: iLine              ! loop index of line in sLines
  integer(i4b)                      :: iunit              ! file unit

  ierr=0; message="read_scls/"
  !allocation
  do iSclass=1,nSclass; allocate(scls2prp(iSclass)%var(nPrpSoil-2)); end do
    
  ! open file (also returns un-used file unit used to open the file)
  call file_open(in_soilTable,iunit,ierr,message)
  ! get a list of character strings from non-comment lines
  call get_vlines(iunit,sLines,ierr,message)
  ! close the file unit
  close(iunit)
  do iLine=1,size(sLines) ! looping through lines in the soil data file 
    ! identify start and end of the name and the data
    ibeg_name = index(sLines(iLine),'<'); if(ibeg_name==0) ierr=20
    iend_name = index(sLines(iLine),'>'); if(iend_name==0) ierr=20
    if(ierr/=0)then; message='problem disentangling sLines(iLine) [string='//trim(sLines(iLine))//']'; return; endif
    ! extract name of the information, and the information itself
    spName = adjustl(sLines(iLine)(ibeg_name:iend_name))
    spData = adjustl(sLines(iLine)(iend_name+1:))
    select case(trim(spName))
      ! define directories 
      case('<sclass>')
        read(spData,*,iostat=ierr) soilClass(:)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of USDA_sclass'; return; endif
      case('<sand_frac>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%sand_frac),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<clay_frac>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%clay_frac),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<bulk_density>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%bulk_density),iSclass=1,nSclass) 
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<soil_density>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%soil_density),iSclass=1,nSclass) 
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<field_capacity>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%field_capacity),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<wilting_point>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%wilting_point),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<porosity>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%porosity),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<ks>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%ks),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<slope_ret_curve>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%slope_ret_curve),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<psi_sat>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%psi_sat),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case('<myu>')
        read(spData,*,iostat=ierr) (scls2prp(iSclass)%var(ixPrpSoil%myu),iSclass=1,nSclass)
        if(ierr/=0)then; message=trim(message)//'problem with internal read of '//trim(spName); return; endif
      case default; message=trim(message)//'spName not recognized'; ierr=35; return
    end select
  end do  ! end of loop thru soil data 

 end subroutine read_scls

! *****
! Subroutine: soil thickness mod
! *********************************************
 subroutine mod_hslyrs(sdata, &         ! input/output: data structure of soil data including soil layer thickness [m] 
                       ierr, message)   ! output: error control
  use globalData, only:parMaster, parSubset
  implicit none 

  ! Input/output
  type(namevar), intent(inout)    :: sdata(:)       ! soil data container
  ! local
  real(dp),allocatable            :: a_h_array(:)   ! temp vector of soil layer multiplier 
  real(dp),allocatable            :: soil_h_mod(:)  ! temp holder of modified soil layer thickness 
  integer(i4b)                    :: nSpoly         ! number of soil polygon 
  integer(i4b)                    :: nSlyrs         ! number of soil layer 
  integer(i4b)                    :: iSpoly         ! Loop index of soil polygon 
  integer(i4b)                    :: iSlyrs         ! Loop index of soil layers 
  integer(i4b), intent(out)       :: ierr           ! error code
  character(*), intent(out)       :: message        ! error message

  ! initialize error control
  ierr=0; message='mod_hslyrs/'
  
  nSpoly=size(sdata(ixVarSoilData%hslyrs)%dvar2,2)
  nSlyrs=size(sdata(ixVarSoilData%hslyrs)%dvar2,1)

  allocate(a_h_array(nSlyrs),stat=ierr); 
  if(ierr/=0)then; message=trim(message)//'problem with allocating a_h_array'; return; endif
  allocate(soil_h_mod(nSlyrs),stat=ierr); 
  if(ierr/=0)then; message=trim(message)//'problem with allocating soil_h_mod'; return; endif
  !if
  !  a_h_array = something 
  !else
  !  a_h_array = parMaster(ixPar%z1gamma1)                                 ! Layer Multiplier 
  !endif
  do iSpoly =1,nSpoly
    soil_h_mod = a_h_array*sdata(ixVarSoilData%hslyrs)%dvar2(:,iSpoly)  ! modified soil layer thickness [m] 
    sdata(ixVarSoilData%hslyrs)%dvar2(:,iSpoly) = soil_h_mod            ! reassign modified layer thickness in data structure
  end do

 end subroutine mod_hslyrs 
 
! *****
! Subroutine: Map soil properties from soil class with lookup table 
! *********************************************
 subroutine map_scls2prp(scls,       &  ! input:  soil class array for soil polygon and layer
                         scls2prp,   &  ! input:  soil class-properties mapping
                         sdata,      &
                         sprpSxySz)     ! output: soil properties array for soil polygon and layer

  implicit none 
  ! Input
  integer(i4b),intent(in)      :: scls(:,:)      ! soil class for polygon and layer 
  type(var_d), intent(in)      :: scls2prp(:)    ! soil properties for polygon and layer 
  type(namevar),intent(in)     :: sdata(:)
  ! Output
  type(namedvar2),intent(inout) :: sprpSxySz(:)   ! soil properties for polygon and layer 
  ! local
  integer(i4b)                :: iPrpSoil       ! Loop index of soil properties 
  integer(i4b)                :: iscls          ! Loop index of soil class 
  integer(i4b)                :: iSpoly         ! Loop index of soil polygon 
  integer(i4b)                :: nSpoly         ! number of soil polygons
  integer(i4b)                :: iSlyrs         ! Loop index of soil layer 
  integer(i4b)                :: nSlyrs         ! number of soil layers
  
  nSpoly=size(scls,2)
  nSlyrs=size(scls,1)
  do iPrpSoil=1,nPrpSoil ! go through soil properties
    do iSpoly = 1,nSpoly
      do iSlyrs = 1,nSlyrs
        if (iPrpSoil==ixPrpSoil%h) then
          sprpSxySz(iPrpSoil)%varData(iSlyrs,iSpoly) = sdata(ixVarSoilData%hslyrs)%dvar2(iSlyrs,iSpoly)
        else if (iPrpSoil==ixPrpSoil%z) then
          sprpSxySz(iPrpSoil)%varData(iSlyrs,iSpoly) = sum(sdata(ixVarSoilData%hslyrs)%dvar2(1:iSlyrs,iSpoly))
        else 
          iscls = scls(isLyrs,iSpoly)
          sprpSxySz(iPrpSoil)%varData(iSlyrs,iSpoly) = scls2prp(iscls)%var(iPrpSoil) ! assign soil properties based on soil class
        endif
      enddo
    enddo
  enddo

 end subroutine map_scls2prp

! *****
! Subroutine: Extract soil variables from soil data into soil data structure 
! *********************************************
 subroutine get_soilinfo(sdata,         &  ! input:  data structure for soil data containing all the data in soil data
                         soil,          &  ! output: topography properties array for soil polygon 
                         ierr, message)    ! output: error control
  implicit none 
  !input
  type(namevar),  intent(in)    :: sdata(:)        ! soil data containing all the variables 
  !output
  type(namedvar2),intent(inout) :: soil(:)         ! storage of topo property values in subset domain-1D(poly)
  integer(i4b), intent(out)     :: ierr            ! error code
  character(*), intent(out)     :: message         ! error message
  !local
  integer(i4b)                  :: iVarSoil        ! Loop index of topo properties 
  
  ! initialize error control
  ierr=0; message='get_soilinfo/'

  soil(1)%varName = 'hslyrs' 
  soil(2)%varName = 'sand_frc' 
  soil(3)%varName = 'silt_frc' 
  soil(4)%varName = 'clay_frc' 
  soil(5)%varName = 'bulk_density' 

  do iVarSoil = 1,nVarSoil
    select case (soil(iVarSoil)%varName)
      case('hslyrs')
        soil(ivarSoil)%varData = sdata(ixVarSoilData%hslyrs)%dvar2
      case('sand_frc')
        soil(iVarSoil)%varData = sdata(ixVarSoilData%sand_frc)%dvar2
      case('silt_frc')
        soil(iVarSoil)%varData = sdata(ixVarSoilData%silt_frc)%dvar2
      case('clay_frc')
        soil(iVarSoil)%varData = sdata(ixVarSoilData%clay_frc)%dvar2
      case('bulk_density')
        soil(iVarSoil)%varData = sdata(ixVarSoilData%bulk_density)%dvar2
      case default; message=trim(message)//'Soil variable name not recognized'; ierr=35; return
    end select
  end do

 end subroutine get_soilinfo

! *****
! Subroutine: Extract topo properties from soil data into topo data structure 
! *********************************************
 subroutine get_topoinfo(sdata,         &  ! input:  data structure for soil data containing all the data in soil data
                         topo,          &  ! output: topography properties array for soil polygon 
                         ierr, message)    ! output: error control
  implicit none 
  !input
  type(namevar),  intent(in)    :: sdata(:)         ! soil data containing all the variables 
  !output
  type(namedvar),intent(inout)  :: topo(:)          ! storage of topo property values in subset domain-1D(poly)
  integer(i4b), intent(out)     :: ierr             ! error code
  character(*), intent(out)     :: message          ! error message
  !local
  integer(i4b)                  :: iVarTopo         ! Loop index of topo properties 
  real(dp), parameter           :: slope_min=0.01_dp ! minimum slope angle [percent]
  real(dp),allocatable          :: slope(:)
  
  ! initialize error control
  ierr=0; message='get_topoinfo/'

  topo(1)%varName = 'ele_mean' 
  topo(2)%varName = 'ele_std' 
  topo(3)%varName = 'slp_mean' 
  
  ! put threshold for lower limit of slope percentage
  allocate(slope(size(sdata(ixVarSoilData%slp_mean)%dvar1)),stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for slope';return;endif  
  slope=sdata(ixVarSoilData%slp_mean)%dvar1
  where ( slope > -0.01_dp .and. slope < slope_min )  !if computed value is less than min, set slope_min
    slope=slope_min
  end where

  do iVarTopo = 1,nVarTopo
    select case (topo(iVarTopo)%varName)
      case('ele_mean')
        topo(iVarTopo)%varData = sdata(ixVarSoilData%ele_mean)%dvar1
      case('ele_std')
        topo(iVarTopo)%varData = sdata(ixVarSoilData%ele_std)%dvar1
      case('slp_mean')
        topo(iVarTopo)%varData = slope 
      case default; message=trim(message)//'topo variable name not recognized'; ierr=35; return
    end select
  end do

 end subroutine get_topoinfo

end module read_soildata
