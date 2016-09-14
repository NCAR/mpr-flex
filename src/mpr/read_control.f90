module read_control 

! Data type definition
 USE nrtype                                       ! Including numerical type definition 
 USE data_type                                    ! Including custum data structure definition
 USE public_var                                   ! Including common constant (physical constant, other e.g., missingVal, etc.)
 USE ascii_util,only:file_open                    ! routine to open text file
 USE ascii_util,only:split_line                   !
 USE ascii_util,only:get_vlines                   !

subroutine get_control(ctlfile, err, message)
  implicit none
  ! input
  character(len=strLen),intent()    :: ctlfile 
  ! output 
  integer(i4b),         intent(out) :: err                     ! error code
  character(len=strLen),intent(out) :: message                 ! error message of downwind routine
  
  ! define directories
  character(len=strLen),save        :: input_dir                ! directory containing input data
  character(len=strLen),save        :: output_dir               ! directory containing output data
  character(len=strLen),save        :: nml_dir                  ! directory containing namelist file 
  ! Soil data stuff 
  character(len=strLen),save        :: fname_soil               ! filename containing soil data (e.g. USDA soil class id for 11 layers) 
  character(len=strLen),save        :: sclass_table             ! filename containing soil class id and soil properties 
  integer(i4b)         ,save        :: nSclass                  ! number of soil classes
  character(len=strLen),save        :: dname_spoly              ! dimension name of the soil polygon 
  character(len=strLen),save        :: dname_slyrs              ! dimension name of the soil layer
  character(len=32)    ,save        :: sprp_comp                ! method to compute soil properties - either "pedotf" or "lookup" 
  ! vege data stuff
  character(len=strLen),save        :: fname_veg                ! filename containing vege data (e.g. land cover id, monthly lai, and monthly green frac) 
  character(len=strLen),save        :: dname_vpoly              ! dimension name of the vege polygon 
  character(len=strLen),save        :: vclass_table             ! filename containing veg class id and veg properties 
  integer(i4b)         ,save        :: nVclass                  ! number of vegetation types
  ! mapping stuff  
  character(len=strLen),save        :: fname_smapping           ! filename containing weight mapping (soil poly -> model hru) 
  character(len=strLen),save        :: fname_vmapping           ! filename containing weight mapping (vege poly -> model hru) 
  character(len=strLen),save        :: dname_overSpoly          ! dimension name of the maping weight - overlapping soil polygon 
  character(len=strLen),save        :: dname_overVpoly          ! dimension name of the maping weight - overlapping vege polygon
  character(len=strLen),save        :: dname_hru                ! dimension name of soil maping weight - hru
  ! define output file
  character(len=strLen),save        :: soilnm_output            ! name of ascii output file for model soil parameters
  character(len=strLen),save        :: soilnc_output            ! name of netCDF output file for model soil parameters
  character(len=strLen),save        :: vegenc_output            ! name of netCDF output file for model soil parameters
  ! Model stuff 
  character(len=strLen),save        :: model_name               ! name of hydrologic model 
  character(len=strLen),save        :: dname_mhru               ! dimension name for model HRU
  character(len=strLen),save        :: dname_mlyr               ! dimension name for model layer
  integer(i4b),         save        :: nMLyr                    ! Number of model soil layers
  ! namelist
  character(len=strLen),save        :: param_nml                ! filename of the namelist containing transfer function parameters 
  ! No saved Miscleneous local
  character(len=strLen)             :: cfile_name               ! name of the control file
  character(len=strLen),allocatable :: cLines(:)                ! vector of character strings
  integer(i4b)                      :: iunit                    ! file unit
  integer(i4b)                      :: iLine                    ! index of line in cLines
  integer(i4b)                      :: ibeg_name                ! start index of variable name in string cLines(iLine)
  integer(i4b)                      :: iend_name                ! end index of variable name in string cLines(iLine)
  integer(i4b)                      :: iend_data                ! end index of data in string cLines(iLine)
  character(len=strLen)             :: cName,cData              ! name and data from cLines(iLine)

  ! get command-line argument defining the full path to the control file
  call getarg(1,cfile_name)
  if(len_trim(cfile_name)==0) call handle_err(50,'need to supply name of the control file as a command-line argument')
  ! *** get a list of character strings from non-comment lines ****
  ! open file (also returns un-used file unit used to open the file)
  call file_open(trim(cfile_name),iunit,ierr,cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
  ! get a list of character strings from non-comment lines
  call get_vlines(iunit,cLines,ierr,cmessage)
  if(ierr/=0) call handle_err(ierr, cmessage)
  ! close the file unit
  close(iunit)
  ! loop through the non-comment lines in the input file, and extract the name and the information
  do iLine=1,size(cLines)
    ! identify start and end of the name and the data
    ibeg_name = index(cLines(iLine),'<'); if(ibeg_name==0) ierr=20
    iend_name = index(cLines(iLine),'>'); if(iend_name==0) ierr=20
    iend_data = index(cLines(iLine),'!'); if(iend_data==0) ierr=20
    if(ierr/=0) call handle_err(60,'problem disentangling cLines(iLine) [string='//trim(cLines(iLine))//']')
    ! extract name of the information, and the information itself
    cName = adjustl(cLines(iLine)(ibeg_name:iend_name))
    cData = adjustl(cLines(iLine)(iend_name+1:iend_data-1))
    print*, trim(cName), ' --> ', trim(cData)
    ! populate variables
    select case(trim(cName))
      ! PART 1: DEFINE DIRECTRIES
      case('<input_dir>');     input_dir       = trim(cData)    ! directory containing input data
      case('<output_dir>');    output_dir      = trim(cData)    ! directory containing output data
      case('<nml_dir>');       nml_dir         = trim(cData)    ! directory containing parameter namelist 
      ! PART 2: SPECIFICATION FOR SOIL DATA
      case('<fname_soil>');    fname_soil      = trim(cData)    ! name of netcdf containing soil data (soil class id) 
      case('<sclass_table>');  sclass_table    = trim(cData)    ! name of text file containing soil class- soil property look up table 
      case('<nsclass>'); read(cData,*,iostat=ierr) nSclass      ! number of soil classes
        if(ierr/=0) call handle_err(70,'problem with internal read of nSclass, read from control file')
      case('<dname_spoly>');   dname_spoly     = trim(cData)    ! dimension name of the soil polygons 
      case('<dname_slyrs>');   dname_slyrs     = trim(cData)    ! dimension name of the soil layer
      case('<sprp_comp>');     sprp_comp       = trim(cData)    ! method to compute soil properties- either "pedotf" or "lookup"
      ! PART 3: SPECIFICATION FOR VEGE DATA
      case('<fname_veg>');     fname_veg       = trim(cData)    ! name of netCDF containing vege data (MODIS) 
      case('<dname_vpoly>');   dname_vpoly     = trim(cData)    ! dimension name of the soil polygons 
      case('<vclass_table>');  vclass_table    = trim(cData)    ! name of text file containing veg class- veg property look up table 
      case('<nvclass>'); read(cData,*,iostat=ierr) nVclass      ! number of soil classes
        if(ierr/=0) call handle_err(70,'problem with internal read of nVclass, read from control file')
      ! PART 4: SPECIFICATION FOR MAPPING FROM GEOPHYSICAL DATA POLY TO MODEL HRU 
      case('<fname_smapping>'); fname_smapping = trim(cData)    ! name of weight netcdf to map from soil-poly to hru 
      case('<fname_vmapping>'); fname_vmapping = trim(cData)    ! name of weight netcdf to map from vege-poly to hru
      case('<dname_overspoly>');dname_overSpoly= trim(cData)    ! 
      case('<dname_overvpoly>');dname_overVpoly= trim(cData)    ! 
      case('<dname_hru>');      dname_hru      = trim(cData)    ! hru dimension name in mapping data 
      ! PART 5: DEFINE SPECIFICATION OF HYDROLOGIC MODEL  
      case('<model_name>');     model_name     = trim(cData)    ! name of hydrologic model
      case('<model_layer>'); read(cData,*,iostat=ierr) nMLyr  ! 
        if(ierr/=0) call handle_err(70,'problem with internal read of nMLyr, read from control file')
      case('<dname_mhru>');    dname_mhru      = trim(cData)    ! dimension name of the HRUs
      case('<dname_mlyr>');    dname_mlyr      = trim(cData)    ! dimension name of the HRUs
      ! PART 6: DEFINE OUTPUT FILE 
      case('<soilnc_output>'); soilnc_output   = trim(cData)    ! netcdf filename for the model output
      case('<vegenc_output>'); vegenc_output   = trim(cData)    ! netcdf filename for the model output
      case('<soilnm_output>'); soilnm_output   = trim(cData)    ! text filename for the model output
      ! PART 7: Namelist file name 
      case('<param_nml>');     param_nml       = trim(cData)    ! name of namelist including TF parameters 
    end select
  end do  ! looping through lines in the control file
  
  return
end subroutine read_control 
