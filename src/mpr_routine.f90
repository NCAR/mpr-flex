module mpr_routine 

  use nrtype                                           ! Including numerical type definition 
  use data_type                                        ! Including custum data structure definition
  use public_var                                       ! Including common constant (physical constant, other e.g., missingVal, etc.)

  implicit none

  private

  public:: run_mpr
  public:: mpr

contains

! ************************************************************************************************
! Public subroutine: run MPR and save estimated parameters in netCDF
! ************************************************************************************************
! this subroutine is used for opt = 2 in namelist (run only mpr and output parameters)
subroutine run_mpr( calParam, restartFile, err, message ) 
  use globalData,    only: betaCalScale, parSubset, betaInGamma, gammaSubset, nBetaGamma
  use model_wrapper, only: read_hru_id
  use write_param_nc,only: defSoilNetCDF, write_vec_ivar, write_array2_dvar
  implicit none
  ! input variables
  real(dp),             intent(in)  :: calParam(:)                   ! parameter in namelist, not necessarily all parameters are calibrated
  character(len=strLen),intent(in)  :: restartFile                   ! name of restart file including iteration, the most recent parameter values 
  ! output variables
  integer(i4b),         intent(out) :: err                           ! error id 
  character(len=strLen),intent(out) :: message                       ! error message
  ! local
  type(var_d)                       :: calParStr(nBetaGamma)         ! parameter storage including perLayr values converted from parameter array 
  type(var_d)                       :: pnormCoef(size(betaCalScale)) ! parameter storage converted from parameter array 
  type(var_d),          allocatable :: paramGammaStr(:)              ! calibratin gamma parameter storage extracted from calParStr
  integer(i4b)                      :: idummy                        ! dummy vaiable
  integer(i4b)                      :: idx                           ! counter 
  integer(i4b)                      :: iLyr                          ! loop index for model layer 
  integer(i4b)                      :: iPar                          ! loop index for parameter 
  integer(i4b)                      :: nVegParModel                  ! Number of model vege parameters associated with calibrating gamma parameter 
  integer(i4b)                      :: nSoilParModel                 ! Number of model soil parameters associated with calibrating gamma parameter 
  logical(lgc),         allocatable :: mask(:)                       ! 1D mask
  integer(i4b)                      :: hruID(nHru)                   ! Hru ID
  real(dp),             allocatable :: hModel(:,:)                   ! storage of model layer thickness at model layer x model hru 
  real(dp),             allocatable :: params(:)                     ! parameter vector that is input into mpr 
  type(namedvar2),      allocatable :: parMxyMz(:)                   ! storage of model soil parameter at model layer x model hru 
  type(namedvar2),      allocatable :: vegParMxy(:)                  ! storage of model vege parameter at model hru
  logical                           :: isExistFile                   ! logical to check if the file exist or not
  character(len=strLen)             :: cmessage                      ! error message from subroutine

  err=0; message='run_mpr/' ! to initialize error control
  if ( idModel/=0 )then;idModel=0;print*,trim(message)//'idModel is set to zero - model inepenedent';endif
  allocate(params, source=calParam) ! copy calParameter default 
  inquire(file=trim(adjustl(restartFile)), exist=isExistFile)
  if ( isExistFile ) then !  if state file exists, update calParam, otherwise use default value
    print*, 'read restart file'
    open(unit=70,file=trim(adjustl(restartFile)), action='read', status = 'unknown')
    read(70,*) idummy ! restart file include iStart  
    read(70,*) (params(iPar),iPar=1,nBetaGamma)    
    close(70)
  endif
  idx=1
  do iPar=1,nBetaGamma
    if (parSubset(iPar)%perLyr)then
      allocate(calParStr(iPar)%var(nLyr))
      calParStr(iPar)%var=params(idx:idx+nLyr-1)
      idx=idx+nLyr
    else
      allocate(calParStr(iPar)%var(1))
      calParStr(iPar)%var=params(idx)
      idx=idx+1
    endif
  end do
  do iPar=1,size(betaCalScale) ! put calpar vector from optimization routine output pnorm coef. data strucure
    pnormCoef(iPar)%var=calParam(idx:idx+1)
    idx=idx+2
  end do
  call read_hru_id(idModel, hruID, err, cmessage)    ! to get hruID
  if (err/=0)then;message=trim(message)//trim(cmessage);return;endif
  if ( any(parSubset(:)%beta /= "beta") )then ! calPar includes gamma parameters to be used for MPR 
    nSoilParModel=size(betaInGamma)           ! number of soil parameters associated with gamma parameters
    nVegParModel=1                            ! number of vege parameters associated with gamma parameters (now 1 temporarily)
    allocate(hModel(nLyr,nHru),stat=err);      if(err/=0)then;message=trim(message)//'error allocating hModel';return;endif
    allocate(parMxyMz(nSoilParModel),stat=err);if(err/=0)then;message=trim(message)//'error allocating parMxyMz';return;endif
    allocate(vegParMxy(nVegParModel),stat=err);if(err/=0)then;message=trim(message)//'error allocating vegParMxy';return;endif
    do iPar=1,nSoilParModel
      allocate(parMxyMz(iPar)%varData(nLyr,nHru),stat=err)
    enddo
    do iPar=1,nVegParModel
      allocate(vegParMxy(iPar)%varData(nMonth,nHru),stat=err)
    enddo
    allocate(mask(nBetaGamma))
    mask=parSubset(:)%beta/="beta"
    allocate(paramGammaStr(count(mask)))
    paramGammaStr=pack(calParStr,mask)
    do iPar=1,size(paramGammaStr)
      if (size(paramGammaStr(iPar)%var)>1)then;message=trim(message)//'gammaParameter should not have perLayer value';return;endif
    enddo
    call mpr(hruID, pnormCoef, paramGammaStr, gammaSubset, hModel, parMxyMz, vegParMxy, err, cmessage) ! to output model layer thickness and model parameter via MPR
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    !! Write parameter derived from MPR in netCDF 
    call defSoilNetCDF(trim(mpr_output_dir)//trim(param_nc),nHru,nLyr,err,cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    call write_vec_ivar(trim(mpr_output_dir)//trim(param_nc),"hruid",hruID,1,err,cmessage) 
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    call write_vec_ivar(trim(mpr_output_dir)//trim(param_nc),"lyr",(/(iLyr,iLyr=1,nLyr)/),1,err,cmessage) 
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    do iPar=1,size(betaInGamma)
      call write_array2_dvar(trim(mpr_output_dir)//trim(param_nc),trim(betaInGamma(iPar)),parMxyMz(iPar)%varData,(/1,1/),(/nLyr,nHru/),err,cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    enddo
  else  
    print*,trim(message)//'there is no gamma pamameters in parameter input file to perform MPR';stop
  endif
  return
end subroutine

! ************************************************************************************************
! public subroutine: mpr 
! ************************************************************************************************
subroutine mpr(hruID,             &     ! input: hruID
               pnormCoefStr,      &     ! input: list of pnorm coefficients
               gammaParStr,       &     ! input: array of gamma parameter 
               gammaParMeta,      &     ! input: array of gamma parameter metadata
               hModel,            &     ! output: Model layer thickness
               parMxyMz,          &     ! output: MPR derived soil parameter
               vegParMxy,         &     ! output: MPR derived veg parameter
               err, message)            ! output: error id and message
  use model_wrapper,        only:read_hru_id
  use popMeta,              only:popMprMeta
  use globalData,           only:parMaster, betaInGamma, betaCalScale !THIS (betaCalScale) may be input 
  use globalData,           only:sdata_meta,vdata_meta,map_meta, vprp_meta
  use get_ixname,           only:get_ixPar
  use tf,                   only:comp_soil_model_param         ! Including Soil model parameter transfer function
  use modelLayer,           only:comp_model_depth              ! Including model layr depth computation routine 
  use modelLayer,           only:map_slyr2mlyr                 ! Including model layr computation routine 
  !use vegtf,                only:comp_veg_model_param          ! Including Veg model parameter transfer function
  use upscaling,            only:aggreg                        ! Including Upscaling operator 
  use read_mapdata,         only:getMapData                    ! routine to read mapping data into data structures
  use read_vegdata,         only:getVegData                    ! routine to read veg data into data structures
  use read_vegdata,         only:getVegClassLookup             ! routine to read veg calss-property lookupu table 
  use read_soildata,        only:getData                       ! routine to read soil data into data structures
  use read_soildata,        only:mod_hslyrs                    ! routine to modify soil layer thickness and updata soil data structure
  use var_lookup,           only:ixPar,nPar                    ! 
  USE var_lookup,           only:ixVarSoilData,nVarSoilData    ! index of soil data variables and number of variables 
  USE var_lookup,           only:ixVarVegData,nVarVegData      ! index of vege data variables and number of variables
  USE var_lookup,           only:ixVarMapData,nVarMapData      ! index of map data variables and number of variables
  USE var_lookup,           only:ixPrpVeg,nPrpVeg              ! index of veg properties and number of properties

   implicit none

  ! input
  integer(i4b),         intent(in)   :: hruID(:)                 ! hruID list
  type(var_d),          intent(in)   :: gammaParStr(:)           ! data structure of gamma parameter value adjusted with calibration
  type(var_d),          intent(in)   :: pnormCoefStr(:)          ! data structure of pnorm coefficient value adjusted with calibration
  type(cpar_meta),      intent(in)   :: gammaParMeta(:)          ! array of calibrating meta data
  ! output
  real(dp),             intent(out)  :: hModel(:,:)              ! Model layer thickness at model layer x model hru 
  type(namedvar2),      intent(inout):: parMxyMz(:)              ! storage of model soil parameter at model layer x model hru 
  type(namedvar2),      intent(inout):: vegParMxy(:)             ! storage of model vege parameter at model hru
  integer(i4b),         intent(out)  :: err                      ! error code
  character(len=strLen),intent(out)  :: message                  ! error message 
  ! local
  character(len=strLen)              :: cmessage                 ! error message from downward subroutine
  integer,     parameter             :: iHruPrint = 1            ! model hru id for which everything is printed for checking
  integer(i4b),parameter             :: nSub=11                  ! max. number of Soil layer within Model layer
  integer(i4b)                       :: iLocal                   ! index of hru array in mapping file that match hru id of interest 
  integer(i4b)                       :: iDummy(1)                ! 1D integer array for temporal storage 
  integer(i4b)                       :: iGamma                   ! index loop
  integer(i4b)                       :: iPoly                    ! Loop index of soil polygon
  integer(i4b)                       :: iSLyr                    ! Loop index of soil layer 
  integer(i4b)                       :: iMLyr                    ! loop index of model soil layer (1,...,n from top to bottom) 
  integer(i4b)                       :: iParm                    ! Loop index of model parameters (e.g., VIC)
  integer(i4b)                       :: iVar                     ! Loop index of miscleneous variables 
  integer(i4b)                       :: iPrpVeg                  ! Loop index of veg properties 
  integer(i4b)                       :: iHru                     ! loop index of hrus 
  integer(i4b)                       :: iSub                     ! Loop index of multiple soi layers in model layer
  logical(lgc),allocatable           :: mask(:)                  ! mask for 1D array 
  logical(lgc),allocatable           :: vmask(:)                 ! mask for vpolIdSub array 
  type(par_meta),allocatable         :: gammaParMasterMeta(:)
  integer(i4b)                       :: nVegParModel             ! Number of model vege parameters
  integer(i4b)                       :: nSoilParModel            ! Number of model soil parameters
  integer(i4b)                       :: nSpoly                   ! number of soil polygon in entire soil data domain 
  integer(i4b)                       :: nSlyrs                   ! number of soil layers
  integer(i4b)                       :: nVpoly                   ! number of vege polygon (grid box) in entire vege data domain 
  integer(i4b)                       :: nShru                    ! number of hrus in soil mapping file) 
  integer(i4b)                       :: nVhru                    ! number of hrus in vege mapping file (should be equal to nVhru)
  integer(i4b)                       :: nOverSpoly               ! number of overlapped soil polygon for each model HRU 
  integer(i4b)                       :: nOverVpoly               ! number of overlapped vege polygon for each model HRU 
  integer(i4b)                       :: nSpolyLocal              ! number of subset overlapped soil polygon for each model HRU 
  integer(i4b)                       :: nVpolyLocal              ! number of subset overlapped vege polygon for each model HRU 
  real(dp)                           :: hmult                    ! mulitplier of soil layer  
  real(dp)                           :: hfrac(nLyr-1)            ! fraction of soil depth for each model layer 
  type(namevar)                      :: sdata(nVarSoilData)      ! soil data container for all the soil polygons
  type(namevar)                      :: sdataLocal(nVarSoilData) ! soil data container for local soil polygon 
  type(namevar)                      :: vdata(nVarVegData)       ! veg data container for all the veg polygons
  integer(i4b),    allocatable       :: vegClass(:)              ! veg class array (e.g., IGBP)
  type(var_d),     allocatable       :: vcls2prp(:)              ! storage of property value for each veg class
  integer(i4b)                       :: iVclass                  ! ID (= index) of vege class 
  real(dp),        allocatable       :: hModelLocal(:,:)         ! Model layer thickness for soil polygon within one hru
  real(dp),        allocatable       :: zModelLocal(:,:)         ! Model layer depth for soil polygon within one hru
  type(namedvar2), allocatable       :: parSxySz(:)              ! storage of model soil parameter for 2D field -soil layer x soil poy 
  type(namedvar2), allocatable       :: parSxyMz(:)              ! storage of model soil parameter for 2D field -model layer x soil poy
  type(namedvar2), allocatable       :: parVxy(:)                ! storage of model vege parameter for 1D or 2D field - vege poly (x month)
  integer(i4b),    allocatable       :: polySub(:)               ! list of ID (=index) of soil polygons contributing model hru
  integer(i4b),    allocatable       :: vPolySub(:)              ! list of ID (=index) of veg polygons contributing model hru 
  type(mapvar)                       :: mapdata(2)               ! map data container for all the soil polygons
  real(dp),        allocatable       :: swgtSub(:)               ! adjusted Areal weight of soil polygon for all model hrus
  real(dp),        allocatable       :: vwgtSub(:)               ! adjusted Areal weight of veg polygon for one model hrus
  type(poly),      allocatable       :: soil2model_map(:)        ! data structure to hold weight and index of contributing soil layer per model layer and per soil polygon
  type(lyr_d),     allocatable       :: paramvec(:)              !
  
  ! initialize error control
  err=0; message='mpr/'
  !(0) Preparation
  allocate(gammaParMasterMeta, source=parMaster) ! copy master parameter metadata
  nSoilParModel=size(betaInGamma)                ! number of soil and vege parameters associated with gamma parameter
  nVegParModel=1
  ! Swap gammaParMasterMeta%val with gammaPar value
  do iGamma=1,size(gammaParStr)
    gammaParMasterMeta(gammaParMeta(iGamma)%ixMaster)%val=gammaParStr(iGamma)%var(1)
  enddo
  do iParm=1,size(pnormCoefStr)
    associate(ix=>get_ixPar(betaCalScale(iParm)%betaname))
    gammaParMasterMeta(ix)%hpnorm=pnormCoefStr(iParm)%var(1)
    gammaParMasterMeta(ix)%vpnorm=pnormCoefStr(iParm)%var(2)
    end associate
  enddo 
  call popMprMeta( err, cmessage)   !for sdata_meta, vdata_meta, map_meta
  if(err/=0)then; message=trim(message)//cmessage; return; endif
  hmult=gammaParMasterMeta(ixPar%z1gamma1)%val
  call pop_hfrac(gammaParStr, gammaParMeta, hfrac, err, cmessage) ! to get hfrac 
  if(err/=0)then; message=trim(message)//cmessage; return; endif
  ! Memory allocation
  allocate(parSxySz(nSoilParModel),stat=err);if(err/=0)then; message=trim(message)//'error allocating parSxySz'; return; endif
  allocate(parSxyMz(nSoilParModel),stat=err);if(err/=0)then; message=trim(message)//'error allocating parSxyMz'; return; endif
  allocate(paramvec(nSoilParModel),stat=err);if(err/=0)then; message=trim(message)//'error allocating paramvec'; return; endif 
  ! (1) Get Geophysical data 
  ! *****
  ! (1.1) Get soil data  
  ! *********************************************
  call getData(trim(mpr_input_dir)//trim(fname_soil),& ! input: soil data (netCDF)
               sdata_meta,                           & ! input: soil data meta
               dname_spoly,                          & ! input: spatial dimension (polygon ID)
               dname_slyrs,                          & ! input: spatial dimension (polygon ID)
               sdata,                                & ! input-output: soil data structure
               nSpoly,                               & ! output: number of dimension (i.e. number of soil polygon)
               nSlyrs,                               & ! output: number of dimension (i.e. number of soil layer)
               err, cmessage)
  if(err/=0)then; message=trim(message)//cmessage; return; endif 
  call mod_hslyrs(sdata,hmult,err,cmessage) ! modify soil layer thickness in sdata data structure
  if(err/=0)then; message=trim(message)//cmessage; return; endif 
  ! *****
  ! (1.2)  veg class - properties lookup table ...
  ! *********************************************
  ! (1.2.1) Read in veg data netCDF...
   allocate(vegClass(nVclass),stat=err);if(err/=0)then;message=trim(message)//'error allocating vegClass'; return; endif 
   allocate(vcls2prp(nVclass),stat=err);if(err/=0)then;message=trim(message)//'error allocating vcls2prp'; return; endif
   do iVclass=1,nVclass
     allocate(vcls2prp(iVclass)%var(nPrpVeg),stat=err);if(err/=0)then;message=trim(message)//'error allocating vcls2prp%var';return;endif
   enddo
   call getVegData(trim(mpr_input_dir)//trim(fname_veg), & ! input: file name
                   vdata_meta,                           & ! input: veg data meta
                   dname_vpoly,                          & ! input: dimension name for veg polygon
                   vdata,                                & ! input-output: veg data structure
                   nVpoly,                               & ! output: number of veg polygon
                   err, cmessage)                          ! output: error control
   if(err/=0)then; message=message//cmessage; return; endif
  ! (1.2.2) Read in veg class-property lookup table 
   call getvegClassLookup(trim(mpr_input_dir)//trim(vclass_table), &
                          nVclass,                                 &
                          vegClass,                                &
                          vcls2prp,                                &
                          err, cmessage)
   if(err/=0)then; message=message//cmessage; return; endif
  ! *****
  ! (2.) Read in mapping netcdf 
  ! *********************************************
  call getMapData(trim(mpr_input_dir)//trim(fname_smapping), &   ! input: file name
                  'soil',                                    &   ! input: geophysical data type 
                  map_meta,                                  &   ! input: map data meta
                  dname_hru,                                 &   ! input: dimension name for hru 
                  dname_overSpoly,                           &   ! input: dimension name for overlap polygon 
                  mapdata,                                   &   ! input-output: map data structure
                  nShru,                                     &   ! output: number of hru 
                  nOverSpoly,                                &   ! output: max number of overlap polygons
                  err,cmessage)                                  ! output: error control
  if (err/=0)then; message=trim(message)//cmessage; return; endif
  call getMapData(trim(mpr_input_dir)//trim(fname_vmapping), &   ! input: file name
                  'veg',                                     &   ! input: geophysical data type 
                  map_meta,                                  &   ! input: map data meta
                  dname_hru,                                 &   ! input: dimension name for hru 
                  dname_overVpoly,                           &   ! input: dimension name for overlap Polygon 
                  mapdata,                                   &   ! input-output: map data structure
                  nVhru,                                     &   ! output: number of hru 
                  nOverVpoly,                                &   ! output: max number of overlap polygon
                  err,cmessage)                                  ! output: error control
  if (err/=0)then; message=message//cmessage; return; endif
  if ( nShru /= nVhru )then;err=10;message=trim(message)//'Different number of hru in vege and soil mapping file';return;endif  
  if ( opt==2 .and. nHru /= nShru )then;err=10;message=trim(message)//'nHru= '//trim(int2str(nShru))//' NOT '//trim(int2str(nHru));return;endif
  associate( hruMap      => mapdata(1)%var(ixVarMapData%hru_id)%ivar1,        &
             hruMapVege  => mapdata(2)%var(ixVarMapData%hru_id)%ivar1,        &
             swgt        => mapdata(1)%var(ixVarMapData%weight)%dvar2,        &
             overVpolyID => mapdata(2)%var(ixVarMapData%overlapPolyId)%ivar2, & 
             overSpolyID => mapdata(1)%var(ixVarMapData%overlapPolyId)%ivar2 )
  if ( minval(abs(hruMap-hruMapVege)) /= 0 )then;err=11;message=trim(message)//'different hru id in vege and soil mapping file';return;endif
  !!! ---------------------------------------------
  !!! Start of model hru loop (from mapping file) !!!
  !!! ---------------------------------------------
  hru: do iHru=1,nHru
    ! Get index (iLocal) of hru id that matches with current hru from hru id array in mapping file
    if (minval(abs(hruMap-hruID(iHru))) /= 0 )then;err=10;message=trim(message)//'hru id does not exist in mapping file';return; endif
    idummy = minloc( abs(hruMap-hruID(iHru)) )
    iLocal=idummy(1)
    ! Select list of soil polygons contributing a current hru
    allocate(mask(nOverSpoly))
    mask = ( overSpolyID(:,iLocal) /= imiss )
    allocate(polySub(count(mask)))
    allocate(swgtSub(count(mask)))
    nSpolyLocal = size(polySub)                    ! number of contributing soil polygons to current hru
    polySub     = pack(overSpolyID(:,iLocal),mask) ! id of soil polygons contributing to current hru
    swgtSub     = pack(swgt(:,iLocal),mask)        ! weight of soil polygons contributing to current hru
    ! allocate memmory
    do iParm=1,nSoilParModel
      allocate(parSxySz(iParm)%varData(nSlyrs,nSpolyLocal),stat=err); if(err/=0)then;message=message//'error allocating parSxySz%varData';return;endif 
      allocate(parSxyMz(iParm)%varData(nLyr,nSpolyLocal),stat=err);   if(err/=0)then;message=message//'error allocating parSxyMz%varData';return;endif
    enddo
    allocate(hModelLocal(nLyr,nSpolyLocal),stat=err);    if(err/=0)then;message=message//'error allocating hModelLocal'; return;endif
    allocate(zModelLocal(nLyr,nSpolyLocal),stat=err);    if(err/=0)then;message=message//'error allocating zModelLocal'; return;endif
    allocate(soil2model_map(nSpolyLocal),stat=err);      if(err/=0)then;message=trim(message)//'error allocating soil2model_map';return;endif
    do iPoly=1,nSpolyLocal
      allocate(soil2model_map(iPoly)%layer(nLyr),stat=err); if(err/=0)then;message=trim(message)//'error allocating soil2model_map%layer';return;endif
      do iMLyr=1,nLyr
        allocate(soil2model_map(iPoly)%layer(iMLyr)%weight(nSub),stat=err);  if(err/=0)then;message=trim(message)//'error allocating lyrmap%layer%weight';return;endif
        allocate(soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(nSub),stat=err);if(err/=0)then;message=trim(message)//'error allocating lyrmap%layer%ixSubLyr';return;endif
      enddo
    enddo
  ! *****
  ! (3.1) Extract soil poly ID, weight polygon , and soil properties for current model hru 
  ! *********************************************************************
    call subSoilData(sdata, polySub, sdataLocal, err, cmessage)
    if(err/=0)then; message=trim(message)//cmessage; return; endif
    if ( iHru == iHruPrint ) then
      print*,' '
      print*,'****************************************************'
      print*,'HRU, hruID = ',iHru,hruID(iHru)
      print*,'****************************************************'
      print*,'(1.1) Print list of soil polygon ID and weigth'
      write(*,"(' polyID = ',100I9)") (polySub(iPoly), iPoly=1,nSpolyLocal)
      write(*,"(' weight = ',100f9.3)") (swgtSub(iPoly), iPoly=1,nSpolyLocal)
      do iVar=1,nVarSoilData
        write(*,"(A12,'= ')") adjustl(sdataLocal(iVar)%varName)
        do iPoly=1,nSpolyLocal
        select case(trim(sdata_meta(iVar)%vartype))
          case('integer')
            select case(trim(sdata_meta(iVar)%vardims))
              case('2D'); write(*,"(20I9)") (sdataLocal(iVar)%ivar2(iSLyr,iPoly), iSlyr=1,nSlyrs)
              case('1D'); write(*,"(I9)") (sdataLocal(iVar)%ivar1(iPoly))
            end select
          case('double')
            select case(trim(sdata_meta(iVar)%vardims))
              case('2D')
                write(*,"(20f9.3)") (sdataLocal(iVar)%dvar2(iSLyr,iPoly), iSlyr=1,nSlyrs)
              case('1D')
                write(*,"(f9.3)") (sdataLocal(iVar)%dvar1(iPoly))
            end select
          end select
        end do
      enddo
    endif
  ! *****
  ! (3.2) Extract veg poly ID, weight polygon, and veg properties for current model hru 
  ! *********************************************************************
    ! Select list of veg polygon ID contributing a current hru 
    allocate(vmask(nOverVpoly),stat=err); if(err/=0)then;message=message//'error allocating vmask';return;endif
    vmask       = ( overVpolyID(:,iLocal) /= imiss )
    allocate(vPolySub(count(vmask)),stat=err); if(err/=0)then;message=message//'error allocating vPolySub';return;endif
    allocate(vwgtSub(count(vmask)),stat=err); if(err/=0)then;message=message//'error allocating vwgtSub';return;endif
    vPolySub    = pack(overVpolyID(:,iLocal),vmask)
    vwgtSub     = pack(overVpolyID(:,iLocal),vmask)
    nVpolyLocal = size(vPolySub)
    if ( iHru == iHruPrint ) then
      print*,' '
      print*,'(1.3) Print list of vege polygon ID and weigth'
        write(*,"(' polyID = ',100I7)") (vPolySub(iPoly), iPoly=1,nVpolyLocal)
        write(*,"(' weight = ',100f7.3)") (vwgtSub(iPoly), iPoly=1,nVpolyLocal)
    endif
  ! *****
  ! (3.4) Compute model soil parameters using transfer function
  ! *********************************************************
    ! compute model soil parameters
    call comp_soil_model_param(parSxySz, sdataLocal, gammaParMasterMeta, nSlyrs, nSpolyLocal, err, cmessage)
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif 
    if ( iHru == iHruPrint ) then
      print*,'(2) Print Model parameter for polygon and layer'
      write(*,"(' Layer       =',20I9)") (iSLyr, iSlyr=1,nSlyrs)
      do iParm=1,nSoilParModel
        do iPoly = 1,nSpolyLocal
          write(*,"(1X,A12,'=',20f9.3)") betaInGamma(iParm), (parSxySz(iParm)%varData(iSLyr,iPoly), iSlyr=1,nSlyrs)
        enddo
      enddo
    endif
  ! *****
  ! (3.5) Compute Model vege parameters using transfer function
  ! *********************************************************
    ! compute model veg parameters
!      allocate(parVxy(nVegParModel),stat=err)
!      do iParm=1,nVegParModel
!        select case ( vpar_meta(iParm)%dims )
!          case( '1D' , 'ST' )
!            allocate(parVxy(iParm)%varData(nVpolyLocal),stat=err)
!        end select
!      enddo
!      call comp_veg_model_param(parVxy, vprpLocal, vpar_meta, model_name)
!      if ( iHru == iHruPrint ) then
!        print*,'(2) Print Model parameter for Vege polygon'
!        do iparm=1,nVegParModel
!          write(*,"(1X,A10,'= ',100f5.2)") vpar_meta(iparm)%parname, (parVxy(iparm)%varData(iPoly), iPoly=1,nVpolyLocal)
!        enddo
!      endif
  ! ***********
  ! (4) Spatial aggregation of Model parameter 
  ! *********************************************************************
    ! **********
    ! (4.1) Aggregate model parameveter vertical direction - soil data layers to model soil layers
    ! *********************************************************************
    ! (4.1.1) Compute Model layer depth 
    call comp_model_depth(hModelLocal, zModelLocal, hfrac, sdataLocal, err, cmessage) 
    if ( iHru == iHruPrint ) then
      print*, '(3.1.1) Print model depth ---' 
      write(*,"(' Layer       =',20I9)") (iMLyr, iMlyr=1,nlyr)
      do iPoly = 1,nSpolyLocal
        write(*,"('z            =',20f9.3)") (zModelLocal(iMLyr,iPoly), iMlyr=1,nLyr)
      enddo 
    endif
    ! (4.1.2) Computing vertical weight of soil layer for each model layer
    second: associate( hSoilLocal =>sdataLocal(ixVarSoilData%hslyrs)%dvar2 ) 
    call map_slyr2mlyr(hSoilLocal, zModelLocal, soil2model_map, err, cmessage)
    end associate second
    if ( iHru == iHruPrint ) then
      print*, '(3.1.1) Print info on mapping from soil to model layer --' 
      do iPoly = 1,nSpolyLocal
        do iMLyr = 1,nLyr
          print*, 'Poly, Layer=',iPoly,iMLyr
          write(*,"('weightSubLyr= ', 11f9.3)") (soil2model_map(iPoly)%layer(iMLyr)%weight(iSub), iSub=1,nSub)
          write(*,"('indexSubLyr = ', 11I7)")   (soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(iSub), iSub=1,nSub)
        enddo 
      enddo
    endif
    ! (4.1.3) Computing vertical weight for each model layer
    do iPoly=1,nSpolyLocal
      do iMLyr = 1,nLyr
        do iParm = 1,nSoilParModel
          allocate( paramvec(iParm)%layer(nSub),stat=err); if(err/=0)then;message=trim(message)//'error allocating paramvec';return;endif 
          paramvec(iParm)%layer = 0._dp
        enddo
        do iSub=1,nSub
          third: associate(iSelect => soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(iSub) )
          if (iSelect > 0 ) then 
            do iparm = 1,nSoilParModel
              paramvec(iParm)%layer(iSub) = parSxySz(iparm)%varData(iSelect,iPoly)
            enddo
          endif
          end associate third
        enddo
        do iParm = 1,nSoilParModel
          forth: associate( ix=>get_ixPar(trim(betaInGamma(iParm))) )
          if ( trim(gammaParMasterMeta(ix)%vups)/='na' )then
            call aggreg(parSxyMz(iParm)%varData(iMLyr,iPoly),         &
                        soil2model_map(iPoly)%layer(iMLyr)%weight(:), &
                        paramvec(iParm)%layer(:),                     &
                        gammaParMasterMeta(ix)%vups,                  &
                        gammaParMasterMeta(ix)%vpnorm,                &
                        err, cmessage)
             if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
          endif
          end associate forth
        enddo
        do iparm = 1,nSoilParModel
          deallocate(paramvec(iParm)%layer,stat=err);if(err/=0)then;message=trim(message)//'error deallocating paramvec%layer';return;endif
        enddo
      enddo ! end of iMLyr (model layer) loop
    enddo ! end of iPoly
    ! ************
    ! (4.2) Aggregate model parameveter horizontally - use spatial weight netCDF 
    ! *********************************************************************
    do iMLyr=1,nLyr
      do iParm = 1,nSoilParModel
        allocate( paramvec(iParm)%layer(nSpolyLocal),stat=err); if(err/=0)then; message=trim(message)//'error allocating paramvec'; return; endif
        paramvec(iParm)%layer = 0._dp
      enddo
      do iPoly=1,nSpolyLocal
        do iParm = 1,nSoilParModel
            paramvec(iParm)%layer(iPoly) = parSxyMz(iparm)%varData(iMLyr,iPoly)
        enddo
      enddo
      call aggreg(hModel(iMLyr,iHru), swgtsub(:), hModelLocal(iMLyr,:),'pnorm', 1.0_dp, err, cmessage)
      do iParm = 1,nSoilParModel
        fifth: associate( ix=>get_ixPar(trim(betaInGamma(iParm))) )
        if ( trim(gammaParMasterMeta(ix)%hups)/='na' )then
          call aggreg(parMxyMz(iParm)%varData(iMLyr,iHru), &
                      swgtsub(:),                          &
                      paramvec(iparm)%layer(:),            &
                      gammaParMasterMeta(ix)%hups,         &
                      gammaParMasterMeta(ix)%hpnorm,       &
                      err, cmessage)
          if ( iHru == iHruPrint ) then
            print*,'-----------------------------------'
            print*,'Aggregated soil parameter '
            write(*,"(1X,A17,'(layer ',I2,') = ',100f9.3)") gammaParMasterMeta(ix)%pname,iMLyr ,parMxyMz(iParm)%varData(iMLyr,iHru)
          endif
        endif
        end associate fifth
      enddo
      do iParm = 1,nSoilParModel
        deallocate(paramvec(iParm)%layer, stat=err); if(err/=0)then; message=trim(message)//'error deallocating paramvec%layer'; return; endif 
      enddo
    enddo
!      do iparm = 1,nVegParModel
!        if (vpar_meta(iparm)%h_agg) then
!          call aggreg(vegParMxy(iparm)%varData(iHru), &
!                      vwgtSub(:),                     &
!                      parVxy(iparm)%varData,          &
!                      dmiss,                          &
!                      vpar_meta(iparm)%haggmethod)
!          if ( iHru == iHruPrint ) then
!            print*,'-----------------------------------'
!            print*,'Aggregated vege parameter '
!            write(*,"(1X,A17,'= ',100f9.3)") vpar_meta(iparm)%parname, (vegParMxy(iparm)%varData(iHru))
!          endif
!        endif
!      enddo
    ! deallocate memmory
    deallocate(mask,stat=err);           if(err/=0)then; message=trim(message)//'error deallocating mask'; return; endif
    deallocate(polySub,stat=err);        if(err/=0)then; message=trim(message)//'error deallocating polyIdSub'; return; endif
    deallocate(swgtSub,stat=err);        if(err/=0)then; message=trim(message)//'error deallocating swgtSub'; return; endif
    deallocate(vmask,stat=err);          if(err/=0)then; message=trim(message)//'error deallocating vmask'; return; endif
    deallocate(vpolySub,stat=err);       if(err/=0)then; message=trim(message)//'error deallocating vpolySub'; return; endif
    deallocate(vwgtSub,stat=err);        if(err/=0)then; message=trim(message)//'error deallocating vwgtSub'; return; endif
    deallocate(soil2model_map,stat=err); if(err/=0)then; message=trim(message)//'error deallocating soil2model_map';return;endif
    deallocate(hModelLocal,stat=err);    if(err/=0)then; message=trim(message)//'error deallocating hModelLocal';return;endif
    deallocate(zModelLocal,stat=err);    if(err/=0)then; message=trim(message)//'error deallocating zModelLocal';return;endif
    do iParm=1,nSoilParModel
      deallocate(parSxySz(iParm)%varData,stat=err);if(err/=0)then; message=message//'error deallocating parSxySz%varData'; return; endif 
      deallocate(parSxyMz(iparm)%varData,stat=err);if(err/=0)then; message=message//'error deallocating parSxyMz%varData'; return; endif
    enddo
  enddo hru
  end associate
  return
end subroutine
 
! private subroutine:
subroutine subSoilData(soilData, subPolyID, soilDataLocal, err, message)
  use var_lookup,           only:ixVarSoilData,nVarSoilData  ! index of soil data variables and number of variables 
  use globalData,           only:sdata_meta
  implicit none
  !input variables
  type(namevar),        intent(in)    :: soilData(:)      ! soil data container for all the soil polygons
  integer(i4b),         intent(in)    :: subPolyID(:)  ! subset of soil polygon id  
  !output variables
  type(namevar),        intent(inout) :: soilDataLocal(:) ! soil data container for local soil polygon 
  integer(i4b),         intent(out)   :: err           ! error code
  character(*),         intent(out)   :: message       ! error message
  !local variables
  integer(i4b),allocatable            :: polyID(:)
  integer(i4b)                        :: iPoly         ! index of soil polygon loop 
  integer(i4b)                        :: iVar          ! index of named variable loop
  integer(i4b)                        :: iLocal        ! index of hru array in mapping file that match hru id of interest 
  integer(i4b)                        :: nPoly         ! number of polygons in subset 
  integer(i4b)                        :: nSlyrs        ! number of soil layers 
  integer(i4b)                        :: iDummy(1)     ! 1D integer array for temporal storage 

  err=0; message='subSoilData/'
  allocate(polyID(size(soilData(ixVarSoilData%polyid)%ivar1)))
  polyID = soilData(ixVarSoilData%polyid)%ivar1 
  nPoly=size(subPolyID)
  do iVar=1,nVarSoilData
    soilDataLocal(ivar)%varName=trim(soilData(ivar)%varName)
    select case(trim(sdata_meta(iVar)%vartype))
      case('integer')
        select case(trim(sdata_meta(iVar)%vardims))
          case('2D')
            nSlyrs=size(soilData(ivar)%ivar2,1) 
            if ( allocated(soilDataLocal(ivar)%ivar2) ) deallocate(soilDataLocal(ivar)%ivar2)
            allocate(soilDataLocal(ivar)%ivar2(nSlyrs,nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 2D int space for soilDataLocal data structure'; return; endif
          case('1D')
            if ( allocated(soilDataLocal(ivar)%ivar1) ) deallocate(soilDataLocal(ivar)%ivar1)
            allocate(soilDataLocal(ivar)%ivar1(nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 1D int space for soilDataLocal data structure'; return; endif
          end select
      case('double')
        select case(trim(sdata_meta(iVar)%vardims))
          case('2D')
            nSlyrs=size(soilData(ivar)%dvar2,1) 
            if ( allocated(soilDataLocal(ivar)%dvar2) ) deallocate(soilDataLocal(ivar)%dvar2)
            allocate(soilDataLocal(ivar)%dvar2(nSlyrs,nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 2D real space for soilDataLocal data structure'; return; endif
          case('1D')
            if ( allocated(soilDataLocal(ivar)%dvar1) ) deallocate(soilDataLocal(ivar)%dvar1)
            allocate(soilDataLocal(ivar)%dvar1(nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 1D real space for soilDataLocal data structure'; return; endif
        end select
    end select 
    do iPoly=1,nPoly
      if (minval( abs(polyID-subPolyID(iPoly)) ) /= 0 )then; err=10; message=trim(message)//'hru id does not exist in mapping file';return; endif
      iDummy = minloc( abs(polyID-subPolyID(iPoly)) )
      iLocal=iDummy(1)
      select case(trim(sdata_meta(iVar)%vartype))
        case('integer')
          select case(trim(sdata_meta(iVar)%vardims))
            case('2D'); soilDataLocal(ivar)%ivar2(:,iPoly) = soilData(ivar)%ivar2(:,iLocal)
            case('1D'); soilDataLocal(ivar)%ivar1(iPoly)   = soilData(ivar)%ivar1(iLocal)
            end select
        case('double')
          select case(trim(sdata_meta(iVar)%vardims))
            case('2D'); soilDataLocal(ivar)%dvar2(:,iPoly) = soilData(ivar)%dvar2(:,iLocal)
            case('1D'); soilDataLocal(ivar)%dvar1(iPoly)   = soilData(ivar)%dvar1(iLocal)
          end select
      end select 
    end do
  end do 
  return 
end subroutine 

! private subroutine:
subroutine pop_hfrac(gammaParStr, gammaParMeta,hfrac, err, message)
  use globalData,   only: gammaSubset
  implicit none
  !input variables
  type(var_d),          intent(in)  :: gammaParStr(:)
  type(cpar_meta),      intent(in)  :: gammaParMeta(:)
  !output variables
  real(dp),             intent(out) :: hfrac(:)
  integer(i4b),         intent(out) :: err         ! error code
  character(*),         intent(out) :: message     ! error message
  !local variables
  real(dp)                          :: dummy(20) 
  logical(lgc)                      :: mask(20) 
  integer(i4b)                      :: i

  ! initialize error control
  err=0; message='pop_hfrac/'
  dummy=-999
  !check h parameters - now can chcek up to 5 layers
  do i=1,size(gammaSubset)
    if (gammaParMeta(i)%pname=="h1gamma1")then;dummy(1)=gammaParStr(i)%var(1);cycle;endif 
    if (gammaParMeta(i)%pname=="h2gamma1")then;dummy(2)=gammaParStr(i)%var(1);cycle;endif
    if (gammaParMeta(i)%pname=="h3gamma1")then;dummy(3)=gammaParStr(i)%var(1);cycle;endif
    if (gammaParMeta(i)%pname=="h4gamma1")then;dummy(4)=gammaParStr(i)%var(1);cycle;endif
  enddo
  mask=(dummy>0)
  if ( count(mask)/=nLyr-1 ) stop 'number of h1gamma prameters mismatch with nLyr'
  hfrac=pack(dummy,mask)
  return
end subroutine 

character(len=20) function int2str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (int2str, *) k
    int2str = adjustl(int2str)
end function

end module mpr_routine
