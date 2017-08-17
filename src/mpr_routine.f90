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
! this subroutine is used for opt=3 in namelist (run only mpr and output parameters)
subroutine run_mpr( calParam, restartFile, err, message ) 
  use globalData,    only: calScaleMeta, calParMeta, calGammaMeta, parArray, parMask, nCalPar, nSoilBetaModel, nVegBetaModel
  use model_wrapper, only: read_hru_id
  use write_param_nc,only: write_nc_soil,write_nc_veg
  implicit none
  ! input variables
  real(dp),             intent(in)  :: calParam(:)                   ! parameter in namelist, not necessarily all parameters are calibrated
  character(len=strLen),intent(in)  :: restartFile                   ! name of restart file including iteration, the most recent parameter values 
  ! output variables
  integer(i4b),         intent(out) :: err                           ! error id 
  character(len=strLen),intent(out) :: message                       ! error message
  ! local
  real(dp),             allocatable :: params(:)                     ! parameter vector that is input into mpr
  type(var_d)                       :: calParStr(nCalPar)            ! parameter storage including perLayr values converted from parameter array 
  type(var_d)                       :: pnormCoef(size(calScaleMeta)) ! parameter storage converted from parameter array 
  type(var_d),          allocatable :: paramGammaStr(:)              ! calibratin gamma parameter storage extracted from calParStr
  integer(i4b)                      :: idx                           ! counter 
  integer(i4b)                      :: iPar                          ! loop index for parameter 
  logical(lgc)                      :: mask(nCalPar)                 ! 1D mask
  integer(i4b)                      :: hruID(nHru)                   ! Hru ID
  real(dp)                          :: hModel(nLyr,nHru)             ! storage of model layer thickness at model layer x model hru 
  type(namedvar2)                   :: parMxyMz(nSoilBetaModel)      ! storage of model soil parameter at model layer x model hru 
  type(namedvar2)                   :: vegParMxy(nVegBetaModel)      ! storage of model vege parameter at month x model hru
  character(len=strLen)             :: cmessage                      ! error message from subroutine

  err=0; message='run_mpr/' ! to initialize error control
  if ( any(calParMeta(:)%beta /= "beta") )then ! calPar need to include min. one gamma parameter to be used for MPR 
    if ( idModel/=0 )then;idModel=0;print*,trim(message)//'idModel is set to zero - model inepenedent';endif
    allocate(params, source=calParam) ! copy calParameter default
    call restartIO( params, restartFile ) ! check restartFile exist, if it does, update params vector, otherwise wite them in restartFile
    ! transform parameter vector to custom data type - calParStr and pnormCoef
    idx=1
    do iPar=1,nCalPar !beta and gamma parameter values
      if (calParMeta(iPar)%perLyr)then
        allocate(calParStr(iPar)%var(nLyr))
        calParStr(iPar)%var=params(idx:idx+nLyr-1)
        idx=idx+nLyr
      else
        allocate(calParStr(iPar)%var(1))
        calParStr(iPar)%var=params(idx)
        idx=idx+1
      endif
    end do
    do iPar=1,size(calScaleMeta) ! scaling parameter for beta parameter estimated with MPR
      allocate(pnormCoef(iPar)%var(2))
      pnormCoef(iPar)%var=params(idx:idx+1)
      idx=idx+2
    end do
    ! Get hruID from mapping file
    call read_hru_id(idModel, hruID, err, cmessage)
    if (err/=0)then;message=trim(message)//trim(cmessage);return;endif
    do iPar=1,nSoilBetaModel
      allocate(parMxyMz(iPar)%varData(nLyr,nHru),stat=err)
    enddo
    do iPar=1,nVegBetaModel
      allocate(vegParMxy(iPar)%varData(nMonth,nHru),stat=err)
    enddo
    mask=calParMeta(:)%beta/="beta"
    allocate(paramGammaStr(count(mask)))
    paramGammaStr=pack(calParStr,mask)
    do iPar=1,size(paramGammaStr)
      if (size(paramGammaStr(iPar)%var)>1)then;message=trim(message)//'gammaParameter should not have perLayer value';return;endif
    enddo
    !perform MPR
    call mpr(hruID, pnormCoef, paramGammaStr, calGammaMeta, hModel, parMxyMz, vegParMxy, err, cmessage) ! to output model layer thickness and model parameter via MPR
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    !! Write parameter derived from MPR in netCDF 
    if (nSoilBetaModel>0_i4b)then
      call write_nc_soil(trim(mpr_output_dir)//trim(soil_param_nc), hruID, hModel, parMxyMz, err, cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    endif
    if (nVegBetaModel>0_i4b)then
      call write_nc_veg(trim(mpr_output_dir)//trim(veg_Param_nc), hruID, vegParMxy, err, cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
    endif
  else  
    print*,trim(message)//'there is no gamma pamameters in parameter input file to perform MPR';stop
  endif
  return

  contains

  subroutine restartIO( params, restartFile )
    implicit none
    character(len=strLen), intent(in)    :: restartFile
    real(dp),              intent(inout) :: params(:)
    logical                              :: isExistFile ! logical to check if the file exist or not
    integer(i4b)                         :: i,j         ! loop index for parameter 
    integer(i4b)                         :: idummy      ! dummy integer vaiable
    logical(lgc)                         :: ldummy      ! dummy logical vaiable
    character(len=strLen)                :: cdummy      ! dummy character vaiable

    inquire(file=trim(restartFile), exist=isExistFile)
    if ( isExistFile ) then !  if state file exists, read it and update params
      print*, 'read restart file'
      open(unit=70,file=trim(adjustl(restartFile)), action='read', status = 'unknown')
      do i=1,nCalPar
        if (calParMeta(i)%perLyr)then
          do j=1,nLyr
            read(70,*) cdummy, params(i), ldummy, cdummy
          end do
        else
          read(70,*) cdummy, params(i), ldummy, cdummy
        endif
      end do
      do i=1,size(calScaleMeta) 
         read(70,*) cdummy, params(nCalPar+2*i-1), ldummy, cdummy
         read(70,*) cdummy, params(nCalPar+2*i), ldummy, cdummy
      end do
      close(70)
    else          !otherwise write out
      print*, 'write restart file'
      open(unit=70,file=trim(adjustl(restartFile)), action='write', status = 'unknown')
      do i=1,nCalPar
        if (calParMeta(i)%perLyr)then
          do j=1,nLyr
            write(70,100) calParMeta(i)%pname(1:20), parArray(i,1), parMask(i), 'Beta-par-perLayer'
            100 format(1X,A,1X,ES17.10,1X,L9,1X,A20)
          end do
        else
          write(70,200) calParMeta(i)%pname(1:20), parArray(i,1), parMask(i), 'Gamma-par'
          200 format(1X,A,1X,ES17.10,1X,L9,1X,A20)
        endif
      enddo  
      do i=1,size(calScaleMeta) 
         write(70,300) calScaleMeta(i)%betaname(1:20), parArray(nCalPar+2*i-1,1), parMask(nCalPar+2*i-1), 'H-scaling-par'
         write(70,300) calScaleMeta(i)%betaname(1:20), parArray(nCalPar+2*i  ,1), parMask(nCalPar+2*i),   'V-scaling-par'
         300 format(1X,A,1X,ES17.10,1X,L9,1X,A20)
      end do
      close(70)
    endif
    return
  end subroutine

end subroutine

! ************************************************************************************************
! public subroutine: mpr 
! ************************************************************************************************
subroutine mpr(hruID,             &     ! input: hruID
               pnormCoefStr,      &     ! input: list of pnorm coefficients
               gammaParStr,       &     ! input: array of gamma parameter (including h and z)
               gammaParMeta,      &     ! input: array of gamma parameter metadata (including h and z)
               hModel,            &     ! output: Model layer thickness
               parMxyMz,          &     ! output: MPR derived soil parameter
               vegParMxy,         &     ! output: MPR derived veg parameter
               err, message)            ! output: error id and message
  use model_wrapper,        only:read_hru_id
  use popMeta,              only:popMprMeta
  use globalData,           only:betaMeta, gammaMeta, soilBetaCalName, vegBetaCalName, calScaleMeta, &
                                 sdata_meta, vdata_meta, map_meta, nSoilBetaModel, nVegBetaModel
  use get_ixname,           only:get_ixGamma, get_ixBeta
  use tf,                   only:comp_model_param              ! Including Soil model parameter transfer function
  use modelLayer,           only:comp_model_depth              ! Including model layr depth computation routine 
  use modelLayer,           only:map_slyr2mlyr                 ! Including model layr computation routine 
  use upscaling,            only:aggreg                        ! Including Upscaling operator 
  use read_mapdata,         only:getMapData                    ! routine to read mapping data into data structures
  use read_vegdata,         only:getVegData                    ! routine to read veg data into data structures
  use read_vegdata,         only:getVegClassLookup             ! routine to read veg calss-property lookupu table 
  use read_soildata,        only:getData                       ! routine to read soil data into data structures
  use read_soildata,        only:mod_hslyrs                    ! routine to modify soil layer thickness and updata soil data structure
  use var_lookup,           only:ixBeta, ixGamma               ! 
  use var_lookup,           only:ixVarSoilData,nVarSoilData    ! index of soil data variables and number of variables 
  use var_lookup,           only:ixVarVegData,nVarVegData      ! index of vege data variables and number of variables
  use var_lookup,           only:ixVarMapData,nVarMapData      ! index of map data variables and number of variables
  use var_lookup,           only:ixPrpVeg,nPrpVeg              ! index of veg properties and number of properties

  implicit none

  ! input
  integer(i4b),         intent(in)   :: hruID(:)                   ! hruID list
  type(var_d),          intent(in)   :: gammaParStr(:)             ! data structure of gamma parameter value adjusted with calibration
  type(var_d),          intent(in)   :: pnormCoefStr(:)            ! data structure of pnorm coefficient value adjusted with calibration
  type(cpar_meta),      intent(in)   :: gammaParMeta(:)            ! array of calibrating gamma parameter meta data
  ! output
  real(dp),             intent(out)  :: hModel(:,:)                ! Model layer thickness at model layer x model hru 
  type(namedvar2),      intent(inout):: parMxyMz(:)                ! storage of model soil parameter at model layer x model hru 
  type(namedvar2),      intent(inout):: vegParMxy(:)               ! storage of model vege parameter at model hru
  integer(i4b),         intent(out)  :: err                        ! error code
  character(len=strLen),intent(out)  :: message                    ! error message 
  ! local
  character(len=strLen)              :: cmessage                   ! error message from downward subroutine
  integer(i2b),parameter             :: iHruPrint = -1              ! model hru id (index) for which everything is printed for checking
  integer(i4b),parameter             :: nSub=11                    ! max. number of Soil layer within Model layer
  integer(i4b)                       :: iLocal                     ! index of hru array in mapping file that match hru id of interest 
  integer(i4b)                       :: iGamma                     ! index loop
  integer(i4b)                       :: iPoly                      ! Loop index of soil polygon
  integer(i4b)                       :: iSLyr                      ! Loop index of soil layer 
  integer(i4b)                       :: iMLyr                      ! loop index of model soil layer (1,...,n from top to bottom) 
  integer(i4b)                       :: iParm                      ! Loop index of model parameters (e.g., VIC)
  integer(i4b)                       :: iVar                       ! Loop index of miscleneous variables 
  integer(i4b)                       :: iMon                       ! Loop index of month 
  integer(i4b)                       :: iHru                       ! loop index of hrus 
  integer(i4b)                       :: iSub                       ! Loop index of multiple soi layers in model layer
  integer(i4b)                       :: ixStart                    ! starting index of subset geophysical polygons in the entire dataset 
  integer(i4b)                       :: ixEnd                      ! ending index of subset geophysical polygons in the entire dataset 
  type(par_meta),allocatable         :: gammaUpdateMeta(:)
  type(par_meta),allocatable         :: betaUpdateMeta(:)
  integer(i4b)                       :: nSpoly                     ! number of soil polygon in entire soil data domain 
  integer(i4b)                       :: nSlyrs                     ! number of soil layers
  integer(i4b)                       :: nVpoly                     ! TO BE DELETED: number of vege polygon (grid box) in entire vege data domain 
  integer(i4b)                       :: nShru                      ! number of hrus in soil mapping file) 
  integer(i4b)                       :: nVhru                      ! TO BE DELETED: number of hrus in vege mapping file (should be equal to nVhru)
  integer(i4b)                       :: nOverSpoly                 ! number of overlapped soil polygon for each model HRU 
  integer(i4b)                       :: nOverVpoly                 ! TO BE DELETED: number of overlapped vege polygon for each model HRU 
  integer(i4b)                       :: nSpolyLocal                ! number of subset overlapped soil polygon for each model HRU 
  integer(i4b)                       :: nVpolyLocal                ! TO BE DELETED: number of subset overlapped vege polygon for each model HRU 
  real(dp)                           :: hfrac(nLyr-1)              ! fraction of soil depth for each model layer 
  type(namevar)                      :: sdata(nVarSoilData)        ! soil data container for all the soil polygons
  type(namevar)                      :: sdataLocal(nVarSoilData)   ! soil data container for local soil polygon 
  type(namevar)                      :: vdata(nVarVegData)         ! veg data container for all the veg polygons
  type(namevar)                      :: vdataLocal(nVarVegData)    ! veg data container for local vege polygon 
  integer(i4b)                       :: vegClass(nVclass)          ! veg class array (e.g., IGBP)
  type(var_d)                        :: vcls2prp(nVclass)          ! storage of property value for each veg class
  integer(i4b)                       :: iVclass                    ! ID (= index) of vege class 
  real(dp),        allocatable       :: hModelLocal(:,:)           ! Model layer thickness for soil polygon within one hru
  real(dp),        allocatable       :: zModelLocal(:,:)           ! Model layer depth for soil polygon within one hru
  type(namedvar2)                    :: parSxySz(nSoilBetaModel)   ! storage of model soil parameter for 2D field -soil layer x soil poy 
  type(namedvar2)                    :: parSxyMz(nSoilBetaModel)   ! storage of model soil parameter for 2D field -model layer x soil poy
  type(namedvar2)                    :: parVxy(nVegBetaModel)      ! storage of model vege parameter for 1D or 2D field - vege poly (x month)
  integer(i4b),    allocatable       :: polySub(:)                 ! list of ID (=index) of soil polygons contributing model hru
  integer(i4b),    allocatable       :: vPolySub(:)                ! TO BE DELETED: list of ID (=index) of veg polygons contributing model hru 
  type(mapvar)                       :: mapdata(2)                 ! map data container for all the soil polygons
  real(dp),        allocatable       :: swgtSub(:)                 ! adjusted Areal weight of soil polygon for all model hrus
  real(dp),        allocatable       :: vwgtSub(:)                 ! TO BE DELETED: adjusted Areal weight of veg polygon for one model hrus
  type(poly),      allocatable       :: soil2model_map(:)          ! data structure to hold weight and index of contributing soil layer per model layer and per soil polygon
  type(var_d)                        :: soilParVec(nSoilBetaModel) ! data structure to hold soil parameter vector e.g., all layers for one polygon
  type(var_d)                        :: vegParVec(nVegBetaModel)    ! data structure to hold veg parameter vector e.g., all months for one polygon
  
  ! initialize error control
  err=0; message='mpr/'
  !(0) Preparation
  allocate(gammaUpdateMeta, source=gammaMeta) ! copy master gamma parameter metadata
  allocate(betaUpdateMeta,  source=betaMeta)  ! copy master beta parameter metadata
  ! Swap gammaUpdateMeta%val with adjusted gammaPar value
  do iGamma=1,size(gammaParStr)
    gammaUpdateMeta(gammaParMeta(iGamma)%ixMaster)%val=gammaParStr(iGamma)%var(1)
  enddo
  ! Swap betaUpdateMeta%hpnorm and vpnorm with adjusted pnorm coefficient value
  do iParm=1,size(pnormCoefStr)
    associate(ix=>get_ixBeta(calScaleMeta(iParm)%betaname))
      betaUpdateMeta(ix)%hpnorm=pnormCoefStr(iParm)%var(1)
      betaUpdateMeta(ix)%vpnorm=pnormCoefStr(iParm)%var(2)
    end associate
  enddo 
  call popMprMeta( err, cmessage)   !for sdata_meta, vdata_meta, map_meta
  if(err/=0)then; message=trim(message)//cmessage; return; endif
  !call pop_hfrac(gammaParStr, gammaParMeta, hfrac, err, cmessage) ! to get hfrac 
  call pop_hfrac( gammaUpdateMeta, hfrac, err, cmessage) ! to get hfrac 
  if(err/=0)then; message=trim(message)//cmessage; return; endif
  ! (1) Get Geophysical data 
  ! *****
  ! (1.1) Get soil data  
  ! *********************************************
  call getData(trim(mpr_input_dir)//trim(fname_soil),& ! input: soil data input name (netCDF)
               sdata_meta,                           & ! input: soil data meta
               dname_spoly,                          & ! input: spatial dimension (polygon ID)
               dname_slyrs,                          & ! input: spatial dimension (polygon ID)
               sdata,                                & ! input-output: soil data structure
               nSpoly,                               & ! output: number of dimension (i.e. number of soil polygon)
               nSlyrs,                               & ! output: number of dimension (i.e. number of soil layer)
               err, cmessage)
  if(err/=0)then; message=trim(message)//cmessage; return; endif 
  call mod_hslyrs(sdata, gammaUpdateMeta(ixGamma%z1gamma1)%val, err,cmessage) ! modify soil layer thickness in sdata data structure
  if(err/=0)then; message=trim(message)//cmessage; return; endif 
  ! *****
  ! (1.2)  veg class - properties lookup table ...
  ! *********************************************
  ! (1.2.1) Read in veg data netCDF...
   do iVclass=1,nVclass
     allocate(vcls2prp(iVclass)%var(nPrpVeg),stat=err);if(err/=0)then;message=trim(message)//'error allocating vcls2prp%var';return;endif
   enddo
   call getVegData(trim(mpr_input_dir)//trim(fname_veg), & ! input: veg data input name (netCDF) 
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
  if ( opt==3 .and. nHru /= nShru )then;err=10;message=trim(message)//'nHru= '//trim(int2str(nShru))//' NOT '//trim(int2str(nHru));return;endif
  associate( hruMap        => mapdata(1)%var(ixVarMapData%hru_id)%ivar1(1:nShru),  &
             hruMapVeg     => mapdata(2)%var(ixVarMapData%hru_id)%ivar1(1:nVhru),  &
             swgt          => mapdata(1)%var(ixVarMapData%weight)%dvar1,           &
             vwgt          => mapdata(2)%var(ixVarMapData%weight)%dvar1,           &
             nSoilOverPoly => mapdata(1)%var(ixVarMapData%overlaps)%ivar1,         & 
             nVegOverPoly  => mapdata(2)%var(ixVarMapData%overlaps)%ivar1,         & 
             overVpolyID   => mapdata(2)%var(ixVarMapData%overlapPolyId)%ivar1,    & 
             overSpolyID   => mapdata(1)%var(ixVarMapData%overlapPolyId)%ivar1 )
  if ( minval(abs(hruMap-hruMapVeg)) /= 0 )then;err=11;message=trim(message)//'different hru id in vege and soil mapping file';return;endif
  !!! ---------------------------------------------
  !!! Start of model hru loop (from mapping file) !!!
  !!! ---------------------------------------------
  hru_loop: do iHru=1,nHru
    ! Get index (iLocal) of hru id that matches with current hru from hru id array in mapping file
    call findix(hruID(iHru), hruMap, iLocal, err, cmessage); if(err/=0)then; message=trim(cmessage)//cmessage; return; endif
    ! Select list of soil polygons contributing a current hru
    ixEnd       = sum(nSoilOverPoly(1:iLocal));
    ixStart     = ixEnd-nSoilOverPoly(iLocal)+1
    nSpolyLocal = nSoilOverPoly(iLocal)
    allocate(polySub(nSpolyLocal),stat=err); if(err/=0)then;message=message//'error allocating polySub';return;endif
    allocate(swgtSub(nSpolyLocal),stat=err); if(err/=0)then;message=message//'error allocating swgtSub';return;endif
    polySub = overSpolyId(ixStart:ixEnd) ! id of soil polygons contributing to current hru
    swgtSub = swgt(ixStart:ixEnd)        ! weight of soil polygons contributing to current hru
    ! Select list of veg polygon ID contributing a current hru 
    ixEnd       = sum(nVegOverPoly(1:iLocal));
    ixStart     = ixEnd-nVegOverPoly(iLocal)+1
    nVpolyLocal = nVegOverPoly(iLocal)
    allocate(vPolySub(nVpolyLocal),stat=err); if(err/=0)then;message=message//'error allocating vPolySub';return;endif
    allocate(vwgtSub(nVpolyLocal),stat=err);  if(err/=0)then;message=message//'error allocating vwgtSub';return;endif
    vPolySub = overVpolyID(ixStart:ixEnd) ! id of soil polygons contributing to current hru
    vwgtSub  = vwgt(ixStart:ixEnd)        ! weight of soil polygons contributing to current hru
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
    ! allocate memmory
    if (nSoilBetaModel>0_i4b)then !if at least one soil parameter is included 
      do iParm=1,nSoilBetaModel
        allocate(parSxySz(iParm)%varData(nSlyrs,nSpolyLocal),stat=err); if(err/=0)then;message=message//'error allocating parSxySz%varData';return;endif 
        allocate(parSxyMz(iParm)%varData(nLyr,nSpolyLocal),stat=err);   if(err/=0)then;message=message//'error allocating parSxyMz%varData';return;endif
      enddo
    endif
    if (nVegBetaModel>0_i4b)then !if at least one veg parameter is included
      do iParm=1,nVegBetaModel
        allocate(parVxy(iParm)%varData(nMonth, nVpolyLocal),stat=err)
      enddo
    endif
  ! *****
  ! (3.1) Extract soil data for current model hru 
  ! *********************************************************************
    call subsetData(sdata, polySub, sdataLocal, 'soil' ,err, cmessage)
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
  ! (3.2) Extract vege data for current model hru 
  ! *********************************************************************
    call subsetData(vdata, vpolySub, vdataLocal, 'veg' ,err, cmessage)
    if(err/=0)then; message=trim(message)//cmessage; return; endif
    if ( iHru == iHruPrint ) then
      print*,' '
      print*,'(1.3) Print list of vege polygon ID and weigth'
        write(*,"(' polyID = ',100I9)") (vPolySub(iPoly), iPoly=1,nVpolyLocal)
        write(*,"(' weight = ',100f9.3)") (vwgtSub(iPoly), iPoly=1,nVpolyLocal)
    endif
  ! *****
  ! (3.4) Compute model parameters using transfer function
  ! *********************************************************
    ! compute model soil parameters
    if (nSoilBetaModel>0_i4b)then
      call comp_model_param(parSxySz, parVxy, sdataLocal, vdataLocal, gammaUpdateMeta, nSlyrs, nSpolyLocal, nVpolyLocal, err, cmessage)
      if(err/=0)then;message=trim(message)//trim(cmessage);return;endif 
      if ( iHru == iHruPrint ) then
        print*,'(2) Print Model parameter at native resolution'
        if (nSoilBetaModel>0_i4b)then
          write(*,"(' Layer       =',20I9)") (iSLyr, iSlyr=1,nSlyrs)
          do iParm=1,nSoilBetaModel
            do iPoly = 1,nSpolyLocal
              write(*,"(1X,A12,'=',20f9.3)") soilBetaCalName(iParm), (parSxySz(iParm)%varData(iSLyr,iPoly), iSlyr=1,nSlyrs)
            enddo
          enddo
        endif
        if (nVegBetaModel>0_i4b)then
          do iParm=1,nVegBetaModel
            do iPoly = 1,nVpolyLocal
              write(*,"(1X,A10,'= ',100f9.3)") vegBetaCalName(iParm), (parVxy(iParm)%varData(iMon, iPoly), iMon=1,nMonth)
            enddo
          enddo
        endif 
      endif
    end if
    ! (4.1.1) Compute Model layer depth 
    call comp_model_depth(hModelLocal, zModelLocal, hfrac, sdataLocal, err, cmessage) 
    if(err/=0)then;message=trim(message)//trim(cmessage);return;endif 
    if ( iHru == iHruPrint ) then
      print*, '(3.1.1) Print model depth ---' 
      write(*,"(' Layer       =',20I9)") (iMLyr, iMlyr=1,nlyr)
      do iPoly = 1,nSpolyLocal
        write(*,"('z            =',20f9.3)") (zModelLocal(iMLyr,iPoly), iMlyr=1,nLyr)
      enddo 
    endif
  ! ***********
  ! (4) Spatial aggregation of Model parameter 
  ! *********************************************************************
    if (nSoilBetaModel>0_i4b)then
      ! **********
      ! (4.1) Aggregate model parameveter vertical direction - soil data layers to model soil layers
      ! *********************************************************************
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
          do iParm = 1,nSoilBetaModel
            allocate( soilParVec(iParm)%var(nSub),stat=err); if(err/=0)then;message=trim(message)//'error allocating soilParVec';return;endif 
            soilParVec(iParm)%var = 0._dp
          enddo
          do iSub=1,nSub
            third: associate(iSelect => soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(iSub) )
            if (iSelect > 0 ) then 
              do iParm = 1,nSoilBetaModel
                soilParVec(iParm)%var(iSub) = parSxySz(iParm)%varData(iSelect,iPoly)
              enddo
            endif
            end associate third
          enddo
          do iParm = 1,nSoilBetaModel
            forth: associate( ix=>get_ixBeta(trim(soilBetaCalName(iParm))) )
            if ( trim(betaUpdateMeta(ix)%vups)/='na' )then
              call aggreg(parSxyMz(iParm)%varData(iMLyr,iPoly),         &
                          soil2model_map(iPoly)%layer(iMLyr)%weight(:), &
                          soilParVec(iParm)%var(:),                     &
                          betaUpdateMeta(ix)%vups,                      &
                          betaUpdateMeta(ix)%vpnorm,                    &
                          err, cmessage)
               if(err/=0)then;message=trim(message)//trim(cmessage);return;endif
            endif
            end associate forth
          enddo
          do iParm = 1,nSoilBetaModel
            deallocate(soilParVec(iParm)%var,stat=err);if(err/=0)then;message=trim(message)//'error deallocating soilParVec%layer';return;endif
          enddo
        enddo ! end of iMLyr (model layer) loop
      enddo ! end of iPoly
    endif 
    ! ************
    ! (4.2) Aggregate model parameveter horizontally - use spatial weight netCDF 
    ! *********************************************************************
    do iMLyr=1,nLyr
      call aggreg(hModel(iMLyr,iHru), swgtsub(:), hModelLocal(iMLyr,:),'pnorm', 1.0_dp, err, cmessage)
    end do
    if (nSoilBetaModel>0_i4b)then
      do iMLyr=1,nLyr
        do iParm = 1,nSoilBetaModel
          allocate( soilParVec(iParm)%var(nSpolyLocal),stat=err); if(err/=0)then; message=trim(message)//'error allocating soilParVec'; return; endif
          soilParVec(iParm)%var = 0._dp
        enddo
        do iParm = 1,nSoilBetaModel
          do iPoly=1,nSpolyLocal
              soilParVec(iParm)%var(iPoly) = parSxyMz(iParm)%varData(iMLyr,iPoly)
          enddo
        enddo
        do iParm = 1,nSoilBetaModel
          fifth: associate( ix=>get_ixBeta(trim(soilBetaCalName(iParm))) )
          if ( trim(betaUpdateMeta(ix)%hups)/='na' )then
            call aggreg(parMxyMz(iParm)%varData(iMLyr,iHru), &
                        swgtsub(:),                          &
                        soilParVec(iParm)%var(:),            &
                        betaUpdateMeta(ix)%hups,             &
                        betaUpdateMeta(ix)%hpnorm,           &
                        err, cmessage)
            if ( iHru == iHruPrint ) then
              print*,'-----------------------------------'
              print*,'Aggregated soil parameter '
              write(*,"(1X,A17,'(layer ',I2,') = ',100f9.3)") betaUpdateMeta(ix)%pname,iMLyr ,parMxyMz(iParm)%varData(iMLyr,iHru)
            endif
          endif
          end associate fifth
        enddo
        do iParm = 1,nSoilBetaModel
          deallocate(soilParVec(iParm)%var, stat=err); if(err/=0)then; message=trim(message)//'error deallocating soilParVec%var'; return; endif 
        enddo
      enddo
      do iParm=1,nSoilBetaModel
        deallocate(parSxySz(iParm)%varData,stat=err);if(err/=0)then; message=message//'error deallocating parSxySz%varData'; return; endif 
        deallocate(parSxyMz(iParm)%varData,stat=err);if(err/=0)then; message=message//'error deallocating parSxyMz%varData'; return; endif
      enddo
    endif
    ! deallocate memmory
    deallocate(soil2model_map,stat=err); if(err/=0)then; message=trim(message)//'error deallocating soil2model_map';return;endif
    deallocate(hModelLocal,stat=err);    if(err/=0)then; message=trim(message)//'error deallocating hModelLocal';return;endif
    deallocate(zModelLocal,stat=err);    if(err/=0)then; message=trim(message)//'error deallocating zModelLocal';return;endif
    if (nVegBetaModel>0_i4b)then
      do iMon=1,nMonth
        do iParm = 1,nVegBetaModel
          allocate( vegParVec(iParm)%var(nVpolyLocal),stat=err); if(err/=0)then; message=trim(message)//'error allocating vegParVec%var'; return; endif
        enddo
        do iParm = 1,nVegBetaModel
          do iPoly=1,nVpolyLocal
              vegParVec(iParm)%var(iPoly) = parVxy(iParm)%varData(iMon, iPoly)
          enddo
        enddo
        do iParm = 1,nVegBetaModel
          sixth: associate( ix=>get_ixBeta(trim(vegBetaCalName(iParm))) )
          call aggreg(vegParMxy(iParm)%varData(iMon, iHru), &
                      vwgtSub(:),                           &
                      vegParVec(iParm)%var,                 &
                      betaUpdateMeta(ix)%hups,              &
                      betaUpdateMeta(ix)%hpnorm,            &
                      err,cmessage)
          if ( iHru == iHruPrint ) then
            print*,'-----------------------------------'
            print*,'Aggregated vege parameter '
            write(*,"(1X,A17,'(Month ',I2,') = ',100f9.3)") betaUpdateMeta(ix)%pname, iMon ,parMxyMz(iParm)%varData(iMon,iHru)
          endif
          end associate sixth
        enddo
        do iParm = 1,nVegBetaModel
          deallocate(vegParVec(iParm)%var, stat=err); if(err/=0)then; message=trim(message)//'error deallocating vegParVec%var'; return; endif 
        enddo
      enddo
      do iParm=1,nVegBetaModel
        deallocate(parVxy(iParm)%varData,stat=err);if(err/=0)then; message=message//'error deallocating parVxy%varData'; return; endif 
      enddo
    endif
    ! deallocate memmory
    deallocate(polySub,stat=err);        if(err/=0)then; message=trim(message)//'error deallocating polyIdSub'; return; endif
    deallocate(swgtSub,stat=err);        if(err/=0)then; message=trim(message)//'error deallocating swgtSub'; return; endif
    deallocate(vpolySub,stat=err);       if(err/=0)then; message=trim(message)//'error deallocating vpolySub'; return; endif
    deallocate(vwgtSub,stat=err);        if(err/=0)then; message=trim(message)//'error deallocating vwgtSub'; return; endif
  enddo hru_loop
  end associate
  return
end subroutine
 
! private subroutine:
subroutine subsetData(entireData, subPolyID, localData, data_name, err, message)
  use var_lookup,           only:ixVarSoilData,nVarSoilData  ! index of soil data variables and number of variables 
  use var_lookup,           only:ixVarVegData,nVarVegData    ! index of vege data variables and number of variables 
  use globalData,           only:sdata_meta,vdata_meta
  implicit none
  !input variables
  type(namevar),        intent(in)    :: entireData(:) ! soil data container for all the soil polygons
  integer(i4b),         intent(in)    :: subPolyID(:)  ! subset of soil polygon id  
  character(*),         intent(in)    :: data_name     ! data_name e.g., soil, vege 
  !output variables
  type(namevar),        intent(inout) :: localData(:) ! soil data container for local soil polygon 
  integer(i4b),         intent(out)   :: err           ! error code
  character(*),         intent(out)   :: message       ! error message
  !local variables
  type(var_meta), allocatable         :: data_meta(:)
  integer(i4b),allocatable            :: polyID(:)
  integer(i4b)                        :: polyIdx       ! index of polyid variable in dataset 
  integer(i4b)                        :: nVarData      ! number of data 
  integer(i4b)                        :: iPoly         ! index of soil polygon loop 
  integer(i4b)                        :: iVar          ! index of named variable loop
  integer(i4b)                        :: iLocal        ! index of hru array in mapping file that match hru id of interest 
  integer(i4b)                        :: nPoly         ! number of polygons in subset 
  integer(i4b)                        :: nDim1         ! number of 1st dimension 

  err=0; message='subsetData/'
  select case(trim(data_name))
    case('soil');nVarData=nVarSoilData; allocate(data_meta,source=sdata_meta); polyIdx=ixVarSoilData%polyid
    case('veg'); nVarData=nVarVegData;  allocate(data_meta,source=vdata_meta); polyIdx=ixVarVegData%polyid
  end select 

  allocate(polyID(size(entireData(polyIdx)%ivar1)))
  polyID = entireData(polyIdx)%ivar1 

  nPoly=size(subPolyID)
  do iVar=1,nVarData
    localData(ivar)%varName=trim(entireData(ivar)%varName)
    select case(trim(data_meta(iVar)%vartype))
      case('integer')
        select case(trim(data_meta(iVar)%vardims))
          case('2D')
            nDim1=size(entireData(ivar)%ivar2,1) 
            if ( allocated(localData(ivar)%ivar2) ) deallocate(localData(ivar)%ivar2)
            allocate(localData(ivar)%ivar2(nDim1,nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 2D int space for localData data structure'; return; endif
          case('1D')
            if ( allocated(localData(ivar)%ivar1) ) deallocate(localData(ivar)%ivar1)
            allocate(localData(ivar)%ivar1(nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 1D int space for localData data structure'; return; endif
          end select
      case('double')
        select case(trim(data_meta(iVar)%vardims))
          case('2D')
            nDim1=size(entireData(ivar)%dvar2,1) 
            if ( allocated(localData(ivar)%dvar2) ) deallocate(localData(ivar)%dvar2)
            allocate(localData(ivar)%dvar2(nDim1,nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 2D real space for localData data structure'; return; endif
          case('1D')
            if ( allocated(localData(ivar)%dvar1) ) deallocate(localData(ivar)%dvar1)
            allocate(localData(ivar)%dvar1(nPoly),stat=err)
            if(err/=0)then; message=trim(message)//'problem allocating 1D real space for localData data structure'; return; endif
        end select
    end select 
    do iPoly=1,nPoly
      iLocal=subPolyID(iPoly)
      select case(trim(data_meta(iVar)%vartype))
        case('integer')
          select case(trim(data_meta(iVar)%vardims))
            case('2D'); localData(ivar)%ivar2(:,iPoly) = entireData(ivar)%ivar2(:,iLocal)
            case('1D'); localData(ivar)%ivar1(iPoly)   = entireData(ivar)%ivar1(iLocal)
            end select
        case('double')
          select case(trim(data_meta(iVar)%vardims))
            case('2D'); localData(ivar)%dvar2(:,iPoly) = entireData(ivar)%dvar2(:,iLocal)
            case('1D'); localData(ivar)%dvar1(iPoly)   = entireData(ivar)%dvar1(iLocal)
          end select
      end select 
    end do
  end do 
  return 
end subroutine 

! private subroutine:
subroutine pop_hfrac( gammaUpdateMeta, hfrac, err, message)
  use var_lookup,           only: ixGamma 
  implicit none
  !input variables
  type(par_meta),       intent(in)  :: gammaUpdateMeta(:)
  !output variables
  real(dp),             intent(out) :: hfrac(:)
  integer(i4b),         intent(out) :: err         ! error code
  character(*),         intent(out) :: message     ! error message
  !local variables
  real(dp)                          :: dummy(20) 
  !logical(lgc)                      :: mask(20) 
  !integer(i4b)                      :: i

  err=0; message='pop_hfrac/'
  dummy=-999
  dummy(1)=gammaUpdateMeta(ixGamma%h1gamma1)%val
  dummy(2)=gammaUpdateMeta(ixGamma%h2gamma1)%val
  hfrac=dummy(1:nLyr-1)
  
  !check h parameters - now can chcek up to 5 layers
  !do i=1,size(calGammaMeta)
  !  if (gammaParMeta(i)%pname=="h1gamma1")then;dummy(1)=gammaParStr(i)%var(1);cycle;endif 
  !  if (gammaParMeta(i)%pname=="h2gamma1")then;dummy(2)=gammaParStr(i)%var(1);cycle;endif
  !  if (gammaParMeta(i)%pname=="h3gamma1")then;dummy(3)=gammaParStr(i)%var(1);cycle;endif
  !  if (gammaParMeta(i)%pname=="h4gamma1")then;dummy(4)=gammaParStr(i)%var(1);cycle;endif
  !enddo
  !mask=(dummy>0)
  !if ( count(mask)/=nLyr-1 ) stop 'number of h1gamma prameters mismatch with nLyr'
  !#hfrac=pack(dummy,mask)
  return
end subroutine 

! private function 
character(len=20) function int2str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (int2str, *) k
    int2str = adjustl(int2str)
end function

! private subroutine:
subroutine findix(scl, vec, iSelect, ierr, message)
  ! Find index where the value match up with scl  
  implicit none
  !input 
  integer(i4b),intent(in)              :: scl
  integer(i4b),intent(in)              :: vec(:)
  integer(i4b),intent(out)             :: iSelect
  integer(i4b), intent(out)            :: ierr      ! error code
  character(*), intent(out)            :: message   ! error message
  integer(i4b)                         :: i(1)

  ! initialize error control
  ierr=0; message='findix/' 

  i = minloc(abs(vec - scl))
  iSelect = i(1)  ! de-vectorize the found index 
  if(vec(iSelect) /= scl)&
    ierr=60; message=trim(message)//'unable to find matched value'; return  
  return
end subroutine

end module mpr_routine
