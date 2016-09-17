module mpr_routine 
  USE nrtype                                           ! Including numerical type definition 
  USE data_type                                        ! Including custum data structure definition
  USE public_var                                       ! Including common constant (physical constant, other e.g., missingVal, etc.)

  implicit none

  PRIVATE
   
  public:: mpr

contains

subroutine mpr(idModel,           &     ! model ID
               gammaPar,          &     ! array of gamma parameter 
               gammaParMeta,      &     ! array of gamma parameter metadata
               err, message)
  use model_wrapper,        only:read_soil_param
  use popMeta,              only:mprData
  use globalData,           only:parMaster, parSubset, gammaSubset, betaInGamma
  use globalData,           only:sdata_meta
  ! MPR core
!   USE soiltf                  only:comp_soil_model_param         ! Including Soil model parameter transfer function
!   USE vegtf,                  only:comp_veg_model_param          ! Including Veg model parameter transfer function
!   USE modelLayer,             only:comp_model_depth              ! Including model layr depth computation routine 
!   USE layerWeight,            only:map_slyr2mlyr                 ! Including model layr computation routine 
!   USE upscaling,              only:aggreg                        ! Including Upscaling operator 
 ! I/O
!   USE read_mapdata                                               ! routine to read mapping data into data structures
!   USE read_vegdata                                               ! routine to read veg data into data structures
  use read_soildata,          only:getData                       ! routine to read soil data into data structures
  use read_soildata,          only:mod_hslyrs                    ! routine to modify soil layer thickness and updata soil data structure
!   USE read_soildata,          only:get_topoinfo                  ! routine to extract topographic variables from soil data structure
!   USE read_ncdata,            only:get_vec_dvar                  !
  ! Named variable index 
  use var_lookup,             only: nPar
  USE var_lookup,             only:ixVarMapData,nVarMapData     ! index of map data variables and number of variables
  USE var_lookup,             only:ixVarVegData,nVarVegData     ! index of vege data variables and number of variables
  USE var_lookup,             only:ixVarSoilData,nVarSoilData   ! index of soil data variables and number of variables 
  USE var_lookup,             only:ixVarHru,nVarHru             ! index of model HRU variables and number of variables
!  USE var_lookup,             only:ixPrpSoil,nPrpSoil           ! index of soil properties and number of properties
!  USE var_lookup,             only:ixPrpVeg,nPrpVeg             ! index of veg properties and number of properties

   implicit none

  ! input
  integer(i4b),         intent(in)   :: idModel 
  real(dp),             intent(in)   :: gammaPar(:)              ! array of parameter value adjusted with calibration
  type(cpar_meta),      intent(in)   :: gammaParMeta(:)          ! array of calibrating meta data
  ! output
  integer(i4b),         intent(out)  :: err                      ! error code
  character(len=strLen),intent(out)  :: message                  ! error message 
  !! local
  character(len=strLen)              :: cmessage                 ! error message from downward subroutine
  ! index for printing (set to negative to supress printing
  integer,     parameter             :: iHruPrint = 1            ! model hru id for which everything is printed for checking
  integer(i4b),parameter             :: nSub=11                  ! max. number of Soil layer within Model layer
  integer(i4b)                       :: iSelect                  ! Array index selected
  integer(i4b)                       :: iGamma, iGammaMaster     ! index loop
  integer(i4b)                       :: iPoly                    ! Loop index of soil polygon
  integer(i4b)                       :: iSLyr                    ! Loop index of soil layer 
  integer(i4b)                       :: iMLyr                    ! loop index of model soil layer (1,...,n from top to bottom) 
  integer(i4b)                       :: iparm                    ! Loop index of model parameters (e.g., VIC)
  logical(lgt),allocatable           :: mask(:)                  ! mask for spolIdSub array 
  ! soil data stuff  
  real(dp), allocatable              :: param(:,:)               ! current soil parameter values for each HRU
  type(par_meta),allocatable         :: gammaParMasterMeta(:)
  integer(i4b)                       :: nSpoly                   ! number of soil polygon in entire soil data domain 
  integer(i4b)                       :: nSlyrs                   ! number of soil layers
  integer(i4b)                      :: iSclass                  ! ID (= index) of soil class 
  integer(i4b)                      :: iVarSoil                 ! Loop index of soil data 
  integer(i4b)                      :: iVarTopo                 ! Loop index of topographic properties 
  type(namevar)                     :: sdata(nVarSoilData)      ! soil data container for all the soil polygons
  type(namevar)                     :: sdataLocal(nVarSoilData) ! soil data container for local soil polygon 
  real(dp),        allocatable      :: mhlyr(:,:)               ! 
  integer(i4b),    allocatable      :: soilClass(:)             ! soil class array (e.g., USDA soil class used for STATSGO)
  type(var_d),     allocatable      :: scls2prp(:)              ! storage of property value for each soil class
  type(namedvar2), allocatable      :: sprp(:)                  ! storage of soil property values for 2D(poly, lyr) 
  type(namedvar) , allocatable      :: topo(:)                  ! storage of topo property values in subset domain-1D(poly)
  type(namedvar2), allocatable      :: sprpLocal(:)             ! storage of soil property values for 2D(poly, lyr) 
  type(namedvar) , allocatable      :: topoLocal(:)             ! storage of topo property values in subset domain-1D(poly)
  integer(i4b),    allocatable      :: sclsLocal(:,:)           ! soil class id for each poly and depth for subset of soil data
  ! vege data stuff
  integer(i4b)                      :: nVpoly                   ! number of vege polygon (grid box) in entire vege data domain 
  integer(i4b),    allocatable      :: vegClass(:)              ! veg class array (e.g., IGBP)
  integer(i4b),    allocatable      :: vclsLocal(:)             ! veg class id for each poly for subset of veg data
  type(var_d),     allocatable      :: vcls2prp(:)              ! storage of property value for each veg class
  type(namedvar),  allocatable      :: vprp(:)                  ! storage of veg property values for 1D(poly) 
  type(namedvar),  allocatable      :: vprpLocal(:)             ! storage of veg property values for 1D(poly) 
  type(namevar)                     :: vdata(nVarVegdata)       ! veg data container for all the veg polygons
  integer(i4b)                      :: iVclass                  ! ID (= index) of vege class 
  integer(i4b)                      :: iPrpVeg                  ! Loop index of veg properties 
  ! Model stuff 
  integer(i4b)                      :: nVegParModel             ! Number of model vege parameters
  integer(i4b)                      :: nSoilParModel            ! Number of model soil parameters
  type(dat_d2d),   allocatable      :: ParSxySz(:)              ! storage of model soil parameter for 2D field (2nd variable)-soil poy x soil layer
  type(dat_d2d),   allocatable      :: ParSxyMz(:)              ! storage of model soil parameter for 2D field (2nd variable)-soil poy x model layer
  type(namedvar),  allocatable      :: vegParMxy(:)             ! storage of model soil parameter for 2D field (2nd variable)-model hru x model layer
  type(namedvar),  allocatable      :: ParVxy(:)                ! storage of model vege parameter for 1D or 2D field - vege poly (x month)
  ! mapping stuff  
  integer(i4b),    allocatable      :: polySub_org(:)           ! list of ID (=index) of soil polygons contributing model hru including missing value
  integer(i4b),    allocatable      :: polySub(:)               ! list of ID (=index) of soil polygons contributing model hru
  integer(i4b),    allocatable      :: vPolySub_org(:)          ! list of ID (=index) of veg polygons contributing model hru including missing value
  integer(i4b),    allocatable      :: vPolySub(:)              ! list of ID (=index) of veg polygons contributing model hru 
  logical(lgt),    allocatable      :: vmask(:)                 ! mask for vpolIdSub array 
  integer(i4b)                      :: nShru                    ! number of hrus in soil mapping file) 
  integer(i4b)                      :: nVhru                    ! number of hrus in vege mapping file (should be equal to nVhru)
  integer(i4b)                      :: nOverSpoly               ! number of overlapped soil polygon for each model HRU 
  integer(i4b)                      :: nOverVpoly               ! number of overlapped vege polygon for each model HRU 
  integer(i4b)                      :: nSpolyLocal              ! number of subset overlapped soil polygon for each model HRU 
  integer(i4b)                      :: nVpolyLocal              ! number of subset overlapped vege polygon for each model HRU 
  type(mapvar)                      :: mapdata(2)               ! map data container for all the soil polygons
  integer(i4b)                      :: iVarHru                  ! Loop index of Grid properties
  integer(i4b)                      :: iHru                     ! loop index of hrus 
  integer(i4b)                      :: iSub                     ! Loop index of multiple soi layers in model layer
  real(dp),        allocatable      :: swgtSub_org(:)           ! Areal weight of soil polygon for all model hrus
  real(dp),        allocatable      :: swgtSub(:)               ! adjusted Areal weight of soil polygon for all model hrus
  real(dp),        allocatable      :: vwgtSub_org(:)           ! Areal weight of veg polygon for one model hrus
  real(dp),        allocatable      :: vwgtSub(:)               ! adjusted Areal weight of veg polygon for one model hrus
  type(poly),      allocatable      :: soil2model_map(:)        ! data structure to hold weight and index of contributing soil layer per model layer and per soil polygon
  type(dat_d1d),   allocatable      :: hruPrp(:)                !
  type(lyr_d),     allocatable      :: paramvec(:)              !
  
  ! initialize error control
  err=0; message='mpr/'
  !(0) Preparation
  ! (0.1) List all the gamma parameter from master pararameter metadata
  allocate(mask(size(parMaster)))
  mask=(parMaster(:)%beta/="beta")
  allocate(gammaParMasterMeta(count(mask)))
  gammaParMasterMeta = pack(parMaster,mask)
  ! (0.2) separate beta parameter into veg and soil 
  ! (0.3) get number of soil and vege parameters 
  nSoilParModel=size(gammaPar)
  ! (0.4) swap gammaParMasterMeta%val with gammaPar value
  do iGamma=1,size(gammaPar)
    do iGammaMaster=1,size(gammaParMasterMeta)
      if ( gammaParMasterMeta(iGammaMaster)%pname==gammaParMeta(iGamma)%pname )then
        gammaParMasterMeta(iGammaMaster)%val=gammaPar(iGamma)
      endif
    enddo
  enddo
  ! (0.5) read parameter file 
  allocate(param(nHru,TotNparVic))
  call read_soil_param(idModel, param, err, message) 
  if(err/=0)then; message=message//cmessage; return; endif

  ! (1) Get Geophysical data 
  ! *****
  ! (1.1) soil data  
  ! *********************************************
  ! (1.1.1) Read in Soil data netCDF...
!  allocate(soilClass(nSclass))
!  call getData(trim(mpr_input_dir)//trim(fname_soil),  & ! input: soil data (netCDF)
!                    sdata_meta,                    & ! input: soil data meta
!                    dname_spoly,                   & ! input: spatial dimension (polygon ID)
!                    dname_slyrs,                   & ! input: spatial dimension (polygon ID)
!                    sdata,                         & ! input-output: soil data structure
!                    nSpoly,                        & ! output: number of dimension (i.e. number of soil polygon)
!                    nSlyrs,                        & ! output: number of dimension (i.e. number of soil layer)
!                    err, cmessage)
!  if(err/=0)then; message=message//cmessage; return; endif 
  ! (1.1.2) modify soil layer thickness in sdata data structure
!  call mod_hslyrs(sdata,err,cmessage)
!  if(err/=0)then; message=message//cmessage; return; endif 
!  print*, sdata(ixVarSoilData%hslyrs)%dvar2(:,10)
!  ! (1.1.2a) Read in soil class lookup table
!   call read_scls(trim(mpr_input_dir)//trim(sclass_table), &
!                  nSclass,                              &
!                  soilClass,                            &
!                  scls2prp,                             &
!                  ierr, cmessage)
!   if(ierr/=0)then; message=message//cmessage; return; endif 
!  ! (1.1.3) Populate soil properties for each layer and polygon 
!  ! Get soil class and associated soil properties for only soil polygons that overlap current model hru
!   allocate(topo(nVarTopo),stat=ierr);   if(ierr/=0) call handle_err(ierr,'error allocating for topo')
!   do iVarTopo=1,nVarTopo
!      allocate(topo(iVarTopo)%varData(nSpoly),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for topo%varData')
!    enddo
!    call get_topoinfo(sdata, topo, ierr, cmessage)
!    if(ierr/=0)then; message=message//cmessage; return; endif
!  ! Use USDA soil class - soil properties look up table
!  ! 1.1 read soil class-soil property look up table 
!  ! 1.2 read USDA soil class for each polygon and layer 
!    allocate(sprp(nPrpSoil),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for sprp')
!    do iPrpSoil=1,nPrpSoil 
!      allocate(sprp(iPrpSoil)%varData(nSlyrs,nSpoly))
!    enddo
!    call map_scls2prp(sdata(ixVarSoilData%soilclass)%ivar2, scls2prp, sdata, sprp) 
!  
!  ! *****
!  ! (1.2) Read in veg class - properties lookup table ...
!  ! *********************************************
!  ! (1.2.1) Read in veg data netCDF...
!   allocate(vegClass(nVclass),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vegClass')
!   allocate(vcls2prp(nVclass),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vcls2prp')
!   do iVclass=1,nVclass
!     allocate(vcls2prp(iVclass)%var(nPrpVeg), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vcls2prp%var') 
!   enddo
!   call getVegData(trim(mpr_input_dir)//trim(fname_veg), & ! input: file name
!                   vdata_meta,                       & ! input: veg data meta
!                   dname_vpoly,                      & ! input: dimension name for veg polygon
!                   vdata,                            & ! input-output: veg data structure
!                   nVpoly,                           & ! output: number of veg polygon
!                   ierr, cmessage)                     ! output: error control
!   if(ierr/=0)then; message=message//cmessage; return; endif
!  ! (1.2.2) Read in veg class-property lookup table 
!   call getvegClassLookup(trim(mpr_input_dir)//trim(vclass_table), &
!                          nVclass,                             &
!                          vegClass,                            &
!                          vcls2prp,                            &
!                          ierr, cmessage)
!   if(ierr/=0)then; message=message//cmessage; return; endif
!  
!  ! (1.2.3) populate veg properties for each polygon 
!   allocate(vprp(nPrpVeg),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for sprp')
!   do iPrpVeg=1,nPrpVeg
!     allocate(vprp(iPrpVeg)%varData(nVpoly), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for sprp%varData')
!   enddo
!   call map_vcls2prp(vdata, vcls2prp, vegClass, vprp, ierr, cmessage)
!   if(ierr/=0)then; message=message//cmessage; return; endif
!  
!  ! *****
!  ! (2.) Read in mapping netcdf 
!  ! *********************************************
!   ! (2.1) mapping soil polygon to model hru 
!    call getMapData(trim(mpr_input_dir)//trim(fname_smapping), &   ! input: file name
!                    'soil',                                &   ! input: geophysical data type 
!                    mapdata_meta,                          &   ! input: map data meta
!                    dname_hru,                             &   ! input: dimension name for hru 
!                    dname_overSpoly,                       &   ! input: dimension name for overlap polygon 
!                    mapdata,                               &   ! input-output: map data structure
!                    nShru,                                 &   ! output: number of hru 
!                    nOverSpoly,                            &   ! output: max number of overlap polygons
!                    ierr,cmessage)                             ! output: error control
!    if (err/=0)then; message=message//cmessage; return; endif
!    ! (2.2) mapping vege polygon to model hru 
!    call getMapData(trim(mpr_input_dir)//trim(fname_vmapping), &   ! input: file name
!                    'veg',                                 &   ! input: geophysical data type 
!                    mapdata_meta,                          &   ! input: map data meta
!                    dname_hru,                             &   ! input: dimension name for hru 
!                    dname_overVpoly,                       &   ! input: dimension name for overlap Polygon 
!                    mapdata,                               &   ! input-output: map data structure
!                    nVhru,                                 &   ! output: number of hru 
!                    nOverVpoly,                            &   ! output: max number of overlap polygon
!                    ierr,cmessage)                             ! output: error control
!    if (err/=0)then; message=message//cmessage; return; endif
!  
!    ! (2.3) Check if "hru" variables in vege and soil mappling netcdf are identifal
!    if ( nShru /= nVhru ) then
!      call handle_err(10,'Different number of hru in vege and soil mapping file')  
!    end if
!    if ( minval(abs(mapdata(1)%var(ixVarMapData%hru_id)%ivar1-mapdata(2)%var(ixVarMapData%hru_id)%ivar1)) /= 0 ) then
!      call handle_err(11,'different hru id in vege and soil mapping file')
!    end if
!  
!    ! (2.4) assign pointers for saving type 
!    associate( nHru => nShru, &
!               hruID => mapdata(1)%var(ixVarMapData%hru_id)%ivar1, 
!               swgt  => mapdata(1)%var(ixVarMapData%weight)%dvar2
!               vwgt  => mapdata(2)%var(ixVarMapData%weight)%dvar2
!               overSpolyID  => mapdata(1)%var(ixVarMapData%overlapPolyId)%ivar2
!               overVpolyID  => mapdata(2)%var(ixVarMapData%overlapPolyId)%ivar2)
!  
!  ! *****
!  ! (3.) Read in hru property data 
!  ! *********************************************
!    allocate(hruPrp(nVarHru),stat=ierr); if(ierr/=0) message=message//cmessage; return
!    do iVarHru=1,nVarHru
!      allocate(hruPrp(iVarHru)%dat(nHru),stat=ierr); if(ierr/=0)  message=message//cmessage; return
!    enddo
!    do iVarHru=1,nVarHru
!      call get_vec_dvar(trim(mpr_input_dir)//trim(fname_smapping),hru_meta(iVarHru)%varname,hruPrp(iVarHru)%dat,1,nHru,ierr,cmessage); call handle_err(ierr,cmessage)
!    enddo
!  
!  ! *****
!  ! (3.) Computing model parameters from geophysical properties (soil, vege) and scale them to model HRU 
!  ! *********************************************************************
!    allocate(ParMxyMz(nSoilParModel),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParMxyMz')
!    do iparm=1,nSoilParModel
!      allocate(ParMxyMz(iparm)%dat(nMLyr,nHru),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParMxyMz%dat')
!    enddo
!    allocate(vegParMxy(nVegParModel),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vegParMxy')
!    do iparm=1,nVegParModel
!      allocate(vegParMxy(iparm)%varData(nHru),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vegParMxy%varData')
!    enddo
!    !!! ---------------------------------------------
!    !!! Start of model hru loop (from mapping file) !!!
!    !!! ---------------------------------------------
!    do iHru=1,nHru
!    ! *****
!    ! (3.1) Extract soil poly ID, weight polygon , and soil properties for current model hru 
!    ! *********************************************************************
!      ! Select list of soil polygon ID contributing a current hru
!      allocate(polySub_org(nOverSpoly),stat=ierr); if(ierr/=0) message=trim(message)//'error allocating for polySub_org'; return
!      allocate(swgtSub_org(nOverSpoly),stat=ierr); if(ierr/=0) message=trim(message)//'error allocating for swgtSub_org'; return
!      allocate(mask(nOverSpoly),stat=ierr);        if(ierr/=0) message=trim(message)//'error allocating for mask'; return
!      polySub_org = OverSpolyID(:,iHru);
!      swgtSub_org = swgt(:,iHru);
!      mask        = ( polySub_org /= imiss )
!  
!      allocate(polySub(count(mask)),stat=ierr); if(ierr/=0) message=trim(message)//'error allocating for polySub'; return
!      allocate(swgtSub(count(mask)),stat=ierr); if(ierr/=0) message=trim(message)//'error allocating for swgtSub'; return
!      polySub   = pack(polySub_org,mask)
!      swgtSub   = pack(swgtSub_org,mask)
!      nSpolyLocal = size(polySub)
!  
!      allocate(sprpLocal(nPrpSoil),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for sprpLocal')
!      do iPrpSoil=1,nPrpSoil 
!        allocate(sprpLocal(iPrpSoil)%varData(nSlyrs,nSpolyLocal))
!      enddo
!      allocate(topoLocal(nVarTopo),stat=ierr);   if(ierr/=0) call handle_err(ierr,'error allocating for topo')
!      do iVarTopo=1,nVarTopo
!        allocate(topoLocal(iVarTopo)%varData(nSpolyLocal),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for topo%varData')
!      enddo
!  
!      do iPoly = 1,nSpolyLocal
!        do iPrpSoil = 1,nPrpSoil
!          sprpLocal(iPrpSoil)%varData(:,iPoly)=sprp(iPrpSoil)%varData(:,polySub(iPoly)) 
!        enddo
!        do iVarTopo=1,nVarTopo
!          topoLocal(iVarTopo)%varData(iPoly)=topo(iVarTopo)%varData(polySub(iPoly))
!        enddo
!      end do 
!      call subSoilData(sdata, nSpolyLocal, sdataLocal, ierr, message)
!  
!      !Check
!      if ( iHru == iHruPrint ) then
!        print*,' '
!        print*,'****************************************************'
!        print*,'HRU, hruID = ',iHru,hruID(iHru)
!        print*,'****************************************************'
!        print*,'(1.1) Print list of soil polygon ID and weigth'
!        write(*,"(' polyID = ',100I9)") (polySub(iPoly), iPoly=1,nSpolyLocal)
!        write(*,"(' weight = ',100f9.3)") (swgtSub(iPoly), iPoly=1,nSpolyLocal)
!        print*,'(1.2) Print soil property for polygon and layer'
!        do iVarTopo=1,nVarTopo
!          write(*,"(1X,A20,'= ',100f9.3)") stopo_meta(iVarTopo)%varname, (topoLocal(iVarTopo)%varData(iPoly), iPoly=1,nSpolyLocal)
!        enddo
!        do iSLyr = 1,nSlyrs
!          write(*,"(' Layer           = ',I3)")(iSLyr)
!          do iPrpSoil=1,nPrpSoil
!            write(*,"(1X,A17,'= ',100f9.3)") sprp_meta(iPrpSoil)%varname, (sprpLocal(iPrpSoil)%varData(iSLyr,iPoly), iPoly=1,nSpolyLocal)
!          enddo
!        enddo
!      endif
!    ! *****
!    ! (3.2) Extract veg poly ID, weight polygon, and veg properties for current model hru 
!    ! *********************************************************************
!      ! Select list of veg polygon ID contributing a current hru 
!      allocate(vPolySub_org(nOverVpoly),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vPolySub_org')
!      allocate(vwgtSub_org(nOverVpoly),stat=ierr);  if(ierr/=0) call handle_err(ierr,'error allocating for vwgtSub_org')
!      allocate(vmask(nOverVpoly),stat=ierr);        if(ierr/=0) call handle_err(ierr,'error allocating for mask')
!      
!      vPolySub_org = overVpolyID(:,iHru);
!      vwgtSub_org  = vwgt(:,iHru);
!      vmask        = ( vPolySub_org /= imiss )
!  
!      allocate(vPolySub(count(vmask)),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vPolySub')
!      allocate(vwgtSub(count(vmask)),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for swgtSub')
!      vPolySub    = pack(vPolySub_org,vmask)
!      vwgtSub     = pack(vwgtSub_org,vmask)
!      nVpolyLocal = size(vPolySub)
!  
!      allocate(vprpLocal(nPrpVeg),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for vprp')
!      do iPrpVeg=1,nPrpVeg 
!        allocate(vprpLocal(iPrpVeg)%varData(nVpolyLocal))
!      enddo
!      allocate(vclsLocal(nVpolyLocal))
!  
!      do iPoly = 1,nVpolyLocal
!        do iPrpVeg = 1,nPrpVeg
!          vprpLocal(iPrpVeg)%varData(iPoly) =vprp(iPrpVeg)%varData(vPolySub(iPoly)) 
!        enddo
!        vclsLocal(iPoly) = vdata(ixVarVegData%vegclass)%ivar1(vPolySub(iPoly))
!      end do
!      !Check
!      if ( iHru == iHruPrint ) then
!        print*,' '
!        print*,'(1.3) Print list of vege polygon ID and weigth'
!        write(*,"(' polyID = ',100I7)") (vPolySub(iPoly), iPoly=1,nVpolyLocal)
!        write(*,"(' weight = ',100f7.3)") (vwgtSub(iPoly), iPoly=1,nVpolyLocal)
!        write(*,"(' vegclass = ',100I5)") (vclsLocal(iPoly), iPoly=1,nVpolyLocal)
!        print*,'(1.3) Print veg property for polygon'
!        do iPrpVeg=1,nPrpVeg
!          write(*,"(1X,A12,'= ',100f7.3)") vprp_meta(iPrpVeg)%varname, (vprpLocal(iPrpVeg)%varData(iPoly), iPoly=1,nVpolyLocal)
!        enddo
!      endif
!    ! *****
!    ! (3.4) Compute model soil parameters using transfer function
!    ! *********************************************************
!      ! compute model soil parameters
!      allocate(ParSxySz(nSoilParModel),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParSxySz')
!      do iparm=1,nSoilParModel
!        allocate(ParSxySz(iparm)%dat(nSlyrs,nSpolyLocal),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ')
!      enddo
!      call comp_soil_model_param(ParSxySz, nSlyrs, nSpolyLocal, sdataLocal, ierr, cmessage)
!      call comp_pedotf(sdata, sprp, ierr, cmessage);  if(ierr/=0) call handle_err(ierr, cmessage)
!      !Check
!      if ( iHru == iHruPrint ) then
!        print*,'(2) Print Model parameter for polygon and layer'
!        do iSLyr = 1,nSlyrs
!          write(*,"(' Layer        = ',I3)")     (iSLyr)
!          do iparm=1,nSoilParModel
!            write(*,"(1X,A17,'= ',100f9.3)") spar_meta(iparm)%parname, (ParSxySz(iparm)%dat(iSLyr,iPoly), iPoly=1,nSpolyLocal)
!          enddo
!        enddo
!      endif
!    ! *****
!    ! (3.5) Compute Model vege parameters using transfer function
!    ! *********************************************************
!      ! compute model veg parameters
!      allocate(ParVxy(nVegParModel),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParVxy')
!      do iparm=1,nVegParModel
!        select case ( vpar_meta(iparm)%dims )
!          case( '1D' , 'ST' )
!            allocate(ParVxy(iparm)%varData(nVpolyLocal),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParVxy%varData')
!        end select
!      enddo
!      call comp_veg_model_param(ParVxy, vprpLocal, vpar_meta, model_name)
!      !Check
!      if ( iHru == iHruPrint ) then
!        print*,'(2) Print Model parameter for Vege polygon'
!        do iparm=1,nVegParModel
!          write(*,"(1X,A10,'= ',100f5.2)") vpar_meta(iparm)%parname, (ParVxy(iparm)%varData(iPoly), iPoly=1,nVpolyLocal)
!        enddo
!      endif
!  ! ***********
!  ! (4) Spatial aggregation of Model parameter 
!  ! *********************************************************************
!    ! **********
!    ! (4.1) Aggregate model parameveter vertical direction - soil data layers to model soil layers
!    ! *********************************************************************
!    ! (4.1.1) Compute Model layer depth 
!      allocate(ParSxyMz(nSoilParModel),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParSxyMz')
!      do iparm=1,nSoilParModel
!        allocate(ParSxyMz(iparm)%dat(nMLyr,nSpolyLocal),stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for ParSxyMz%dat')
!      enddo
!      call comp_model_depth(ParSxyMz, nMLyr, nSpolyLocal, sprpLocal, sclsLocal, model_name, ierr, cmessage) 
!      if(ierr/=0) call handle_err(ierr,cmessage)
!  
!      !Check
!      if ( iHru == iHruPrint ) then
!        print*, '(3.1.1) Print model depth ---' 
!        do iMLyr = 1,nMLyr
!        print*, 'Layer = ',iMLyr 
!          write(*,"('z= ', 100f9.3)") (ParSxyMz(ixModelDepth)%dat(iMLyr,iPoly), iPoly=1,nSpolyLocal)
!        enddo 
!      endif
!      
!    ! (4.1.2) Computing vertical weight of soil layer for each model layer
!      allocate(soil2model_map(nSpolyLocal), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for soil2model_map')
!      do iPoly=1,nSpolyLocal
!        allocate(soil2model_map(iPoly)%layer(nMLyr), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for soil2model_map%layer')
!        do iMLyr=1,nMLyr
!          allocate(soil2model_map(iPoly)%layer(iMLyr)%weight(nSub), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for soil2model_map%layer%weight')
!          soil2model_map(iPoly)%layer(iMLyr)%weight = dmiss
!          allocate(soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(nSub), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for soil2model_map%layer%ixSubLyr')
!          soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr = imiss
!        enddo
!        
!        second: associate( zSoil =>sprpLocal(ixPrpSoil%z)%varData(:,iPoly)
!                           hSoil =>sprpLocal(ixPrpSoil%h)%varData(:,iPoly)
!                           zModel=>ParSxyMz(ixModelDepth)%dat(:,iPoly) )
!        call map_slyr2mlyr(zSoil, hSoil, zModel, soil2model_map(iPoly)%layer, ierr, cmessage)
!        call handle_err(ierr,cmessage)
!        end associate second
!      enddo
!  
!      !Check
!      if ( iHru == iHruPrint ) then
!        print*, '(3.1.1) Print info on mapping from soil to model layer --' 
!        do iPoly = 1,nSpolyLocal
!          do iMLyr = 1,nMLyr
!            print*, 'Poly, Layer=',iPoly,iMLyr
!            write(*,"('weightSubLyr= ', 11f9.3)") (soil2model_map(iPoly)%layer(iMLyr)%weight(iSub), iSub=1,nSub)
!            write(*,"('indexSubLyr = ', 11I7)") (soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(iSub), iSub=1,nSub)
!          enddo 
!        enddo
!      endif
!  
!    ! (4.1.3) Computing vertical weight for each model layer
!      do iPoly=1,nSpolyLocal
!        do iMLyr = 1,nMLyr
!          !!-- Start of Memory allocation and initialization
!          allocate( paramvec(nSoilParModel), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for paramvec(nSub)') 
!          do iparm = 1,nSoilParModel
!            allocate( paramvec(iparm)%layer(nSub), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for paramvec(nSub)') 
!            paramvec(iparm)%layer = 0._dp
!          enddo
!          !!-- End of Memory allocation and initialization
!          do iSub=1,nSub
!            if (soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(iSub) > 0 ) then 
!              ! Get index of the Layer ID from STATSGO that match up with the ID specified  
!              iSelect = soil2model_map(iPoly)%layer(iMLyr)%ixSubLyr(iSub)
!              do iparm = 1,nSoilParModel
!                paramvec(iparm)%layer(iSub) = ParSxySz(iparm)%dat(iSelect,iPoly)
!              enddo ! End of iparm loop
!            endif
!          enddo ! End of iSub loop 
!          ! Aggregation per iPoly and iMLayer
!          do iparm = 1,nSoilParModel
!            if (spar_meta(iparm)%v_agg) then
!              call aggreg(ParSxyMz(iparm)%dat(iMLyr,iPoly),             &
!                          soil2model_map(iPoly)%layer(iMLyr)%weight(:), &
!                          paramvec(iparm)%layer(:),                     &
!                          dmiss,                                        &
!                          spar_meta(iparm)%vaggmethod)
!            endif
!          enddo
!          deallocate(paramvec, stat=ierr); if(ierr/=0) call handle_err(ierr,'error deallocating space for paramvec')
!        enddo ! end of iMLyr (model layer) loop
!      enddo ! end of iPoly
!  
!    ! ************
!    ! (4.2) Aggregate model parameveter horizontally - use spatial weight netCDF 
!    ! *********************************************************************
!      do iMLyr=1,nMLyr
!        ! Aggregation per grid box and all MLayers                    
!        allocate( paramvec(nSoilParModel), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for paramvec(nSub)') 
!        do iparm = 1,nSoilParModel
!          allocate( paramvec(iparm)%layer(nSpolyLocal), stat=ierr); if(ierr/=0) call handle_err(ierr,'error allocating for paramvec(nSub)') 
!          paramvec(iparm)%layer = 0._dp
!        enddo
!        do iPoly=1,nSpolyLocal
!          do iparm = 1,nSoilParModel
!              paramvec(iparm)%layer(iPoly) = ParSxyMz(iparm)%dat(iMLyr,iPoly)
!          enddo
!        enddo
!        do iparm = 1,nSoilParModel
!          if (spar_meta(iparm)%h_agg) then
!            call aggreg(ParMxyMz(iparm)%dat(iMLyr,iHru), &
!                        swgtsub(:),                      &
!                        paramvec(iparm)%layer(:),        &
!                        dmiss,                           &
!                        spar_meta(iparm)%haggmethod)
!          if ( iHru == iHruPrint ) then
!            print*,'-----------------------------------'
!            print*,'Aggregated soil parameter '
!            write(*,"(1X,A17,'= ',100f9.3)") spar_meta(iparm)%parname, (ParMxyMz(iparm)%dat(iMLyr,iHru))
!          endif
!          endif
!        enddo
!        ! Deallocate memory for next grid box iteration
!        deallocate(paramvec, stat=ierr); if(ierr/=0) call handle_err(ierr,'error deallocating for paramvec')
!      enddo
!     
!      ! Aggregation veg parameter per model hru 
!      do iparm = 1,nVegParModel
!        if (vpar_meta(iparm)%h_agg) then
!          call aggreg(vegParMxy(iparm)%varData(iHru), &
!                      vwgtSub(:),                     &
!                      ParVxy(iparm)%varData,          &
!                      dmiss,                          &
!                      vpar_meta(iparm)%haggmethod)
!          if ( iHru == iHruPrint ) then
!            print*,'-----------------------------------'
!            print*,'Aggregated vege parameter '
!            write(*,"(1X,A17,'= ',100f9.3)") vpar_meta(iparm)%parname, (vegParMxy(iparm)%varData(iHru))
!          endif
!        endif
!      enddo
!  
!      ! deallocate memory 
!      deallocate(polySub_org,stat=ierr);      if(ierr/=0) call handle_err(ierr,'error deallocating for polyIdSub_org')
!      deallocate(swgtSub_org,stat=ierr);      if(ierr/=0) call handle_err(ierr,'error deallocating for swgtSub_org')
!      deallocate(mask,stat=ierr);             if(ierr/=0) call handle_err(ierr,'error deallocating for mask')
!      deallocate(polySub,stat=ierr);          if(ierr/=0) call handle_err(ierr,'error deallocating for polyIdSub')
!      deallocate(swgtSub,stat=ierr);          if(ierr/=0) call handle_err(ierr,'error deallocating for swgtSub')
!      deallocate(vpolySub_org,stat=ierr);     if(ierr/=0) call handle_err(ierr,'error deallocating for vpolyIdSub_org')
!      deallocate(vwgtSub_org,stat=ierr);      if(ierr/=0) call handle_err(ierr,'error deallocating for vwgtSub_org')
!      deallocate(vmask,stat=ierr);            if(ierr/=0) call handle_err(ierr,'error deallocating for vmask')
!      deallocate(vpolySub,stat=ierr);         if(ierr/=0) call handle_err(ierr,'error deallocating for vpolyIdSub')
!      deallocate(vwgtSub,stat=ierr);          if(ierr/=0) call handle_err(ierr,'error deallocating for vwgtSub')
!      deallocate(sclsLocal,stat=ierr);        if(ierr/=0) call handle_err(ierr,'error deallocating for sclsLocal')
!      deallocate(vclsLocal,stat=ierr);        if(ierr/=0) call handle_err(ierr,'error deallocating for vclsLocal')
!      deallocate(sprpLocal,stat=ierr);        if(ierr/=0) call handle_err(ierr,'error deallocating for sprpLocal')
!      deallocate(topoLocal,stat=ierr);        if(ierr/=0) call handle_err(ierr,'error deallocating for topoLocal')
!      deallocate(vprpLocal,stat=ierr);        if(ierr/=0) call handle_err(ierr,'error deallocating for sprpLocal')
!      deallocate(soil2model_map, stat=ierr);  if(ierr/=0) call handle_err(ierr,'error deallocating for soil2model_map')
!  
!    enddo 
    !!! End of model hru loop !!!
!    end associate
!    print*, hruPrp, ParMxyMz
    return
end subroutine mpr

!! routine to output parameter in netCDF
!subroutine output_nc( ncname, ierr, cmessage)
!  use write_param_nc,only:defSoilNetCDF                 ! routine to define soil parameter variables in netCDF output
!  use write_param_nc,only:defVegNetCDF                  ! routine to define soil parameter variables in netCDF output
!  use write_param_nc,only:write_vec_ivar                ! routine to write 1D integer variable
!  use write_param_nc,only:write_vec_dvar                ! routine to write 1D real variable
!  use write_param_nc,only:write_array2_ivar             ! routine to write 2D integer variable
!  use write_param_nc,only:write_array2_dvar             ! routine to write 2D real variable 
!  implicit none
!
!  ! Create NetCDF file and define variables to be written
!  call defSoilNetCDF(trim(output_dir)//trim(soilnc_output), & ! input: output file name
!                     nHru,                                  & ! input: number of Hru 
!                     nMLyr,                                 & ! input: number of model layers
!                     spar_meta,                             & ! input: model parameter metadata 
!                     hru_meta,                              & ! input: model grid metadata 
!                     ierr, cmessage)                          ! output: error control
!  call handle_err(ierr, cmessage)
!  ! write hru dimension, (input = filename, variable name, variable vector, start index; output = error control)
!  call write_vec_ivar(trim(output_dir)//trim(soilnc_output),dname_mhru,hruID,1,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write layer dimension 
!  call write_vec_ivar(trim(output_dir)//trim(soilnc_output),dname_mlyr,(/(iMLyr,iMLyr=1,nMLyr)/),1,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write longitude 
!  call write_vec_dvar(trim(output_dir)//trim(soilnc_output),hru_meta(ixVarHru%lon)%varname,hruPrp(ixVarHru%lon)%dat,1,nHru,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write latitude 
!  call write_vec_dvar(trim(output_dir)//trim(soilnc_output),hru_meta(ixVarHru%lat)%varname,hruPrp(ixVarHru%lat)%dat,1,nHru,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write VIC expt parameter 
!  do iparm =1,nSoilParModel
!    if (spar_meta(iparm)%varout) then
!      select case ( spar_meta(iparm)%dims )
!        case('2D')
!          call write_vec_dvar(trim(output_dir)//trim(soilnc_output),spar_meta(iparm)%parname,ParMxyMz(iparm)%dat(1,:),1,nHru,ierr,cmessage)
!          call handle_err(ierr,cmessage)
!        case('3D')
!          call write_array2_dvar(trim(output_dir)//trim(soilnc_output),spar_meta(iparm)%parname,ParMxyMz(iparm)%dat,(/1,1/),(/nMLyr,nHru/),ierr,cmessage)
!          call handle_err(ierr,cmessage)
!      end select
!    endif     
!  enddo
!  
!  ! Create NetCDF file and define variables to be written
!  call defVegNetCDF(trim(output_dir)//trim(vegenc_output),  & ! input: output file name
!                     nHru,                                  & ! input: number of Hru 
!                     vpar_meta,                             & ! input: model parameter metadata 
!                     hru_meta,                              & ! input: model grid metadata 
!                     ierr, cmessage)                          ! output: error control
!  call handle_err(ierr, cmessage)
!  ! write hru dimension, (input = filename, variable name, variable vector, start index; output = error control)
!  call write_vec_ivar(trim(output_dir)//trim(vegenc_output),dname_mhru,hruID,1,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write longitude 
!  call write_vec_dvar(trim(output_dir)//trim(vegenc_output),hru_meta(ixVarHru%lon)%varname,hruPrp(ixVarHru%lon)%dat,1,nHru,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write latitude 
!  call write_vec_dvar(trim(output_dir)//trim(vegenc_output),hru_meta(ixVarHru%lat)%varname,hruPrp(ixVarHru%lat)%dat,1,nHru,ierr,cmessage); call handle_err(ierr,cmessage)
!  ! write VIC expt parameter 
!  do iparm =1,nVegParModel
!    if (vpar_meta(iparm)%varout) then
!      select case ( vpar_meta(iparm)%dims )
!        case('1D')
!            call write_vec_dvar(trim(output_dir)//trim(vegenc_output),vpar_meta(iparm)%parname,vegParMxy(iparm)%varData,1,nHru,ierr,cmessage)
!            call handle_err(ierr,cmessage)
!        case default; call handle_err(55,'variable dimension not supported')
!      end select
!    endif     
!  enddo
!  return
!end subroutine output_nc

end module mpr_routine
