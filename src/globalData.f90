module globalData

 use nrtype
 use public_var
 use data_type,  only: var_meta, par_meta, cpar_meta, beta_meta, scale_meta, input_meta
 use var_lookup, only: nBeta, nGamma, nVarSoilData, nVarTopoData, nVarVegData, nVarMapData, nPrpVeg

implicit none

private

! master meta data including all available parameters
type(par_meta),       save,            public  :: gammaMeta(nGamma)         ! meta data for all the available gamma parameters
type(par_meta),       save,            public  :: betaMeta(nBeta)           ! meta data for all the available beta parameters
type(beta_meta),      save,            public  :: betaAncilMeta(nBeta)      ! meta data for parameter dependency for all the beta parameters 
! Based on a list of beta parameters to be estimated from "inParList" nml input
type(cpar_meta),      save,allocatable,public  :: calParMeta(:)             ! meta data for calibrating beta and gamma parameters based on 'inParList' nml input
type(cpar_meta),      save,allocatable,public  :: calGammaMeta(:)           ! subset of calParMeta including only gamma parameters. 
type(input_meta),     save,allocatable,public  :: inParMeta(:)              ! meta data for speficification of beta parameteters in "inParList" nml input
type(scale_meta),     save,allocatable,public  :: calScaleMeta(:)           ! meta data for scaling operator of beta paramets 
integer(i4b),         save,            public  :: calBetaOrderIdx(nBeta)    ! index of beta parameters sorted in computing order (only beta parameter estimated with MPR) 
character(len=strLen),save,allocatable,public  :: calBetaName(:)            ! name of beta parameters to be estimated with MPR based on 'inParList' nml input
character(len=strLen),save,allocatable,public  :: soilBetaCalName(:)        ! subset of "calBetaName"  including only soil parameters 
character(len=strLen),save,allocatable,public  :: vegBetaCalName(:)         ! subset of "calBetaName" including only vegetation parameters
real(dp),             save,allocatable,public  :: parArray(:,:)             ! calibrating parameter array - input for optimization routine
logical(lgc),         save,allocatable,public  :: parMask(:)                ! calibrating parameter mask vector - input for optimization routine
! Number of parameters 
integer(i4b),         save,            public  :: nCalBetaDir               ! number of beta parameters to be directly calibrated
integer(i4b),         save,            public  :: nCalGamma                 ! number of gamma parameters to be calibrated
integer(i4b),         save,            public  :: nCalPar                   ! sum of beta parameter directly calibrated and gamma parameters calibrated 
integer(i4b),         save,            public  :: nCalParSum                ! Total number of calibrating parameters = n(gamma)+m1(beta)+m2(beta_per_layer)*nLyr(m1+m2=nCalBetaDir)
integer(i4b),         save,            public  :: nBetaNeed                 ! number of beta parameters computed in MPR including calibrating beta and their dependent beta 
integer(i4b),         save,            public  :: nSoilBetaModel            ! number of soil parameters to be estimated with MPR
integer(i4b),         save,            public  :: nVegBetaModel             ! number of vegetation parameters to be estimated with MPR
integer(i4b),         save,            public  :: nSnowBetaModel            ! number of snow parameters to be estimated with MPR
! meta data for input data 
type(var_meta),       save,            public  :: map_meta  (nVarMapData)   ! mapping data
type(var_meta),       save,            public  :: sdata_meta(nVarSoilData)  ! soil data
type(var_meta),       save,            public  :: tdata_meta(nVarTopoData)  ! topographic data
type(var_meta),       save,            public  :: vdata_meta(nVarVegData)   ! vegetation data
type(var_meta),       save,            public  :: vprp_meta (nPrpVeg)       ! mapping data

end module globalData
