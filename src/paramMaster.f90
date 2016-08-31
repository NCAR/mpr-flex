module paramMaster

implicit none

private

public::popMeta

contains

  subroutine popMeta(err,message)
  
    use nrtype
    use data_type,  only:par_meta
    use var_lookup, only:ixGamma
    use var_lookup, only:ixBeta
    use globalData, only:gammaMaster
    use globalData, only:betaMaster
    
    implicit none
  
    !output variable
    integer(i4b),intent(out)      :: err     ! error code
    character(*),intent(out)      :: message ! error message
    
    ! initialize error control
    err=0; message='popMeta/'
    
    ! -----
    !  Master list of gamma parameters  
    ! -----------------------
    !ks transfer function
    gammaMaster(ixGamma%ks1gamma1)       = par_meta('ks1gamma1'      ,   -0.6_dp ,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%ks1gamma2)       = par_meta('ks1gamma2'      , 0.0126_dp ,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%ks1gamma3)       = par_meta('ks1gamma3'      ,-0.0064_dp ,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%ks2gamma1)       = par_meta('ks2gamma1'      ,   54.0_dp ,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%ks2gamma2)       = par_meta('ks2gamma2'      ,  -0.07_dp ,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%ks2gamma3)       = par_meta('ks2gamma3'      , -0.167_dp ,-0.8_dp ,-0.4_dp, .False.)
    !pororsity transfer function
    gammaMaster(ixGamma%phi1gamma1)      = par_meta('phi1gamma1'     ,   50.5_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi1gamma2)      = par_meta('phi1gamma2'     , -0.142_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi1gamma3)      = par_meta('phi1gamma3'     , -0.037_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi2gamma1)      = par_meta('phi2gamma1'     ,   0.76_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi2gamma2)      = par_meta('phi2gamma2'     , 0.0009_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi2gamma3)      = par_meta('phi2gamma3'     , -0.264_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi2gamma4)      = par_meta('phi2gamma4'     ,   0.89_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi2gamma5)      = par_meta('phi2gamma5'     , -0.001_dp, -0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%phi2gamma6)      = par_meta('phi2gamma6'     , -0.324_dp, -0.8_dp ,-0.4_dp, .False.)
    !field capacity transfer function
    gammaMaster(ixGamma%fc1gamma1)       = par_meta('fc1gamma1'      ,    1.0_dp, 0.8_dp ,  1.2_dp, .False.)
    !wilting point transfer function
    gammaMaster(ixGamma%wp1gamma1)       = par_meta('wp1gamma1'      ,    1.0_dp, 0.8_dp ,  1.2_dp, .False.)
    !b transfer function
    gammaMaster(ixGamma%b1gamma1)        = par_meta('b1gamma1'       ,    3.1_dp,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%b1gamma2)        = par_meta('b1gamma2'       ,  0.157_dp,-0.8_dp ,-0.4_dp, .False.)
    gammaMaster(ixGamma%b1gamma3)        = par_meta('b1gamma3'       , -0.003_dp,-0.8_dp ,-0.4_dp, .False.)
    !saturation matric potential transfer function
    gammaMaster(ixGamma%psis1gamma1)     = par_meta('psis1gamma1'    ,   1.54_dp,  0.8_dp ,  0.4_dp, .False.)
    gammaMaster(ixGamma%psis1gamma2)     = par_meta('psis1gamma2'    ,-0.0095_dp, -0.8_dp , -0.4_dp, .False.)
    gammaMaster(ixGamma%psis1gamma3)     = par_meta('psis1gamma3'    , 0.0063_dp,  0.8_dp ,  0.4_dp, .False.)
    !specific yield transfer function
    gammaMaster(ixGamma%myu1gamma1)      = par_meta('myu1gamma1'     ,    3.5_dp,  0.8_dp ,  0.4_dp, .False.)
    gammaMaster(ixGamma%myu1gamma2)      = par_meta('myu1gamma2'     ,   1.66_dp,  0.8_dp ,  0.4_dp, .False.)
    ! total depth multiplier
    gammaMaster(ixGamma%z1gamma1)        = par_meta('z1gamma1'       ,    1.0_dp,  0.1_dp ,  0.4_dp, .False.)
    ! layer fractions
    gammaMaster(ixGamma%h1gamma1)        = par_meta('h1gamma1'       ,   0.05_dp, 0.01_dp ,  0.1_dp, .False.)
    gammaMaster(ixGamma%h1gamma2)        = par_meta('h1gamma2'       ,    0.3_dp, 0.12_dp ,  0.5_dp, .False.)
    ! transfer function
    gammaMaster(ixGamma%binfilt1gamma1)  = par_meta('binfilt1gamma1' ,     0.0_dp, -2.0_dp , 1.0_dp, .False.)
    gammaMaster(ixGamma%binfilt1gamma2)  = par_meta('binfilt1gamma1' ,    -0.6_dp,  0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%D11gamma1)       = par_meta('D11gamma1'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%D21gamma1)       = par_meta('D21gamma1'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%D31gamma1)       = par_meta('D31gamma1'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%D41gamma1)       = par_meta('D41gamma1'      ,     2.0_dp , 1.2_dp , 2.5_dp, .False.)
    gammaMaster(ixGamma%exp1gamma1)      = par_meta('exp1gamma1'     ,     3.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%exp1gamma1)      = par_meta('exp1gamma1'     ,     2.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%ksat1gamma1)     = par_meta('ksat1gamma1'    ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%bbl1gamma1)      = par_meta('bbl1gamma1'     ,    0.32_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%bbl1gamma2)      = par_meta('bbl1gamma2'     ,     4.2_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%BD1gamma1)       = par_meta('BD1gamma2'      ,     1.0_dp , 0.9_dp , 1.1_dp, .False.)
    gammaMaster(ixGamma%SD1gamma1)       = par_meta('SD1gamma2'      ,     1.0_dp , 0.9_dp , 1.1_dp, .False.)
    gammaMaster(ixGamma%WcrFrac1gamma1)  = par_meta('WcrFrac1gamma1' ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    gammaMaster(ixGamma%WpwpFrac1gamma1) = par_meta('WpwpFrac1gamma1',     1.0_dp , 0.8_dp , 1.2_dp, .False.)

    ! -----
    !  Master list of beta parameters  
    ! -----------------------
    betaMaster(ixBeta%h1)                = par_meta('h1'      ,     1.0_dp,  0.8_dp , 1.0_dp, .False.)
    betaMaster(ixBeta%h2)                = par_meta('h2'      ,     1.0_dp,  0.8_dp , 1.0_dp, .False.)
    betaMaster(ixBeta%h3)                = par_meta('h3'      ,     1.0_dp,  0.8_dp , 1.0_dp, .False.)
    betaMaster(ixBeta%h4)                = par_meta('h4'      ,     1.0_dp,  0.8_dp , 1.0_dp, .False.)
    betaMaster(ixBeta%h5)                = par_meta('h5'      ,     1.0_dp,  0.8_dp , 1.0_dp, .False.)
    betaMaster(ixBeta%binfilt)           = par_meta('binfilt' ,     1.0_dp,  0.8_dp , 1.0_dp, .False.)
    betaMaster(ixBeta%D1)                = par_meta('D1'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%D2)                = par_meta('D2'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%D3)                = par_meta('D3'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%D4)                = par_meta('D4'      ,     1.0_dp , 0.8_dp , 2.5_dp, .False.)
    betaMaster(ixBeta%expt)              = par_meta('expt'    ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%ks)                = par_meta('ks'      ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%bbl)               = par_meta('bbl'     ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%BD)                = par_meta('BD'      ,     1.0_dp , 0.9_dp , 1.1_dp, .False.)
    betaMaster(ixBeta%SD)                = par_meta('SD'      ,     1.0_dp , 0.9_dp , 1.1_dp, .False.)
    betaMaster(ixBeta%WcrFrac)           = par_meta('WcrFrac' ,     1.0_dp , 0.8_dp , 1.2_dp, .False.)
    betaMaster(ixBeta%WpwpFrac)          = par_meta('WpwpFrac',     1.0_dp , 0.8_dp , 1.2_dp, .False.)

  end subroutine popMeta

end module paramMaster 
