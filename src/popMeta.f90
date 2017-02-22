module popMeta

implicit none

private

public::paramMaster
public::popMprMeta

contains

subroutine paramMaster(err,message)
  use nrtype
  use data_type,  only:par_meta
  use var_lookup, only:ixPar
  use globalData, only:parMaster
  implicit none
  !output variable
  integer(i4b),intent(out)      :: err     ! error code
  character(*),intent(out)      :: message ! error message
  
  ! initialize error control
  err=0; message='popMeta/'
  ! -----
  !  Public subroutine: Master list of gamma parameters  
  ! -----------------------
  !                                                        name,    default, lwr bound, upr bound,parent beta,   TF,    type,    mask, h-scale,     p-norm,     v scale,    p-norm,  perLyr
  !ks transfer function
  parMaster(ixPar%ks1gamma1)       = par_meta('ks1gamma1'      ,   -0.6_dp ,  -0.66_dp,  -0.54_dp,"ks"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%ks1gamma2)       = par_meta('ks1gamma2'      , 0.0126_dp , 0.0113_dp, 0.0139_dp,"ks"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%ks1gamma3)       = par_meta('ks1gamma3'      ,-0.0064_dp ,-0.0058_dp,-0.0070_dp,"ks"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%ks2gamma1)       = par_meta('ks2gamma1'      ,   54.0_dp ,   48.6_dp,   59.4_dp,"ks"       ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%ks2gamma2)       = par_meta('ks2gamma2'      ,  -0.07_dp , -0.077_dp, -0.063_dp,"ks"       ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%ks2gamma3)       = par_meta('ks2gamma3'      , -0.167_dp ,   -0.8_dp,   -0.4_dp,"ks"       ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  !pororsity transfer function
  parMaster(ixPar%phi1gamma1)      = par_meta('phi1gamma1'     ,    50.5_dp,   45.5_dp,   55.5_dp,"phi"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi1gamma2)      = par_meta('phi1gamma2'     ,  -0.142_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi1gamma3)      = par_meta('phi1gamma3'     ,  -0.037_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi2gamma1)      = par_meta('phi2gamma1'     ,    0.76_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi2gamma2)      = par_meta('phi2gamma2'     ,  0.0009_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi2gamma3)      = par_meta('phi2gamma3'     ,  -0.264_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi2gamma4)      = par_meta('phi2gamma4'     ,    0.89_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi2gamma5)      = par_meta('phi2gamma5'     ,  -0.001_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%phi2gamma6)      = par_meta('phi2gamma6'     ,  -0.324_dp,   -0.8_dp,   -0.4_dp,"phi"      ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  !field capacity transfer function
  parMaster(ixPar%fc1gamma1)       = par_meta('fc1gamma1'      ,     1.0_dp,    0.8_dp,    1.2_dp,"fc"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  !wilting point transfer function
  parMaster(ixPar%wp1gamma1)       = par_meta('wp1gamma1'      ,     1.0_dp,    0.8_dp,    1.2_dp,"wp"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  !b transfer function
  parMaster(ixPar%b1gamma1)        = par_meta('b1gamma1'       ,     3.1_dp,   -0.8_dp,   -0.4_dp,"b"        ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%b1gamma2)        = par_meta('b1gamma2'       ,   0.157_dp,   -0.8_dp,   -0.4_dp,"b"        ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%b1gamma3)        = par_meta('b1gamma3'       ,  -0.003_dp,   -0.8_dp,   -0.4_dp,"b"        ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  !saturation matric potential transfer function
  parMaster(ixPar%psis1gamma1)     = par_meta('psis1gamma1'    ,    1.54_dp,    0.8_dp,    0.4_dp,"psis"     ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%psis1gamma2)     = par_meta('psis1gamma2'    , -0.0095_dp,   -0.8_dp,   -0.4_dp,"psis"     ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%psis1gamma3)     = par_meta('psis1gamma3'    ,  0.0063_dp,    0.8_dp,    0.4_dp,"psis"     ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  !specific yield transfer function
  parMaster(ixPar%myu1gamma1)      = par_meta('myu1gamma1'     ,     3.5_dp,   3.85_dp,   3.15_dp,"myu"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%myu1gamma2)      = par_meta('myu1gamma2'     ,    1.66_dp,   1.83_dp,   1.50_dp,"myu"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  ! total depth multiplier
  parMaster(ixPar%z1gamma1)        = par_meta('z1gamma1'       ,     1.0_dp,    0.1_dp,    4.0_dp,"z"        ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  ! layer fractions
  parMaster(ixPar%h1gamma1)        = par_meta('h1gamma1'       ,    0.05_dp,   0.01_dp,    0.1_dp,"h1"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%h2gamma1)        = par_meta('h2gamma1'       ,     0.3_dp,   0.12_dp,    0.5_dp,"h2"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  ! transfer function
  parMaster(ixPar%binfilt1gamma1)  = par_meta('binfilt1gamma1' ,     0.0_dp,   -2.0_dp,    1.0_dp,"binfilt"  ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%binfilt1gamma2)  = par_meta('binfilt1gamma2' ,     1.0_dp,    0.8_dp,    1.2_dp,"binfilt"  ,     2, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%D11gamma1)       = par_meta('D11gamma1'      ,     2.0_dp,   1.75_dp,    3.5_dp,"D1"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%D21gamma1)       = par_meta('D21gamma1'      ,     2.0_dp,   1.75_dp,    3.5_dp,"D2"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%D31gamma1)       = par_meta('D31gamma1'      ,     1.0_dp,  0.001_dp,    2.0_dp,"D3"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%D41gamma1)       = par_meta('D41gamma1'      ,     2.0_dp,    1.2_dp,    2.5_dp,"D4"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%exp1gamma1)      = par_meta('exp1gamma1'     ,     3.0_dp,    0.8_dp,    1.2_dp,"exp"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%exp1gamma2)      = par_meta('exp1gamma2'     ,     2.0_dp,    0.8_dp,    1.2_dp,"exp"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%bbl1gamma1)      = par_meta('bbl1gamma1'     ,    0.32_dp,    0.8_dp,    1.2_dp,"bbl"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%bbl1gamma2)      = par_meta('bbl1gamma2'     ,     4.2_dp,    0.8_dp,    1.2_dp,"bbl"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%bd1gamma1)       = par_meta('bd1gamma1'      ,     1.0_dp,    0.9_dp,    1.1_dp,"bd"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%sd1gamma1)       = par_meta('sd1gamma1'      ,     1.0_dp,    0.9_dp,    1.1_dp,"sd"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%WcrFrac1gamma1)  = par_meta('WcrFrac1gamma1' ,     1.0_dp,    0.8_dp,    1.2_dp,"WcrFrac"  ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%WpwpFrac1gamma1) = par_meta('WpwpFrac1gamma1',     1.0_dp,    0.8_dp,    1.2_dp,"WcrFrac"  ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%fsm1gamma1)      = par_meta('fsm1gamma1',          1.0_dp,    0.8_dp,    1.2_dp,"fwm"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%zk1gamma1)       = par_meta('zk1gamma1',           1.6_dp,    1.4_dp,    1.8_dp,"zk"       ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%zsk1gamma1)      = par_meta('zsk1gamma1',          1.6_dp,    1.4_dp,    1.8_dp,"zsk"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%zsk1gamma2)      = par_meta('zsk1gamma2',          2.0_dp,    1.8_dp,    2.2_dp,"zsk"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%zpk1gamma1)      = par_meta('zpk1gamma1',          1.0_dp,    0.8_dp,    1.2_dp,"zpk"      ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%pfree1gamma1)    = par_meta('pfree1gamma1',        1.6_dp,    1.4_dp,    1.8_dp,"pfree"    ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  parMaster(ixPar%rexp1gamma1)     = par_meta('rexp1gamma1',        0.03_dp,   0.02_dp,   0.04_dp,"rexp"     ,     1, "soil", .False.,      "na", -999.0_dp,       "na", -999.0_dp,.False.)
  ! -----
  !  Master list of beta parameters  
  ! -----------------------
  !                                                        name,    default, lwr bound, upr bound,parent beta,    TF,    type,    mask, h-scale,    p-norm,   v scale,    p-norm,  perLyr
  parMaster(ixPar%uhshape)         = par_meta('uhshape'        ,     1.0_dp,    0.1_dp,    3.0_dp,     "beta",  -999, "route", .False.,     "na", -999.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%uhscale)         = par_meta('uhscale'        ,     1.0_dp,    0.5_dp,    3.0_dp,     "beta",  -999, "route", .False.,     "na", -999.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%ks)              = par_meta('ks'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .True.)
  parMaster(ixPar%bd)              = par_meta('bd'             ,     1.0_dp,    0.9_dp,    1.1_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%sd)              = par_meta('sd'             ,     1.0_dp,    0.9_dp,    1.1_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%psis)            = par_meta('psis'           ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%b)               = par_meta('b'              ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%phi)             = par_meta('phi'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%fc)              = par_meta('fc'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%wp)              = par_meta('wp'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%myu)             = par_meta('myu'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%binfilt)         = par_meta('binfilt'        ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .False.)
  parMaster(ixPar%D1)              = par_meta('D1'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .False.)
  parMaster(ixPar%D2)              = par_meta('D2'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .False.)
  parMaster(ixPar%D3)              = par_meta('D3'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .False.)
  parMaster(ixPar%D4)              = par_meta('D4'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .False.)
  parMaster(ixPar%c)               = par_meta('c'              ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .False.)
  parMaster(ixPar%Dsmax)           = par_meta('Dsmax'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .False.)
  parMaster(ixPar%Ds)              = par_meta('Ds'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .False.)
  parMaster(ixPar%Ws)              = par_meta('Ws'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .False.)
  parMaster(ixPar%expt)            = par_meta('expt'           ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%bbl)             = par_meta('bbl'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%h1)              = par_meta('h1'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .True.)
  parMaster(ixPar%h2)              = par_meta('h2'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .True.)
  parMaster(ixPar%h3)              = par_meta('h3'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .True.)
  parMaster(ixPar%h4)              = par_meta('h4'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .True.)
  parMaster(ixPar%h5)              = par_meta('h5'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .True.)
  parMaster(ixPar%z)               = par_meta('z'              ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .True.)
  parMaster(ixPar%WcrFrac)         = par_meta('WcrFrac'        ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%WpwpFrac)        = par_meta('WpwpFrac'       ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%twm)             = par_meta('twm'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%fwm)             = par_meta('fwm'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%fsm)             = par_meta('fsm'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%fpm)             = par_meta('fpm'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%zk)              = par_meta('zk'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%zsk)             = par_meta('zsk'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .True.)
  parMaster(ixPar%zpk)             = par_meta('zpk'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",   -1.0_dp,   "pnorm",   -1.0_dp, .True.)
  parMaster(ixPar%pfree)           = par_meta('pfree'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%zperc)           = par_meta('zperc'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%rexp)            = par_meta('rexp'           ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "soil", .False.,  "pnorm",    1.0_dp,   "pnorm",    1.0_dp, .True.)
  parMaster(ixPar%rmin)            = par_meta('rmin'           ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "vege", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%lai)             = par_meta('lai'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "vege", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%scf)             = par_meta('scf'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%mfmax)           = par_meta('mfmax'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%mfmin)           = par_meta('mfmin'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%uadj)            = par_meta('uadj'           ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%si)              = par_meta('si'             ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%pxtemp)          = par_meta('pxtemp'         ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%nmf)             = par_meta('nmf'            ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%tipm)            = par_meta('tipm'           ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%plwhc)           = par_meta('plwhc'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)
  parMaster(ixPar%daygm)           = par_meta('daygm'          ,     1.0_dp,    0.8_dp,    1.2_dp,     "beta",  -999,  "snow", .False.,  "pnorm",    1.0_dp,      "na", -999.0_dp, .False.)

end subroutine 

! -----
!  Public subroutine: Populate metadata for MPR infor 
! -----------------------
subroutine popMprMeta(err,message)
  use nrtype
  use data_type,  only:var_meta
  use var_lookup, only:ixVarMapData
  use var_lookup, only:ixVarSoilData
  use var_lookup, only:ixVarVegData
  use globalData, only:map_meta
  use globalData, only:sdata_meta
  use globalData, only:vdata_meta
  
  implicit none
  !output variable
  integer(i4b),intent(out)      :: err     ! error code
  character(*),intent(out)      :: message ! error message
  
  ! initialize error control
  err=0; message='popMprMeta/'
  ! Mapping data meta
  map_meta(ixVarMapData%hru_id)          = var_meta('hru_id'       ,"hru id"                            , "-"              ,"1D", "integer")
  map_meta(ixVarMapData%weight)          = var_meta('weight'       ,"areal weight of intersect polygon" , "-"              ,"1D", "integer")
  map_meta(ixVarMapData%overlapPolyId)   = var_meta('overlapPolyId',"id of intersect polygo"            , "-"              ,"1D", "integer")
  ! Soil data variables
  sdata_meta(ixVarSoilData%polyid)       = var_meta('polyid'       , "soil polygon id"                  ,  "-"             ,"1D", "integer")
  sdata_meta(ixVarSoilData%soilclass)    = var_meta('soilclass'    , "USDA soil class"                  ,  "-"             ,"2D", "integer")
  sdata_meta(ixVarSoilData%hslyrs)       = var_meta('hslyrs'       , "soil layer thickness"             ,  "m"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%sand_frc)     = var_meta('sand_frc'     , "sand percentage"                  ,  "%"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%silt_frc)     = var_meta('silt_frc'     , "silt percentage"                  ,  "%"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%clay_frc)     = var_meta('clay_frc'     , "clay percentage"                  ,  "%"             ,"2D", "double" )
  sdata_meta(ixVarSoilData%bulk_density) = var_meta('bulk_density' , "bulk density"                     ,  "kg/m^3"        ,"2D", "double" )
  sdata_meta(ixVarSoilData%ele_mean)     = var_meta('ele_mean'     , "mean elevation"                   ,  "m"             ,"1D", "double" )
  sdata_meta(ixVarSoilData%ele_std)      = var_meta('ele_std'      , "std elevation"                    ,  "m"             ,"1D", "double" )
  sdata_meta(ixVarSoilData%slp_mean)     = var_meta('slp_mean'     , "mean slope"                       ,  "-"             ,"1D", "double" )
  ! Vege data variables 
  vdata_meta(ixVarVegData%polyid)        = var_meta('polyid'       , "vege polygon id"                  ,  "-"             ,"1D", "integer")
  vdata_meta(ixVarVegData%vegclass)      = var_meta('vegclass'     , "vegetation class"                 ,  "-"             ,"1D", "double" )
  vdata_meta(ixVarVegData%grnfrc)        = var_meta('grnfrc'       , "green fraction"                   ,  "-"             ,"1D", "double" )
  vdata_meta(ixVarVegData%lai)           = var_meta('lai'          , "monthly lai"                      ,  "m^2/m^2"       ,"2D", "double" )

end subroutine

end module popMeta 
