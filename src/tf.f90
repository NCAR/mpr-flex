module tf
! Library of model parameter transfer function 
use nrtype                                        ! variable types, etc.
use data_type                                     ! Including custum data structure definition
use public_var                                    ! Including common constant (physical constant, other e.g., missingVal, etc.)
use var_lookup, only:ixVarSoilData, ixVarVegData, ixBeta, ixGamma, nBeta

implicit none

private

public :: comp_model_param
public :: betaDependency
public :: betaCompOrder 

contains

! ********************************************************************************************
! Public subroutine: Execute computation of model-dependant soil parameter transfer function  
! *********************************************************************************************
subroutine comp_model_param(parSxySz,          &  ! in/output: soil parameter values for all L0 polygons that are included
                            parVxy,            &  ! in/output: veg parameter values for all L0 polygons that are included 
                            sdata,             &  ! input: soil data
                            vdata,             &  ! input: vege data
                            gammaParMasterMeta,&  ! input: gamma parameter meta file - val of calibrating parammeter is adjusted via calibration 
                            nSLyr,             &  ! input: number of soil layers
                            nSPoly,            &  ! input: number of soil polygons
                            nVPoly,            &  ! input: number of vege polygons (to be deleted as nVpoly=nSpoly)
                            err,message) 

  use globalData, only:betaMaster, betaOrder, soilBetaInGamma, vegBetaInGamma
  use get_ixname, only:get_ixBeta
  implicit none
  ! in/out
  type(namedvar2),      intent(inout) :: parSxySz(:)            ! soil parameter values for ParSxySz(:)%varDat(lyr,poly) 
  type(namedvar2),      intent(inout) :: parVxy(:)              ! veg parameter values for ParVxy(:)%varDat(lyr,poly) 
  ! input
  type(namevar),        intent(in)    :: sdata(:)               ! storage of soil data strucuture
  type(namevar),        intent(in)    :: vdata(:)               ! storage of veg data strucuture
  type(par_meta)                      :: gammaParMasterMeta(:)
  integer(i4b),         intent(in)    :: nSLyr                  ! number of soil layer
  integer(i4b),         intent(in)    :: nSPoly                 ! number of soil polygons
  integer(i4b),         intent(in)    :: nVPoly                 ! number of vege polygons (to be deleted)
  ! output
  integer(i4b),         intent(out)   :: err                    ! error code
  character(len=strLen),intent(out)   :: message                ! error message for current routine
  ! Local 
  type(namedvar2)                     :: parTemp(nBeta)         ! soil parameter values for ParSxySz(:)%dat(lyr,poly) 
  integer(i4b)                        :: ix                     ! index of gamma parameter 
  integer(i4b)                        :: idBeta                 ! id of beta parameter array 
  integer(i4b)                        :: iParm                  ! Loop index of model parameters (e.g., VIC)

  err=0; message="comp_model_param/"
  first: associate(gammaPar=> gammaParMasterMeta(:)%val)
  do iParm = 1,size(betaOrder)
    ix = betaOrder(iParm) 
    if (ix/=-999) then
      if (trim(betaMaster(ix)%ptype)=='soil')then
        allocate(parTemp(ix)%varData(nSLyr,nSPoly) ,stat=err); if(err/=0)then;message=trim(message)//'error allocating parTemp';stop;endif
      elseif (betaMaster(ix)%ptype=='veg')then 
        allocate(parTemp(ix)%varData(nMonth,nVPoly) ,stat=err); if(err/=0)then;message=trim(message)//'error allocating parTemp';stop;endif
      endif
      second: associate (xPar => parTemp(ix)%varData, &
                         tfid => betaMaster(ix)%tftype)
      if (tfid==-999_i4b) tfid=1_i4b
      select case(ix)
        case(ixBeta%ks)
          call ks( err, message, sdata=sdata, gammaPar=gammaPar, ks_out=xPar, opt=tfid )
        case(ixBeta%bd)
          call bd( err, message, sdata=sdata, gammaPar=gammaPar, bd_out=xPar, opt=tfid )
        case(ixBeta%phi)
          call phi( err, message, sdata=sdata, gammaPar=gammaPar, phi_out=xPar, opt=tfid ) 
        case(ixBeta%b)
          call retcurve( err, message, sdata=sdata, gammaPar=gammaPar, retcurve_out=xPar, opt=tfid )
        case(ixBeta%psis)
          call psis( err, message, sdata=sdata, gammaPar=gammaPar, psis_out=xPar, opt=tfid )
        case(ixBeta%fc)
          call fc( err, message, sdata=sdata, phi_in=parTemp(ixBeta%phi)%varData, psis_in=parTemp(ixBeta%psis)%varData, b_in=parTemp(ixBeta%b)%varData, gammaPar=gammaPar, fc_out=xPar, opt=tfid )
        case(ixBeta%wp)
          call wp( err, message, phi_in=parTemp(ixBeta%phi)%varData, psis_in=parTemp(ixBeta%psis)%varData, b_in=parTemp(ixBeta%b)%varData, gammaPar=gammaPar, wp_out=xPar, opt=tfid ) 
        case(ixBeta%myu)
          call myu( err, message, phi_in=parTemp(ixBeta%phi)%varData, fc_in=parTemp(ixBeta%fc)%varData, gammaPar=gammaPar, myu_out=xPar, opt=tfid ) 
        case(ixBeta%binfilt)
          call binfilt(err, message, sdata=sdata, gammaPar=gammaPar, binfilt_out=xPar, opt=tfid )
        case(ixBeta%D1)
          call D1( err, message, sdata=sdata, ks_in=parTemp(ixBeta%ks)%varData, phi_in=parTemp(ixBeta%phi)%varData, gammaPar=gammaPar, D1_out=xPar, opt=tfid ) 
        case(ixBeta%D2)
          call D2( err, message, sdata=sdata, ks_in=parTemp(ixBeta%ks)%varData, D4_in=parTemp(ixBeta%D4)%varData, gammaPar=gammaPar, D2_out=xPar, opt=tfid ) 
        case(ixBeta%D3)
          call D3( err, message, sdata=sdata, fc_in=parTemp(ixBeta%fc)%varData, gammaPar=gammaPar, D3_out=xPar, opt=tfid ) 
        case(ixBeta%D4)
          call D4( err, message, gammaPar=gammaPar, D4_out=xPar, opt=tfid ) 
        case(ixBeta%Ds)
          call Ds( err, message, D1_in=parTemp(ixBeta%D1)%varData, D3_in=parTemp(ixBeta%D3)%varData, Dsmax_in=parTemp(ixBeta%Dsmax)%varData, gammaPar=gammaPar, Ds_out=xPar, opt=tfid )
        case(ixBeta%c)
          call cexpt( err, message, D4_in=parTemp(ixBeta%D4)%varData, gammaPar=gammaPar, cexpt_out=xPar, opt=tfid )
        case(ixBeta%sd)
          call sd( err, message, gammaPar=gammaPar, sd_out=xPar, opt=tfid )
        case(ixBeta%expt)
          call expt( err, message, b_in=parTemp(ixBeta%b)%varData, gammaPar=gammaPar, expt_out=xPar, opt=tfid )
        case(ixBeta%Dsmax)
          call Dsmax( err, message, &
                      sdata=sdata, D1_in=parTemp(ixBeta%D1)%varData, D2_in=parTemp(ixBeta%D2)%varData, D3_in=parTemp(ixBeta%D3)%varData, c_in=parTemp(ixBeta%c)%varData, phi_in=parTemp(ixBeta%phi)%varData, &
                      gammaPar=gammaPar, Dsmax_out=xPar ,opt=tfid )          
        case(ixBeta%bbl)
          call bubble( err, message, expt_in=parTemp(ixBeta%expt)%varData, gammaPar=gammaPar, bubble_out=xPar, opt=tfid )
        case(ixBeta%WcrFrac)
          call WcrFrac( err, message, fc_in=parTemp(ixBeta%fc)%varData, phi_in=parTemp(ixBeta%phi)%varData, gammaPar=gammaPar, WcrFrac_out=xPar, opt=tfid )
        case(ixBeta%WpwpFrac)
          call WpwpFrac( err, message, wp_in=parTemp(ixBeta%wp)%varData, phi_in=parTemp(ixBeta%phi)%varData, gammaPar=gammaPar, WpwpFrac_out=xPar, opt=tfid )
        case(ixBeta%Ws)
          call Ws( err, message, sdata=sdata, D3_in=parTemp(ixBeta%D3)%varData, phi_in=parTemp(ixBeta%phi)%varData, gammaPar=gammaPar, Ws_out=xPar, opt=tfid)          
        case(ixBeta%twm)
          call twm( err, message, sdata=sdata, fc_in=parTemp(ixBeta%fc)%varData, wp_in=parTemp(ixBeta%wp)%varData, gammaPar=gammaPar, twm_out=xPar, opt=tfid )
        case(ixBeta%fwm)
          call fwm( err, message, sdata=sdata, phi_in=parTemp(ixBeta%phi)%varData, fc_in=parTemp(ixBeta%fc)%varData, gammaPar=gammaPar, fwm_out=xPar, opt=tfid )
        case(ixBeta%fsm)
          call fsm( err, message, fwm_in=parTemp(ixBeta%fwm)%varData, phi_in=parTemp(ixBeta%phi)%varData, wp_in=parTemp(ixBeta%wp)%varData, gammaPar=gammaPar, fsm_out=xPar, opt=tfid )
        case(ixBeta%fpm)
          call fpm( err, message, fwm_in=parTemp(ixBeta%fwm)%varData, fsm_in=parTemp(ixBeta%fsm)%varData, gammaPar=gammaPar, fpm_out=xPar, opt=tfid )
        case(ixBeta%zk)
          call zk( err, message, phi_in=parTemp(ixBeta%phi)%varData, fc_in=parTemp(ixBeta%fc)%varData, gammaPar=gammaPar, zk_out=xPar, opt=tfid )
        case(ixBeta%zsk)
          call zsk( err, message, phi_in=parTemp(ixBeta%phi)%varData, fc_in=parTemp(ixBeta%fc)%varData, wp_in=parTemp(ixBeta%wp)%varData, gammaPar=gammaPar, zsk_out=xPar, opt=tfid )
        case(ixBeta%zpk)
          call zpk( err, message, sdata=sdata, ks_in=parTemp(ixBeta%ks)%varData, myu_in=parTemp(ixBeta%myu)%varData, gammaPar=gammaPar, zpk_out=xPar, opt=tfid )
        case(ixBeta%pfree)
          call pfree( err, message, phi_in=parTemp(ixBeta%phi)%varData, wp_in=parTemp(ixBeta%wp)%varData, gammaPar=gammaPar, pfree_out=xPar, opt=tfid )
        case(ixBeta%zperc)
          call zperc( err, message, twm_in=parTemp(ixBeta%twm)%varData, fsm_in=parTemp(ixBeta%fsm)%varData, zsk_in=parTemp(ixBeta%zsk)%varData, fpm_in=parTemp(ixBeta%fpm)%varData, zpk_in=parTemp(ixBeta%zsk)%varData, gammaPar=gammaPar, zperc_out=xPar, opt=tfid )
        case(ixBeta%rexp)
          call rexp( err, message, wp_in=parTemp(ixBeta%wp)%varData, gammaPar=gammaPar, rexp_out=xPar, opt=tfid )
        case(ixBeta%lai)
          call lai( err, message, vdata=vdata, gammaPar=gammaPar, lai_out=xPar, opt=tfid ) 
      end select ! end of parameter case
      end associate second
    endif
  end do ! end of parameter loop
  end associate first
  ! extract beta parameters in 'CalPar' list
  do iParm=1,size(soilBetaInGamma)
    idBeta=get_ixBeta(trim(soilBetaInGamma(iParm)))
    parSxySz(iParm)%varData=parTemp(idBeta)%varData 
  enddo
  do iParm=1,size(vegBetaInGamma)
    idBeta=get_ixBeta(trim(vegBetaInGamma(iParm)))
    parVxy(iParm)%varData=parTemp(idBeta)%varData 
  enddo
  return
end subroutine 

! ***********
! public subroutine: Get list of dependant parameter ids given a parameter 
! *********************************************************************
subroutine betaDependency( err, message )
  use get_ixname, only:get_ixBeta
  use globalData, only:betaMaster, beta
  implicit none
  ! output 
  integer(i4b),         intent(out)   :: err                     ! output: error id 
  character(len=strLen),intent(out)   :: message                 ! output: error message   
  ! local 
  integer(i4b),allocatable            :: ixDepend(:)             ! id list of dependent beta parameters 
  integer(i4b)                        :: iParm                   ! loop index 
  character(len=strLen)               :: cmessage                ! error message from subroutine

  err=0; message="parDependcy/"
  do iParm=1,nBeta
    select case(iParm)
      case(ixBeta%ks);      call ks      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%bd);      call bd      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%sd);      call sd      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%psis);    call psis    (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif 
      case(ixBeta%b);       call retcurve(err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif 
      case(ixBeta%phi);     call phi     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif 
      case(ixBeta%fc);      call fc      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%wp);      call wp      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%myu);     call myu     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%binfilt); call binfilt (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%D1);      call D1      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%D2);      call D2      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%D3);      call D3      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%D4);      call D4      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%Ds);      call Ds      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%Dsmax);   call Dsmax   (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%Ws);      call Ws      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%c);       call cexpt   (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%expt);    call expt    (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%bbl);     call bubble  (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%WcrFrac); call WcrFrac (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%WpwpFrac);call WpwpFrac(err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%twm);     call twm     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%fwm);     call fwm     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%fsm);     call fsm     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%fpm);     call fpm     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%zk);      call zk      (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%zsk);     call zsk     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%zpk);     call zpk     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%pfree);   call pfree   (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%zperc);   call zperc   (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%rexp);    call rexp    (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%lai);     call lai     (err, cmessage, ixDepend=ixDepend); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      case(ixBeta%z); allocate(ixDepend(1), stat=err); ixDepend=-999_i4b
      case(ixBeta%h1);allocate(ixDepend(1), stat=err); ixDepend=-999_i4b
      case(ixBeta%h2);allocate(ixDepend(1), stat=err); ixDepend=-999_i4b
      case(ixBeta%h3);allocate(ixDepend(1), stat=err); ixDepend=-999_i4b
      case(ixBeta%h4);allocate(ixDepend(1), stat=err); ixDepend=-999_i4b
      case(ixBeta%h5);allocate(ixDepend(1), stat=err); ixDepend=-999_i4b
    end select
    if ( allocated(ixDepend) )then
      allocate(beta(iParm)%depend(size(ixDepend)),stat=err);if(err/=0)then;message=trim(message)//'error allocating beta%ixDepend for '//trim(betaMaster(iParm)%pname);return;endif 
      beta(iParm)%depend = ixDepend
      deallocate(ixDepend)
    endif
  enddo
  return
end subroutine

! ***********
! public subroutine: Determin computing order of beta parameters including depedent parameters 
! *********************************************************************
subroutine betaCompOrder( betaList, err, message)
  use globalData, only:beta, betaOrder, nBetaNeed
  use get_ixname, only:get_ixBeta
  implicit none
  ! input
  character(len=strLen)             :: betaList(:)        ! output: error message   
  ! output 
  integer(i4b),         intent(out) :: err                ! output: error id 
  character(len=strLen),intent(out) :: message            ! output: error message   
  ! local
  integer(i4b)                      :: iParm,jParm,kParm  ! loop index for betar parameter 
  logical(lgc)                      :: parFlag(nBeta)     ! logical array to check order assignment 
  integer(i4b)                      :: nassign            ! counter for order assigned 
  integer(i4b)                      :: idx,jdx            ! parameter index 
  integer(i4b)                      :: iDeps              ! loop index for dependent parameters 
  integer(i4b)                      :: nDeps              ! number of dependant parameters of the parameter 

  err=0; message="betaOrder/"
  parFlag = .false.
  nassign = 0
  betaOrder=-999
  nBetaNeed=0
  do  ! do until all beta parameters are assigned
    nassign = 0
    do iParm=1,size(betaList)
      jdx=get_ixBeta(betaList(iParm))
      ! check if the parameter is assigned yet
      if(parFlag(jdx))then; nassign = nassign + 1; cycle; endif
      ! climb up parameter network as far as possible
      jParm = jdx    ! the first parameter under investigation 
      do  ! do until get to a "most uplevel" parameters that is not assigned
        if ( beta(jParm)%depend(1) .eq. -999_i4b ) then     ! (if nDeps = 0, then the parameter is independent of any others )
          ! assign jParm
          nBetaNeed=nBetaNeed+1
          parFlag(jParm) = .true.
          beta(jParm)%order = nBetaNeed 
          betaOrder(nBetaNeed)=jParm
          exit 
        else    ! if the parameter has any dependent parameters 
          kParm = jParm   ! the parameter under investigation 
          nDeps = size(beta(jParm)%depend)    ! get number of number of dependent parameters
          do iDeps=1,nDeps
            idx = beta(jParm)%depend(iDeps)   
            ! check if the parameter is NOT assigned
            if (.not.parFlag(idx))then; jParm = idx; exit; endif
          end do  ! (looping through dependent parameters)
          ! check if all upstream reaches are already assigned (only if kParm=jParm)
          if (jParm .eq. kParm) then 
            ! assign jParm
            nBetaNeed=nBetaNeed+1
            parFlag(jParm) = .true.
            beta(jParm)%order = nBetaNeed
            betaOrder(nBetaNeed)=jParm
            exit 
          endif  
          cycle   ! if jrch changes, keep looping (move upstream)
        endif 
      end do  !  climbing upstream (do-forever)
    end do   ! looping through beta parameters 
    if (nassign.eq.size(betaList)) exit 
  end do  ! do forever (do until all parameters are assigned)
  return
end subroutine

! *********************************************************************
! Library of hydrologic model transfer functions
! *********************************************************************
! infilt parameter 
! *********************************************************************
subroutine binfilt(err, message, ixDepend, sdata, gammaPar, binfilt_out, opt )
  ! Use Arno scheme 
  ! b = (sigma_elev-sigma_min)/(sigma_ele+sigma_max)
  ! where
  ! sigma_ele: standard deviation of elevation in polygon
  ! sigma_min: parameter - depending on spatial scale of operation 
  ! sigma_max: parameter - depending on spatial scale of operation 
  !
  ! Reference: 
  ! van den Hurk, B. J. J. M., and P. Viterbo, 2003: The Torne-Kalix PILPS 2(e) experiment  
  ! as a test bedfor modifications to the ECMWF land surface scheme. Global Planet. Change, 38, 165?173
  !
  ! Dümenil, L., and E. Todini, 1992: A rainfall-runoff scheme for use in the Hamburg climate model. 
  ! Advances in Theoretical Hydrology: A Tribute to James Dooge. J. P. O?Kane, Ed., 
  ! European Geophysical Society Series on Hydrological Sciences, Vol. 1, Elsevier, 129-157. 

  implicit none
  ! input
  type(namevar),            optional,intent(in)   :: sdata(:)           ! input(optional): storage of soil data strucuture
  real(dp),                 optional,intent(in)   :: gammaPar(:)        ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)   :: opt                ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                      intent(out)  :: err                ! output: error id 
  character(len=strLen),             intent(out)  :: message            ! output: error message   
  integer(i4b), allocatable,optional,intent(out)  :: ixDepend(:)        ! output(optional): id of dependent beta parameters 
  real(dp),                 optional,intent(out)  :: binfilt_out(:,:)   ! input(optional): computed residual moisture [mm]
  ! local 
  integer(i4b)                                    :: tfopt              ! option for transfer function form used
  integer(i4b),parameter                          :: nDepend=0          ! binfilt parameter depends on no beta parameters 
  integer(i4b)                                    :: n1               ! number of 1st dimension 
  integer(i4b)                                    :: n2               ! number of 1st dimension 
  real(dp),allocatable                            :: elestd2d(:,:)
  real(dp),     parameter                         :: infilt_min=0.03_dp
  real(dp),     parameter                         :: infilt_max=0.50_dp
  
  message="infilt/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(sdata) .and. present(gammaPar) .and. present(binfilt_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    associate(g1=>gammaPar(ixGamma%binfilt1gamma1), &
              g2=>gammaPar(ixGamma%binfilt1gamma2), &
              elestd_in => sdata(ixVarSoilData%ele_std)%dvar1)
    n1=size(sdata(ixVarSoilData%hslyrs)%dvar2,1)
    n2=size(sdata(ixVarSoilData%hslyrs)%dvar2,2)
    allocate(elestd2d(n1,n2)) 
    elestd2d=spread(elestd_in,1,n1)
    select case(tfopt)
      case(1); 
        where ( elestd2d /= dmiss ) 
          binfilt_out = (log(elestd2d+verySmall)-g1)/(log(elestd2d+verySmall)+g2*10) !transfer function 
        else where
          binfilt_out = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
    !cap value with upper and lower bounds 
    where ( binfilt_out > infilt_max ) binfilt_out=infilt_max 
    where ( binfilt_out /= dmiss .and. binfilt_out < infilt_min ) binfilt_out=infilt_min 
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! residual_moist parameter 
! *********************************************************************
subroutine residMoist(err, message, ixDepend, gammaPar, residMoist_out, opt)
  implicit none
  ! input
  real(dp),                 optional,intent(in)   :: gammaPar(:)         ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)   :: opt                 ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                      intent(out)  :: err                 ! output: error id 
  character(len=strLen),             intent(out)  :: message             ! output: error message   
  integer(i4b), allocatable,optional,intent(out)  :: ixDepend(:)         ! output(optional): id of dependent beta parameters 
  real(dp),                 optional,intent(out)  :: residMoist_out(:,:) ! input(optional): computed residual moisture [mm]
  ! local 
  integer(i4b)                                    :: tfopt               ! option for transfer function form used
  integer(i4b),        parameter                  :: nDepend=0           ! residMoist parameter depends on no beta parameters 
 
  err=0;message="residMoist/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(gammaPar) .and. present(residMoist_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    select case(tfopt)
      case(1); 
        residMoist_out = 0._dp
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine
         
! ***********
!  Nijssen basiflow D1 parameter
! *********************************************************************
subroutine D1( err, message, ixDepend, sdata, ks_in, phi_in, gammaPar, D1_out, opt )
  implicit none
  ! input
  type(namevar),            optional,intent(in)   :: sdata(:)         ! input(optional): storage of soil data strucuture
  real(dp),                 optional,intent(in)   :: ks_in(:,:)       ! input(optional): porosity [fraction]  
  real(dp),                 optional,intent(in)   :: phi_in(:,:)      ! input(optional): porosity [fraction]  
  real(dp),                 optional,intent(in)   :: gammaPar(:)      ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)   :: opt              ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                      intent(out)  :: err              ! output: error id 
  character(len=strLen),             intent(out)  :: message          ! output: error message   
  integer(i4b), allocatable,optional,intent(out)  :: ixDepend(:)      ! output(optional): id of dependent beta parameters 
  real(dp),                 optional,intent(out)  :: D1_out(:,:)      ! output(optional): computed D1 parameter [day-1]
  ! local 
  integer(i4b)                                    :: tfopt            ! option for transfer function form used
  integer(i4b),        parameter                  :: nDepend=2        ! D1 parameter depends on no beta parameters 
  real(dp),allocatable                            :: slope2d(:,:)
  integer(i4b)                                    :: n1               ! number of 1st dimension 
  integer(i4b)                                    :: n2               ! number of 1st dimension 
  real(dp), allocatable                           :: S(:,:)           ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp), parameter                             :: D1_min=0.0001_dp
  real(dp), parameter                             :: D1_max=1.0_dp
   
  err=0;message="D1/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%ks,ixBeta%phi/)
  elseif ( present(sdata) .and. present(phi_in) .and. present(ks_in) .and. present(gammaPar) .and. present(D1_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! local variable allocation
    associate(g1       => gammaPar(ixGamma%D11gamma1), &
              h_in     => sdata(ixVarSoilData%hslyrs)%dvar2, &
              slope_in => sdata(ixVarSoilData%slp_mean)%dvar1 )
    n1=size(D1_out,1)
    n2=size(D1_out,2)
    allocate(S(n1,n2)) 
    allocate(slope2d(n1,n2)) 
    S=phi_in*h_in
    S=1.0_dp
    slope2d=spread(slope_in,1,n1)
    ! compute parameters 
    select case(tfopt)
      case(1); 
        where ( slope2d /= dmiss .and. ks_in /= dmiss )
          D1_out = S**(-1)*10**(-1*g1)*ks_in*(60*60*24)*(slope2d*0.01)
        else where
          D1_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate 
    ! cap value with upper and lower bounds 
    where ( D1_out > D1_max ) D1_out=D1_max 
    where ( D1_out > 0._dp .and. D1_out < D1_min ) D1_out=D1_min 
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
!  Arno basiflow Ds parameter
! *********************************************************************
subroutine Ds( err, message, ixDepend, D1_in, D3_in, Dsmax_in, gammaPar, Ds_out, opt )
  implicit none
  ! input
  real(dp),                 optional,intent(in)  :: D1_in(:,:)        ! input(optional): Nijssen baseflow D1 parameter [day^-1]
  real(dp),                 optional,intent(in)  :: D3_in(:,:)        ! input(optional): Nijssen baseflow D3 parameter 
  real(dp),                 optional,intent(in)  :: Dsmax_in(:,:)     ! input(optional): Dsmax parameter [mm]
  real(dp),                 optional,intent(in)  :: gammaPar(:)       ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)  :: opt               ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err                ! output: error id 
  character(len=strLen),            intent(out) :: message            ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)        ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: Ds_out(:,:)        ! output(optional) Ds parameter [day^-D4]
  ! local 
  integer(i4b)                                  :: tfopt              ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=3          ! Ds parameter depends on three beta parameters (D1, D3, and Dsmax)
  real(dp),    parameter                        :: Ds_min=0.0001_dp
  real(dp),    parameter                        :: Ds_max=1.0_dp

  err=0;message="Ds/"
  if (present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%D1, ixBeta%D3, ixBeta%Dsmax/)
  elseif ( present(D1_in) .and. present(D3_in) .and. present(Dsmax_in) .and. present(gammaPar) .and. present(Ds_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    select case(tfopt)
      case(1); 
        where ( D1_in /= dmiss .and. D3_in /= dmiss .and. Dsmax_in /= dmiss )
          Ds_out = D1_in * D3_in / Dsmax_in
        else where
          Ds_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    ! cap value with upper and lower bounds 
    where ( Ds_out > Ds_max ) Ds_out=Ds_max 
    where ( Ds_out > 0._dp .and. Ds_out < Ds_min ) Ds_out=Ds_min 
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
! Nijssen baseflow D2 parameter
! *********************************************************************
subroutine D2( err, message, ixDepend, sdata, ks_in, D4_in, gammaPar, D2_out, opt ) 
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: sdata(:)    ! input(optional): storage of soil data strucuture
  real(dp),                optional,intent(in)  :: Ks_in(:,:)  ! input(optional): ksat [mm/s]
  real(dp),                optional,intent(in)  :: D4_in(:,:)  ! input(optional): VIC D4 parameter [-]
  real(dp),                optional,intent(in)  :: gammaPar(:) ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt         ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err         ! output: error id 
  character(len=strLen),            intent(out) :: message     ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:) ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: D2_out(:,:) ! output(optional): D2 parameter [day^-D4]
  ! local 
  integer(i4b)                       :: tfopt                  ! option for transfer function form used
  integer(i4b),         parameter    :: nDepend=2              ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp),allocatable               :: slope2d(:,:)
  integer(i4b)                       :: n1                     ! number of element for 1st dimension
  integer(i4b)                       :: n2                     ! number of element for 2nd dimension
  real(dp), allocatable              :: S(:,:)                 ! length scaling term [mm]: 1, Max. soil storage etc
  real(dp),             parameter    :: D2_min=0.0001_dp
  real(dp),             parameter    :: D2_max=1.0_dp
  
  err=0; message="D2/"
  if (present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%ks,ixBeta%D4/)
  elseif ( present(sdata) .and. present(ks_in) .and. present(D4_in) .and. present(gammaPar) .and. present(D2_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    associate(g1=>gammaPar(ixGamma%D21gamma1), &
              slope_in => sdata(ixVarSoilData%slp_mean)%dvar1 )
    n1=size(D2_out,1)
    n2=size(D2_out,2)
    allocate(S(n1,n2))
    allocate(slope2d(n1,n2)) 
    S=1.0_dp
    slope2d=spread(slope_in,1,n1)
    select case(tfopt)
      case(1); 
         where ( slope2d /= dmiss .and. Ks_in /= dmiss )
           D2_out = S**(-1*D4_in)*10**(-1*g1)*Ks_in*(60*60*24)*(slope2d*0.01)
         else where
           D2_out = dmiss
         end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
    ! cap value with upper and lower bounds 
    where ( D2_out > D2_max ) D2_out=D2_max
    where ( D2_out > 0._dp .and. D2_out < D2_min ) D2_out=D2_min
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! ***********
! Arno baseflow Dsmax parameter
! *********************************************************************
subroutine Dsmax( err, message, ixDepend, sdata, D1_in, D2_in, D3_in, c_in, phi_in, gammaPar, Dsmax_out ,opt)          
  implicit none
  ! input
  type(namevar),            optional,intent(in)  :: sdata(:)          ! input(optional): storage of soil data strucuture
  real(dp),                 optional,intent(in)  :: D1_in(:,:)        ! input(optional): Nijssen baseflow D1 parameter [day^-1]
  real(dp),                 optional,intent(in)  :: D2_in(:,:)        ! input(optional): Nijssen baseflow D2 parameter [day^-1 mm^-(c-1)]
  real(dp),                 optional,intent(in)  :: D3_in(:,:)        ! input(optional): Nijssen baseflow D3 parameter 
  real(dp),                 optional,intent(in)  :: c_in(:,:)         ! input(optional): c parameter [mm]
  real(dp),                 optional,intent(in)  :: phi_in(:,:)       ! input(optional): porosity [fraction]  
  real(dp),                 optional,intent(in)  :: gammaPar(:)       ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)  :: opt               ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                      intent(out) :: err               ! output: error id 
  character(len=strLen),             intent(out) :: message           ! output: error message   
  integer(i4b),allocatable, optional,intent(out) :: ixDepend(:)       ! output(optional): id of dependent beta parameters 
  real(dp),                 optional,intent(out) :: Dsmax_out(:,:)    ! output(optional): [mm]
  ! local 
  integer(i4b)                                   :: tfopt             ! option for transfer function form used
  integer(i4b),parameter                         :: nDepend=5         ! Dsmax parameter depends on five beta parameters (D1,D2,D3,cexpt,and phi)
  real(dp),    parameter                         :: Dsmax_min=0.1_dp
  real(dp),    parameter                         :: Dsmax_max=30.0_dp

  err=0;message="Dsmax/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%D1, ixBeta%D2, ixBeta%D3, ixBeta%c, ixBeta%phi/)
  elseif ( present(sdata) .and. present(D1_in) .and. present(D2_in) .and. present(D3_in) .and. present(c_in) .and. present(phi_in) .and. present(gammaPar) .and. present(Dsmax_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    associate(h_in => sdata(ixVarSoilData%hslyrs)%dvar2 )
    select case(tfopt)
      case(1); 
        where ( phi_in /= dmiss .and. h_in /= dmiss )
          Dsmax_out = D2_in*(phi_in*h_in*1000-D3_in)**c_in+D1_in*(phi_in*h_in*1000)
        else where
          Dsmax_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    ! cap value with upper and lower bounds 
    where ( Dsmax_out > Dsmax_max ) Dsmax_out=Dsmax_max
    where ( Dsmax_out > 0._dp .and. Dsmax_out < Dsmax_min ) Dsmax_out=Dsmax_min
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
!  Nijssen baseflow D3 parameter
! *********************************************************************
subroutine D3( err, message, ixDepend, sdata, fc_in, gammaPar, D3_out, opt ) 
  implicit none
  ! input
  type(namevar),            optional,intent(in)  :: sdata(:)         ! input(optional): storage of soil data strucuture
  real(dp),                 optional,intent(in)  :: fc_in(:,:)       ! input(optional): field capacity [frac]
  real(dp),                 optional,intent(in)  :: gammaPar(:)      ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)  :: opt              ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                      intent(out) :: err              ! output: error id 
  character(len=strLen),             intent(out) :: message          ! output: error message   
  integer(i4b),allocatable, optional,intent(out) :: ixDepend(:)      ! output(optional): id of dependent beta parameters 
  real(dp),                 optional,intent(out) :: D3_out(:,:)      ! output(optional): [mm]
  ! local 
  integer(i4b)                                   :: tfopt              ! option for transfer function form used
  integer(i4b),         parameter                :: nDepend=1         ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp), parameter                            :: D3_min=0.0001_dp
  real(dp), parameter                            :: D3_max=1000.0_dp

  err=0; message="D3/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%fc/)
  elseif ( present(sdata) .and. present(fc_in) .and. present(gammaPar) .and. present(D3_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    associate(g1=>gammaPar(ixGamma%D31gamma1), &
              h_in => sdata(ixVarSoilData%hslyrs)%dvar2 )
    select case(tfopt)
      case(1); 
        where ( fc_in /= dmiss .and. h_in /= dmiss ) 
          D3_out = g1* fc_in * (h_in*1000)
        else where
          D3_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
 !   cap value with upper and lower bounds 
    where ( D3_out > D3_max ) D3_out=D3_max
    where ( D3_out > 0._dp .and. D3_out < D3_min ) D3_out=D3_min 
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
!  Arno baseflow Ws parameter (conversion equation)
! *********************************************************************
subroutine Ws( err, message, ixDepend, sdata, D3_in, phi_in, gammaPar, Ws_out, opt)          
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: sdata(:)      ! input(optional): storage of soil data strucuture
  real(dp),                optional,intent(in)  :: D3_in(:,:)    ! input(optional): 
  real(dp),                optional,intent(in)  :: phi_in(:,:)   ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: Ws_out(:,:)   ! output(optional): 
  ! local 
  integer(i4b)                                  :: tfopt         ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2     ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp),    parameter                        :: Ws_min=0.05_dp
  real(dp),    parameter                        :: Ws_max=1.0_dp

  err=0;message="Ws/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%D3, ixBeta%phi/)
  elseif ( present(sdata) .and. present(D3_in) .and. present(phi_in) .and. present(gammaPar) .and. present(Ws_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(h_in => sdata(ixVarSoilData%hslyrs)%dvar2)
    select case(tfopt)
      case(1); 
        where ( phi_in /= dmiss .and. h_in /= dmiss .and. D3_in /= dmiss ) 
          Ws_out = D3_in / phi_in / (h_in*1000)
        else where
          Ws_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    ! cap value with upper and lower bounds 
    where ( Ws_out > Ws_max ) Ws_out=Ws_max
    where ( Ws_out > 0._dp .and. Ws_out < Ws_min ) Ws_out=Ws_min 
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
!  Nijssen baseflow D4 parameter
! *********************************************************************
subroutine D4( err, message, ixDepend, gammaPar, D4_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: D4_out(:,:)
  ! local 
  integer(i4b)                                  :: tfopt         ! option for transfer function form used
  integer(i4b),         parameter               :: nDepend=0     ! D2 parameter depends on two beta parameters (ks and D4)
  
  err=0;message="D4/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(gammaPar) .and. present(D4_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%D41gamma1))
    select case(tfopt)
      case(1); 
        D4_out = g1 
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
!  c parameter
! *********************************************************************
subroutine cexpt( err, message, ixDepend, D4_in, gammaPar, cexpt_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: D4_in(:,:)      ! input(optional):
  real(dp),                optional,intent(in)  :: gammaPar(:)     ! input(optional):  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt             ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err             ! output: error id 
  character(len=strLen),            intent(out) :: message         ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)     ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: cexpt_out(:,:)  ! output(optional): cexpt parameter 
  ! local 
  integer(i4b)                                  :: tfopt         ! option for transfer function form used
  integer(i4b),         parameter               :: nDepend=1     ! D2 parameter depends on two beta parameters (ks and D4)
  
  err=0;message="cexpt/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%D4/)
  elseif ( present(D4_in) .and. present(gammaPar) .and. present(cexpt_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    select case(tfopt)
      case(1); 
        cexpt_out = D4_in
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
! computing expt parameter 
! *********************************************************************
subroutine expt( err, message, ixDepend, b_in, gammaPar, expt_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: b_in(:,:)      ! input: slope of cambell retention curve 
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: expt_out(:,:) ! output(optional): exponent in campbel equation [-]
  ! local 
  integer(i4b)                                  :: tfopt         ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=1     ! rexp parameter depends on one beta parameters (wp)

  err=0;message="expt/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%b/)
  elseif ( present(b_in) .and. present(gammaPar) .and. present(expt_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%exp1gamma1), &
              g2=>gammaPar(ixGamma%exp1gamma2))
    select case(tfopt)
      case(1); 
        where ( b_in /= dmiss ) 
          expt_out = g1+g2*b_in
        else where
          expt_out = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ************
! computing init_moist parameter  
! *********************************************************************
subroutine initMoist( err, message, ixDepend, sdata, phi_in, gammaPar, initMoist_out, opt )
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: sdata(:)           ! input(optional): storage of soil data strucuture
  real(dp),                optional,intent(in)  :: phi_in(:,:)        ! input: porosity [fraction]  
  real(dp),                optional,intent(in)  :: gammaPar(:)        ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt                ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err                ! output: error id 
  character(len=strLen),            intent(out) :: message            ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)        ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: initMoist_out(:,:)
  ! local  
  integer(i4b)                                  :: tfopt              ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=1          ! rexp parameter depends on one beta parameters (wp)
  
  err=0;message="initMoist/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi/)
  elseif ( present(sdata) .and. present(phi_in) .and. present(gammaPar) .and. present(initMoist_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(h_in => sdata(ixVarSoilData%hslyrs)%dvar2)
    select case(tfopt)
      case(1); 
        where ( phi_in /= dmiss ) 
          initMoist_out = phi_in*(h_in*1000.0_dp)
        else where
          initMoist_out = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
! bubble parameter 
! *********************************************************************
subroutine bubble( err, message, ixDepend, expt_in, gammaPar, bubble_out, opt)
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: expt_in(:,:)   ! input(optiona):  expt parameter
  real(dp),                optional,intent(in)  :: gammaPar(:)    ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt            ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err             ! output: error id 
  character(len=strLen),            intent(out) :: message         ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)     ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: bubble_out(:,:) ! output(optional): bubbling pressure parameter [h]
  ! local 
  integer(i4b)                                  :: tfopt           ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=1       ! rexp parameter depends on one beta parameters (wp)
  
  err=0;message="bubble/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%expt/)
  elseif ( present(expt_in) .and. present(gammaPar) .and. present(bubble_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%bbl1gamma1), &
              g2=>gammaPar(ixGamma%bbl1gamma2))
    select case(tfopt)
      case(1); 
        where ( expt_in /= dmiss ) 
          bubble_out = g1*expt_in+g2
        else where
          bubble_out = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
! soil_density parameter 
! *********************************************************************
subroutine sd( err, message, ixDepend, gammaPar, sd_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: gammaPar(:)          ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt                  ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err                  ! output: error id 
  character(len=strLen),            intent(out) :: message              ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)          ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: sd_out(:,:)
  ! local 
  real(dp)                                      :: srho=2685_dp         ! mineral density kg/cm^3
  integer(i4b)                                  :: tfopt                ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=0            ! rexp parameter depends on one beta parameters (wp)
  
  err=0;message="sd/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(gammaPar) .and. present(sd_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%sd1gamma1))
    select case(tfopt)
      case(1); 
        sd_out = g1*srho
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ***********
! WcrFrac parameter  
! *********************************************************************
subroutine WcrFrac( err, message, ixDepend, fc_in, phi_in, gammaPar, WcrFrac_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: fc_in(:,:)       ! input(optional): field capacity [fraction]  
  real(dp),                optional,intent(in)  :: phi_in(:,:)      ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)  :: gammaPar(:)      ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt              ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err              ! output: error id 
  character(len=strLen),            intent(out) :: message          ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)      ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: WcrFrac_out(:,:)
  ! local 
  integer(i4b)                                  :: tfopt            ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2        ! rexp parameter depends on one beta parameters (wp)
  
  err=0;message="WcrFrac/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi, ixBeta%fc/)
  elseif ( present(phi_in) .and. present(fc_in) .and. present(gammaPar) .and. present(WcrFrac_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%WcrFrac1gamma1))
    select case(tfopt)
      case(1); 
        where ( fc_in /= dmiss .and. phi_in /= dmiss ) 
          WcrFrac_out = g1*fc_in/phi_in
        else where
          wcrFrac_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ************
! WpwpFrac parameter  
! *********************************************************************
subroutine WpwpFrac( err, message, ixDepend, wp_in, phi_in, gammaPar, WpwpFrac_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: wp_in(:,:)        ! input(optional): wilting point [fraction]  
  real(dp),                optional,intent(in)  :: phi_in(:,:)       ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)  :: gammaPar(:)       ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt               ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err               ! output: error id 
  character(len=strLen),            intent(out) :: message           ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)       ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: WpwpFrac_out(:,:) ! output(optional): [frac] 
  ! local 
  integer(i4b)                                  :: tfopt             ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2         ! rexp parameter depends on one beta parameters (wp)
  
  err=0;message="WpwFrac/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi, ixBeta%wp/)
  elseif ( present(phi_in) .and. present(wp_in) .and. present(gammaPar) .and. present(WpwpFrac_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%WpwpFrac1gamma1))
    select case(tfopt)
      case(1); 
        where ( wp_in /= dmiss .and. phi_in /= dmiss ) 
          WpwpFrac_out = g1*wp_in/phi_in
        else where
          wpwpFrac_out = dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! ************
! TWM parameter (UZTWM for top layer and LZTWM for bottom layer)
! *********************************************************************
subroutine twm( err, message, ixDepend, sdata, fc_in, wp_in, gammaPar, twm_out, opt )
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: sdata(:)      ! input(optional): storage of soil data strucuture
  real(dp),                optional,intent(in)  :: fc_in(:,:)    ! input(optional): field capacity [-]           
  real(dp),                optional,intent(in)  :: wp_in(:,:)    ! input(optional): Wilting point
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: twm_out(:,:)  ! output(optional): tension water maximum [mm]
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2      ! twm parameter depends on one beta parameters (wp)

  err=0;message="twm/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%fc, ixBeta%wp/)
  elseif ( present(sdata) .and. present(fc_in) .and. present(wp_in) .and. present(gammaPar) .and. present(twm_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(h_in => sdata(ixVarSoilData%hslyrs)%dvar2)
    select case(tfopt)
      case(1); 
        where ( wp_in/=dmiss .and. fc_in/=dmiss ) 
          twm_out=(fc_in-wp_in)*h_in*1000.0_dp  ! convert m to mm 
        else where
          twm_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! FWM parameter (UZFWM for top layer and LZFWM for bottom layer)
! *********************************************************************
subroutine fwm( err, message, ixDepend, sdata, phi_in, fc_in, gammaPar, fwm_out, opt )
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: sdata(:)      ! input(optional): storage of soil data strucuture
  real(dp),                optional,intent(in)  :: phi_in(:,:)   ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)  :: fc_in(:,:)    ! input(optional): field capacity [-]           
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  !output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: fwm_out(:,:)  ! output(optional): free water maximum [mm]
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2      ! fwm parameter depends on one beta parameters (wp)

  err=0;message="fwm/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi, ixBeta%fc/)
  elseif ( present(phi_in) .and. present(fc_in) .and. present(gammaPar) .and. present(fwm_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(h_in => sdata(ixVarSoilData%hslyrs)%dvar2)
    select case(tfopt)
      case(1); 
        where ( phi_in/=dmiss .and. fc_in/=dmiss .and. h_in/=dmiss ) 
          fwm_out=(phi_in-fc_in)*h_in*1000.0_dp ! convert m to mm
        else where
          fwm_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! FSM parameter (LZFSM for bottom layer)
! *********************************************************************
subroutine fsm( err, message, ixDepend, fwm_in, phi_in, wp_in, gammaPar, fsm_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: fwm_in(:,:)   ! input(optional): Free water maximum 
  real(dp),                optional,intent(in)  :: phi_in(:,:)   ! input(optional): input: porosity [fraction]  
  real(dp),                optional,intent(in)  :: wp_in(:,:)    ! input(optional): Wilting point
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: fsm_out(:,:)  ! output(optional): supplementary free water maximum [mm]
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=3      ! fwm parameter depends on one beta parameters (wp)

  err=0;message="fsm/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%fwm, ixBeta%phi, ixBeta%wp/)
  elseif ( present(fwm_in) .and.  present(phi_in) .and. present(wp_in) .and. present(gammaPar) .and. present(fsm_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%fsm1gamma1))
    select case(tfopt)
      case(1); 
        where ( phi_in/=dmiss .and. wp_in/=dmiss .and. fwm_in/=dmiss ) 
          fsm_out=fwm_in*(wp_in/phi_in)**g1
        else where
          fsm_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! FPM parameter (LZFPM for bottom layer)
! *********************************************************************
subroutine fpm( err, message, ixDepend, fwm_in, fsm_in, gammaPar, fpm_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: fwm_in(:,:)   ! input(optional:  Free water maximum [mm]
  real(dp),                optional,intent(in)  :: fsm_in(:,:)   ! input(optional): supplementary Free water maximum [mm]
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: fpm_out(:,:)  ! output(optional): fpm - primary free water maximum [mm]
  ! local  
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2      ! fpm parameter depends on one beta parameters (wp)

  err=0;message="fpm/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%fwm, ixBeta%fsm/)
  elseif ( present(fwm_in) .and. present(fsm_in) .and. present(gammaPar) .and. present(fpm_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    select case(tfopt)
      case(1); 
        where ( fwm_in/=dmiss .and. fsm_in/=dmiss )
          fpm_out=fwm_in-fsm_in
        else where
          fpm_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! ZK parameter (UZK for top layer)
! *********************************************************************
subroutine zk( err, message, ixDepend, phi_in, fc_in, gammaPar, zk_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: phi_in(:,:)   ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)  :: fc_in(:,:)    ! input(optional): field capacity [frac]
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  !output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: zk_out(:,:)   ! output(optional): computed zk - draw coefficient from free water content [/day] 
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2      ! zk parameter depends on one beta parameters (wp)

  err=0;message="zk/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi, ixBeta%fc/)
  elseif ( present(phi_in) .and. present(fc_in) .and. present(gammaPar) .and. present(zk_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%zk1gamma1))
    select case(tfopt)
      case(1); 
        where ( phi_in/=dmiss .and. fc_in/=dmiss ) 
          zk_out=1.0_dp-(fc_in/phi_in)**g1
        else where
          zk_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! ZSK parameter (LZSK for top layer)
! *********************************************************************
subroutine zsk( err, message, ixDepend,  phi_in, fc_in, wp_in, gammaPar, zsk_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: phi_in(:,:)   ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)  :: fc_in(:,:)    ! input(optional): field capacity [frac]
  real(dp),                optional,intent(in)  :: wp_in(:,:)    ! input(optional): Wilting point [frac]
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: zsk_out(:,:)  ! output(optional): computed zsk draw coefficient from supplementary free water content [/day]
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=3      ! zsk parameter depends on one beta parameters (wp)

  err=0;message="zsk/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi, ixBeta%wp, ixBeta%fc/)
  elseif ( present(phi_in) .and. present(wp_in) .and. present(fc_in) .and. present(gammaPar) .and. present(zsk_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%zsk1gamma1),&
              g2=>gammaPar(ixGamma%zsk1gamma2))
    select case(tfopt)
      case(1); 
        where ( phi_in/=dmiss .and. fc_in/=dmiss .and. wp_in/=dmiss ) 
          zsk_out=(1.0_dp-(fc_in/phi_in)**g1)/(1.0_dp+g2*(1.0_dp-wp_in)) 
        else where
          zsk_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! ZPK parameter (LZPK for top layer)
! *********************************************************************
subroutine zpk( err, message, ixDepend, sdata, ks_in, myu_in, gammaPar, zpk_out, opt )
  implicit none
  ! input
  type(namevar),            optional,intent(in) :: sdata(:)      ! input(optional): storage of soil data strucuture
  real(dp),                 optional,intent(in) :: ks_in(:,:)    ! input(optional): porosity [fraction]  
  real(dp),                 optional,intent(in) :: myu_in(:,:)   ! input(optional): specific yield [-]
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: zpk_out(:,:)
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2      ! zpk parameter depends on one beta parameters (wp)
  real(dp)                                      :: dt             ! time step of the simulation [hr]

  err=0;message="zpk/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%ks, ixBeta%myu/)
  elseif ( present(sdata) .and. present(ks_in) .and. present(myu_in) .and. present(gammaPar) .and. present(zpk_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    dt=24.0_dp ! unit: hr
    associate(g1=>gammaPar(ixGamma%zpk1gamma1),&
              h_in     => sdata(ixVarSoilData%hslyrs)%dvar2)
    select case(tfopt)
      case(1); 
        where ( ks_in/=dmiss .and. myu_in/=dmiss .and. h_in/= dmiss ) 
          zpk_out=1.0_dp-exp(-1.0_dp*g1**2.0_dp*pi*ks_in*60.0_dp*h_in*1000.0_dp*dt/myu_in)
        else where
          zpk_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

! *********************************************************************
! PFREE parameter 
! *********************************************************************
subroutine pfree( err, message, ixDepend, phi_in, wp_in, gammaPar, pfree_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: phi_in(:,:)    ! input(optional): porosity [frac] 
  real(dp),                optional,intent(in)  :: wp_in(:,:)     ! input(optional): wilting point [fract] 
  real(dp),                optional,intent(in)  :: gammaPar(:)    ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt            ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err            ! output: error id 
  character(len=strLen),            intent(out) :: message        ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)    ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: pfree_out(:,:) ! output(optional): ratio of max and min percolation rates [day-1]
  ! local 
  integer(i4b)                                  :: tfopt          ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=2      ! pfree parameter depends on one beta parameters (wp)

  err=0;message="pfree/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%phi, ixBeta%wp/)
  elseif ( present(phi_in) .and. present(wp_in) .and. present(gammaPar) .and. present(pfree_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt)) tfopt=opt
    associate(g1=>gammaPar(ixGamma%pfree1gamma1))
    select case(tfopt)
      case(1); 
        where ( phi_in/=dmiss .and. wp_in/=dmiss ) 
          pfree_out=(wp_in/phi_in)**g1
        else where 
          pfree_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! ZPERC parameter 
! *********************************************************************
subroutine zperc( err, message, ixDepend, twm_in, fsm_in, zsk_in, fpm_in, zpk_in, gammaPar, zperc_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: twm_in(:,:)   ! input: tension water maximum [mm] 
  real(dp),                optional,intent(in)  :: fsm_in(:,:)   ! input: supplemental free water maximum [mm] 
  real(dp),                optional,intent(in)  :: zsk_in(:,:)   ! input: flow rate from supplemental free water maximum [day-1] 
  real(dp),                optional,intent(in)  :: fpm_in(:,:)   ! input: primary free water maximum [mm] 
  real(dp),                optional,intent(in)  :: zpk_in(:,:)   ! input: flow rate from primary free water maximum [day-1] 
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! Optional output: id of dependent beta parameters 
  real(dp),                optional,intent(out) :: zperc_out(:,:)    ! output: ratio of max and min percolation rates [day-1]
  ! local 
  integer(i4b)                                  :: tfopt         ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=5     ! zperc parameter depends on one beta parameters (wp)
                       
  err=0;message="zperc/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%twm, ixBeta%fsm, ixBeta%zsk, ixBeta%fpm, ixBeta%zpk/)
  elseif ( present(twm_in) .and. present(fsm_in) .and. present(zsk_in) .and. present(fpm_in) .and. present(zpk_in) .and. present(gammaPar) .and. present(zperc_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    select case(tfopt)
      case(1); 
        where ( twm_in/=dmiss .and. fsm_in/=dmiss .and. zsk_in/=dmiss .and. fpm_in/=dmiss .and. zpk_in/=dmiss ) 
          zperc_out=(twm_in+fsm_in*(1.0_dp-zsk_in)+fsm_in*(1.0_dp-zpk_in))/ &
                (fsm_in*zsk_in+fpm_in*zpk_in)
        else where
          zperc_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! REXP parameter 
! *********************************************************************
subroutine rexp( err, message, ixDepend, wp_in, gammaPar, rexp_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)  :: wp_in(:,:)    ! input(optional): willting point [fraction]  
  real(dp),                optional,intent(in)  :: gammaPar(:)   ! input(optional:  gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt           ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err           ! output: error id 
  character(len=strLen),            intent(out) :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)   ! Optional output: id of dependent beta parameters 
  real(dp),                optional,intent(out) :: rexp_out(:,:) ! output: tension water maximum [mm]
  ! local 
  integer(i4b)                                  :: tfopt         ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=1     ! rexp parameter depends on one beta parameters (wp)
  
  err=0;message="rexp/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%wp/)
  elseif ( present(wp_in) .and. present(gammaPar) .and. present(rexp_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    associate(g1=>gammaPar(ixGamma%rexp1gamma1))
    select case(tfopt)
      case(1); 
        where ( wp_in/=dmiss ) 
          rexp_out=sqrt(wp_in/(g1-0.001))
        else where
          rexp_out=dmiss
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for saturated hydraulic conductivity (ks)
! *********************************************************************
subroutine ks( err, message, ixDepend, sdata, gammaPar, ks_out, opt )
  implicit none
  ! input
  type(namevar),optional,intent(in) :: sdata(:)         ! storage of soil data strucuture
  real(dp),    optional,intent(in)  :: gammaPar(:)          ! input: gamma parameter array 
  integer(i4b),optional,intent(in)  :: opt                  ! input: option for transfer function form
  ! output
  integer(i4b),         intent(out) :: err                ! output: error id 
  character(len=strLen),intent(out) :: message            ! output: error message   
  integer(i4b), allocatable, optional,intent(out)  :: ixDepend(:)      ! Optional output: id of dependent beta parameters 
  real(dp),    optional,intent(out) :: ks_out(:,:)           ! output: mm/s
  ! local 
  integer(i4b)                      :: tfopt              ! option for transfer function form used
  integer(i4b),         parameter   :: nDepend=0          ! D2 parameter depends on two beta parameters (ks and D4)

  err=0;message="ks/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(sdata) .and. present(gammaPar) .and. present(ks_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Cosby et al. WRR 1984
    ! opt 2: campbell & shiozawa 1994 
    associate(g1=>gammaPar(ixGamma%ks1gamma1), &
              g2=>gammaPar(ixGamma%ks1gamma2), &
              g3=>gammaPar(ixGamma%ks1gamma3), &
              g4=>gammaPar(ixGamma%ks2gamma1), &
              g5=>gammaPar(ixGamma%ks2gamma2), &
              g6=>gammaPar(ixGamma%ks2gamma3), & 
              sand_in => sdata(ixVarSoilData%sand_frc)%dvar2, & 
              clay_in => sdata(ixVarSoilData%clay_frc)%dvar2 )
    select case(tfopt)
      case(1); 
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          ks_out = g1 + g2*sand_in + g3*clay_in
          ks_out = (10**ks_out)*25.4/60/60   ! 25.4 mm/inch. Cosby give Ksat in inch/hr 
        else where
          ks_out = dmiss 
        end where
      case(2); 
        where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
          ks_out = g4*exp(g5*sand_in+g6*clay_in)
        else where
          ks_out = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized'; stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for bulk density 
! *********************************************************************
subroutine bd( err, message, ixDepend, sdata, gammaPar, bd_out, opt )
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: sdata(:)         ! input(optional): storage of soil data strucuture
  real(dp),                optional,intent(in)  :: gammaPar(:)      ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt              ! input(optional): option for transfer function form
  ! output 
  integer(i4b),                     intent(out) :: err              ! output: error id 
  character(len=strLen),            intent(out) :: message          ! output: error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)      ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: bd_out(:,:)      ! output(optional): computed bd
  ! local 
  integer(i4b)                                  :: tfopt            ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=0        ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp),    parameter                        :: bd_min=805.0_dp  !
  real(dp),    parameter                        :: bd_max=1880.0_dp !
  real(dp),allocatable                          :: bdslope(:,:)
  real(dp),allocatable                          :: bd_temp(:,:)
  integer(i4b)                                  :: n1               ! number of 1st dimension 
  integer(i4b)                                  :: n2               ! number of 2nd dimension 

  err=0;message="bd/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(sdata) .and. present(gammaPar) .and. present(bd_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    associate ( g1    => gammaPar(ixGamma%bd1gamma1), &
                bd_in => sdata(ixVarSoilData%bulk_density)%dvar2 )
    n1=size(bd_in,1)
    n2=size(bd_in,2)
    allocate(bdslope(n1,n2))
    allocate(bd_temp(n1,n2))
    bdslope=0.0_dp
    bd_temp=0.0_dp
    select case(tfopt)
      case(1);  ! 
        where ( bd_in /= dmiss ) 
          bd_temp = g1*bd_in
          bdslope=(bd_temp-bd_min)/(bd_max-bd_min)
          where ( bdslope > 1.0_dp) bdslope=1.0_dp
          where ( bdslope < 0.0_dp) bdslope=0.0_dp
          bd_out = bdslope*(bd_max-bd_min)+bd_min
        else where
          bd_out = dmiss 
        end where
      case default;print*,trim(message)//'opt not recognized';stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for porosity 
! *********************************************************************
subroutine phi( err, message, ixDepend, sdata, gammaPar, phi_out, opt ) 
  implicit none
  ! input
  type(namevar),            optional,intent(in)   :: sdata(:)       ! input(optional): storage of soil data strucuture
  real(dp),                 optional,intent(in)   :: gammaPar(:)    ! input(optional): gamma parameter array 
  integer(i4b),             optional,intent(in)   :: opt            ! input(optional): id for transfer function form
  ! output 
  integer(i4b),                      intent(out)  :: err            ! output:           error id 
  character(len=strLen),             intent(out)  :: message        ! output:           error message   
  integer(i4b),allocatable, optional,intent(out)  :: ixDepend(:)    ! output(optional): id of dependent beta parameters 
  real(dp),                 optional,intent(out)  :: phi_out(:,:)   ! output(optional): estimated porosity [fraction]
  ! local 
  integer(i4b)                                    :: tfopt          ! option for transfer function form used
  integer(i4b),         parameter                 :: nDepend=0      ! D2 parameter depends on two beta parameters (ks and D4)

  err=0;message="phi/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(sdata) .and. present(gammaPar) .and. present(phi_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Cosby et al. WRR 1984
    ! opt 2: Zacharias & Wessolek 2007
    associate(g1      => gammaPar(ixGamma%phi1gamma1),          &
              g2      => gammaPar(ixGamma%phi1gamma2),          &
              g3      => gammaPar(ixGamma%phi1gamma3),          &
              g4      => gammaPar(ixGamma%phi2gamma1),          &
              g5      => gammaPar(ixGamma%phi2gamma2),          &
              g6      => gammaPar(ixGamma%phi2gamma3),          &
              g7      => gammaPar(ixGamma%phi2gamma4),          &
              g8      => gammaPar(ixGamma%phi2gamma5),          &
              g9      => gammaPar(ixGamma%phi2gamma6),          &
              sand_in => sdata(ixVarSoilData%sand_frc)%dvar2, & 
              clay_in => sdata(ixVarSoilData%clay_frc)%dvar2, & 
              bd_in   => sdata(ixVarSoilData%bulk_density)%dvar2 ) 
      select case(tfopt)
        case(1);  ! Cosby
          where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
            phi_out = ( g1+g2*sand_in+g3*clay_in )/100.0_dp
          else where
            phi_out = dmiss 
          end where
        case(2);  ! Zacharias & Wessolek 2007  
          where ( sand_in /= dmiss .and. clay_in /= dmiss .and. bd_in /= dmiss ) 
            where ( sand_in < 66.5_dp) 
              phi_out = g4+g5*clay_in+g6*bd_in/1000._dp
            else where
              phi_out = g7+g8*clay_in+g9*bd_in/1000._dp
            end where 
          else where
            phi_out = dmiss 
          end where
        case default;print*,trim(message)//'opt not recognized';stop
      end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for field capacity 
! *********************************************************************
subroutine fc( err, message, ixDepend, sdata, phi_in, psis_in, b_in, gammaPar, fc_out, opt )
  implicit none
  ! input
  type(namevar),optional,intent(in)  :: sdata(:)         ! storage of soil data strucuture
  real(dp),     optional,intent(in)  :: phi_in(:,:)    ! input: porosity [fraction]  
  real(dp),     optional,intent(in)  :: psis_in(:,:)   ! input: saturation matric potential [kPa] 
  real(dp),     optional,intent(in)  :: b_in(:,:)      ! input: slope of cambell retention curve 
  real(dp),     optional,intent(in)  :: gammaPar(:)    ! input: gamma parameter array 
  integer(i4b), optional,intent(in)  :: opt            ! id for transfer function form
  ! output 
  integer(i4b),         intent(out)  :: err            ! output: error id 
  character(len=strLen),intent(out)  :: message        ! output: error message   
  integer(i4b), allocatable, optional,intent(out)  :: ixDepend(:)      ! Optional output: id of dependent beta parameters 
  real(dp),    optional,intent(out)  :: fc_out(:,:)    ! estimated field capacity [fraction]
  ! local 
  integer(i4b)                       :: tfopt            ! option for transfer function form used
  integer(i4b),         parameter    :: nDepend=3        ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp),allocatable               :: psi_fc(:,:)      ! matric potential at field capacity [kPa]  
  integer(i4b)                       :: nSpoly           ! number of soil polygon 
  integer(i4b)                       :: nSLyr            ! number of soil layer 

  err=0;message="fc/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%psis, ixBeta%phi, ixBeta%b/)
  elseif ( present(sdata) .and. present(phi_in) .and. present(psis_in) .and. present(b_in) .and. present(gammaPar) .and. present(fc_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Campbell 1974 
    nSpoly=size(phi_in,1)
    nSLyr=size(phi_in,2)
    allocate(psi_fc(nSpoly,nSLyr))
    psi_fc(:,:)=-20
    associate( g1      => gammaPar(ixGamma%fc1gamma1), &
               sand_in => sdata(ixVarSoilData%sand_frc)%dvar2 )
    where (sand_in > 69) psi_fc=-10
    select case(tfopt)
      case(1);  !campbell
        where ( phi_in /= dmiss .and. sand_in /= dmiss .and. b_in /= dmiss .and. psis_in /= dmiss ) 
          fc_out = g1*phi_in*(psi_fc/psis_in)**(-1/b_in)
        else where
          fc_out = dmiss 
        end where
      case default; message=trim(message)//'opt not recognized' 
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for wilting point 
! *********************************************************************
subroutine wp( err, message, ixDepend, phi_in, psis_in, b_in, gammaPar, wp_out, opt ) 
  implicit none
  ! input
  real(dp),    optional,intent(in)    :: phi_in(:,:)    ! input: porosity [fraction]  
  real(dp),    optional,intent(in)    :: psis_in(:,:)   ! input: saturation matric potential [kPa]  
  real(dp),    optional,intent(in)    :: b_in(:,:)      ! input: slope of cambell retention curve [-]
  real(dp),    optional,intent(in)    :: gammaPar(:)    ! input: gamma parameter array 
  integer(i4b),optional,intent(in)    :: opt            ! input: option for transfer function form
  ! output
  integer(i4b),         intent(out)   :: err            ! output: error id 
  character(len=strLen),intent(out)   :: message        ! output: error message   
  integer(i4b), allocatable, optional,intent(out)  :: ixDepend(:)      ! Optional output: id of dependent beta parameters 
  real(dp),    optional,intent(out)   :: wp_out(:,:)    ! estimated field capacity [frac]
  ! local 
  integer(i4b)                        :: tfopt          ! option for transfer function form used
  integer(i4b),         parameter     :: nDepend=3      ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp),allocatable                :: psi_wp(:,:)    ! matric potential at wilting point [kPa]  
  integer(i4b)                        :: nSpoly         ! number of soil polygon 
  integer(i4b)                        :: nSLyr          ! number of soil layer 

  err=0;message="wp/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%ks,ixBeta%phi,ixBeta%b/)
  elseif ( present(phi_in) .and. present(psis_in) .and. present(b_in) .and. present(gammaPar) .and. present(wp_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Campbell 1974 
    nSpoly=size(phi_in,1)
    nSLyr=size(phi_in,2)
    allocate(psi_wp(nSpoly,nSLyr))
    psi_wp(:,:)=-1500
    associate(g1=>gammaPar(ixGamma%wp1gamma1))
    select case(tfopt)
      case(1);  !Cosby et al. 
        where ( phi_in /= dmiss .and. b_in /= dmiss .and. psis_in /= dmiss ) 
          wp_out = g1*phi_in*(psi_wp/psis_in)**(-1/b_in)
        else where
          wp_out = dmiss 
        end where
      case default; print*,trim(message)//'opt not recognized';stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for b (slope of retention curve in log space)
! *********************************************************************
subroutine retcurve( err, message, ixDepend, sdata, gammaPar, retcurve_out, opt )
  implicit none
  ! input
  type(namevar),optional,intent(in)   :: sdata(:)         ! storage of soil data strucuture
  real(dp),    optional,intent(in)    :: gammaPar(:)        ! input: gamma parameter array 
  integer(i4b),optional,intent(in)    :: opt                ! input: option for transfer function form
  ! output
  integer(i4b),         intent(out)   :: err                ! output: error id 
  character(len=strLen),intent(out)   :: message            ! output: error message   
  integer(i4b), allocatable, optional,intent(out)  :: ixDepend(:)      ! Optional output: id of dependent beta parameters 
  real(dp),    optional,intent(out)   :: retcurve_out(:,:)  ! computed [-] 
  ! local 
  integer(i4b)                        :: tfopt              ! option for transfer function form used
  integer(i4b),         parameter     :: nDepend=0          ! D2 parameter depends on two beta parameters (ks and D4)

  err=0;message="retcurve/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(sdata) .and. present(gammaPar) .and. present(retcurve_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Cosby et al. WRR 1984
    associate(g1      => gammaPar(ixGamma%b1gamma1), &
              g2      => gammaPar(ixGamma%b1gamma2), &
              g3      => gammaPar(ixGamma%b1gamma2), &
              sand_in => sdata(ixVarSoilData%sand_frc)%dvar2, & 
              clay_in => sdata(ixVarSoilData%clay_frc)%dvar2 )
      select case(tfopt)
        case(1); 
          where ( sand_in /= dmiss .and. clay_in /= dmiss ) 
            retcurve_out = g1+g2*sand_in+g3*clay_in
          else where
            retcurve_out = dmiss 
          end where
        case default; print*,trim(message)//'opt not recognized';stop
      end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for saturation matric potential 
! *********************************************************************
subroutine psis( err, message, ixDepend, sdata, gammaPar, psis_out, opt )
  implicit none
  ! input
  type(namevar),optional,intent(in)   :: sdata(:)         ! storage of soil data strucuture
  real(dp),     optional,intent(in)   :: gammaPar(:)    ! input: gamma parameter array 
  integer(i4b), optional,intent(in)   :: opt            ! input: option for transfer function form
  ! output
  integer(i4b),         intent(out)   :: err            ! output: error id 
  character(len=strLen),intent(out)   :: message        ! output: error message   
  integer(i4b), allocatable, optional,intent(out)  :: ixDepend(:)      ! Optional output: id of dependent beta parameters 
  real(dp),    optional,intent(out)   :: psis_out(:,:)  ! output: saturation matric potential [kPa]  
  ! local 
  integer(i4b)                        :: tfopt          ! option for transfer function form used
  integer(i4b),         parameter     :: nDepend=0      ! D2 parameter depends on two beta parameters (ks and D4)

  err=0;message="psis/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(sdata) .and. present(gammaPar) .and. present(psis_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Cosby et al. WRR 1984
    associate(g1      => gammaPar(ixGamma%psis1gamma1), &
              g2      => gammaPar(ixGamma%psis1gamma2), &
              g3      => gammaPar(ixGamma%psis1gamma2), &
              sand_in => sdata(ixVarSoilData%sand_frc)%dvar2, & 
              silt_in => sdata(ixVarSoilData%silt_frc)%dvar2 )
      select case(tfopt)
        case(1);  !Cosby et al. 
          where ( sand_in /= dmiss .and. silt_in /= dmiss ) 
            psis_out = g1 + g2*sand_in + g3*silt_in
            psis_out = -1.0_dp*10**psis_out*0.0980665_dp        ! 0.0980665 kPa/cm-H2O. Cosby give psi_sat in cm of water (cm-H2O)
          else where
            psis_out = dmiss 
          end where
        case default; print*,trim(message)//'opt not recognized';stop
      end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! pedo-transfer function for specific yield  
! *********************************************************************
subroutine myu( err, message, ixDepend, phi_in, fc_in, gammaPar, myu_out, opt )
  implicit none
  ! input
  real(dp),                optional,intent(in)   :: phi_in(:,:)   ! input(optional): porosity [fraction]  
  real(dp),                optional,intent(in)   :: fc_in(:,:)    ! input(optional): saturation matric potential [kPa]  
  real(dp),                optional,intent(in)   :: gammaPar(:)   ! input(optional): gamma parameter array 
  integer(i4b),            optional,intent(in)   :: opt           ! input(optional): option for transfer function form
  ! output
  integer(i4b),                     intent(out)  :: err           ! output: error id 
  character(len=strLen),            intent(out)  :: message       ! output: error message   
  integer(i4b),allocatable,optional,intent(out)  :: ixDepend(:)   ! output(optional): id of dependent beta parameters 
  real(dp),                optional,intent(out)  :: myu_out(:,:)  ! output(optional): specific yield [-]  
  ! local 
  integer(i4b)                                   :: tfopt         ! option for transfer function form used
  integer(i4b),parameter                         :: nDepend=2     ! myu parameter depends on two beta parameters (phi and fc)

  err=0;message="myu/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(nDepend),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=(/ixBeta%fc,ixBeta%phi/)
  elseif ( present(phi_in) .and. present(fc_in) .and. present(gammaPar) .and. present(myu_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if (present(opt) ) tfopt=opt
    ! opt 1: Koren et al. 2003
    associate(g1=>gammaPar(ixGamma%myu1gamma1), &
              g2=>gammaPar(ixGamma%myu1gamma2))
    select case(tfopt)
      case(1);  ! koren
        where ( phi_in /= dmiss .and. fc_in /= dmiss ) 
          myu_out = g1*(phi_in-fc_in)**g2 
        else where
          myu_out = dmiss 
        end where
      case default;print*,trim(message)//'opt not recognized';stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine

! *********************************************************************
! Monthly LAI 
! *********************************************************************
subroutine lai( err, message, ixDepend, vdata, gammaPar, lai_out, opt )
  implicit none
  ! input
  type(namevar),           optional,intent(in)  :: vdata(:)         ! input(option):  storage of soil data strucuture
  real(dp),                optional,intent(in)  :: gammaPar(:)      ! inputoption):   gamma parameter array 
  integer(i4b),            optional,intent(in)  :: opt              ! input(option):  option for transfer function form
  ! output
  integer(i4b),                     intent(out) :: err              ! output:         error id 
  character(len=strLen),            intent(out) :: message          ! output:         error message   
  integer(i4b),allocatable,optional,intent(out) :: ixDepend(:)      ! output(option): id of dependent beta parameters 
  real(dp),                optional,intent(out) :: lai_out(:,:)     ! output(option): saturation matric potential [kPa]  
  ! local 
  integer(i4b)                                  :: tfopt            ! option for transfer function form used
  integer(i4b),parameter                        :: nDepend=0        ! D2 parameter depends on two beta parameters (ks and D4)
  real(dp),parameter                            :: lai_min=0.0_dp   ! minimum plausible LAI (0 m2/m2 for bare ground)
  real(dp),parameter                            :: lai_max=10.0_dp  ! maximum plausible LAI (10 m2/m2 for dense conifer forest)
  real(dp),allocatable                          :: lai_temp(:,:)
  real(dp),allocatable                          :: laislope(:,:)
  integer(i4b)                                  :: n1               ! number of 1st dimension 
  integer(i4b)                                  :: n2               ! number of 2nd dimension 

  err=1;message="lai/"
  if ( present(ixDepend) ) then ! setup dependency
    allocate(ixDepend(1),stat=err); if(err/=0)then;message=trim(message)//'error allocating ixDepend';return;endif
    ixDepend=-999_i4b
  elseif ( present(vdata) .and. present(gammaPar) .and. present(lai_out) )then ! compute parameters with TF 
    tfopt=1_i4b
    if ( present(opt) ) tfopt=opt
    associate(g1=>gammaPar(ixGamma%lai1gamma1), &
              lai_in => vdata(ixVarVegData%lai)%dvar2 )
    n1=size(lai_in,1)
    n2=size(lai_in,2)
    allocate(laislope(n1,n2))
    allocate(lai_temp(n1,n2))
    laislope=0.0_dp
    lai_temp=0.0_dp
    select case(tfopt)
      case(1);  ! 
        where ( lai_in /= dmiss ) 
          lai_temp = g1*lai_in*0.1
          laislope=(lai_temp-lai_min)/(lai_max-lai_min)
          where ( laislope > 1.0_dp) laislope=1.0_dp
          where ( laislope < 0.0_dp) laislope=0.0_dp
          lai_out = laislope*(lai_max-lai_min)+lai_min
        else where
          lai_out = dmiss 
        end where
      case default;print*,trim(message)//'opt not recognized';stop
    end select
    end associate
  else
    err=10;message=trim(message)//'WrongOptionalInputs'; return 
  endif
  return
end subroutine 

end module tf 
