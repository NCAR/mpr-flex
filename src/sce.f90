module sce

  use nrtype 

  implicit none

  private

  public::sceua

contains

  subroutine sceua(obj_func, pini, prange, maxiter, kstop, pcento, seed, ngs, tmp_file, mask)
  
  use public_var
  use mo_xor4096, only: xor4096, xor4096g

  !  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
  !     -- VERSION 2.1
  !
  !  BY QINGYUN DUAN
  !  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
  !  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
  !  (602) 621-9360, EMAIL: DUAN@HWR.ARIZONA.EDU
  !
  !  WRITTEN IN OCTOBER 1990.
  !  REVISED IN AUGUST 1991 C  REVISED IN APRIL 1992
  !  RE-WRITTEN in F90 in JULY 2017 
  !
  !  STATEMENT BY AUTHOR:
  !  --------------------
  !
  !     THIS GENERAL PURPOSE GLOBAL OPTIMIZATION PROGRAM IS DEVELOPED AT
  !     THE DEPARTMENT OF HYDROLOGY & WATER RESOURCES OF THE UNIVERSITY
  !     OF ARIZONA.  FURTHER INFORMATION REGARDING THE SCE-UA METHOD CAN
  !     BE OBTAINED FROM DR. Q. DUAN, DR. S. SOROOSHIAN OR DR. V.K. GUPTA
  !     AT THE ADDRESS AND PHONE NUMBER LISTED ABOVE.  WE REQUEST ALL
  !     USERS OF THIS PROGRAM MAKE PROPER REFERENCE TO THE PAPER ENTITLED
  !     'EFFECTIVE AND EFFICIENT GLOBAL OPTIMIZATION FOR CONCEPTUAL
  !     RAINFALL-RUNOFF MODELS' BY DUAN, Q., S. SOROOSHIAN, AND V.K. GUPTA,
  !     WATER RESOURCES RESEARCH, VOL 28(4), PP.1015-1031, 1992.
  !
  !  LIST OF INPUT ARGUEMENT VARIABLES
  !
  !     obj_func    = object function
  !     pini(:)     = Initial parameter set 
  !     prange(:,:) = Min/Max values for parameter set
  !
  !     -- LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:
  !
  !     ngs = NUMBER OF COMPLEXES IN THE INITIAL POPULATION
  !
  !     -- CONVERGENCE CHECK PARAMETERS
  !
  !     maxiter = MAX NO. OF TRIALS ALLOWED BEFORE OPTIMIZATION IS TERMINATED
  !     kstop  = NUMBER OF SHUFFLING LOOPS IN WHICH THE CRITERION VALUE MUST
  !              CHANG BY THE GIVEN PERCENTAGE BEFORE OPTIMIZATION IS TERMINATED
  !     pcento = PERCENTAGE BY WHICH THE CRITERION VALUE MUST CHANGE IN
  !              GIVEN NUMBER OF SHUFFLING LOOPS
  !
  !  LIST OF LOCAL VARIABLES
  !     INIFLG = FLAG ON WHETHER TO INCLUDE THE INITIAL POINT IN POPULATION
  !         = 0, NOT INCLUDED
  !         = 1, INCLUDED
  !     IPRINT = FLAG FOR CONTROLLING PRINT-OUT AFTER EACH SHUFFLING LOOP
  !         = 0, PRINT INFORMATION ON THE BEST POINT OF THE POPULATION
  !         = 1, PRINT INFORMATION ON EVERY POINT OF THE POPULATION
  !     pnum = NUMBER OF PARAMETERS TO BE OPTIMIZED
  !     npg = NUMBER OF POINTS IN EACH COMPLEX
  !     npt = TOTAL NUMBER OF POINTS IN INITIAL POPULATION (npt=ngs*npg)
  !     nps = NUMBER OF POINTS IN A SUB-COMPLEX
  !     nspl = NUMBER OF EVOLUTION STEPS ALLOWED FOR EACH COMPLEX BEFORE
  !         COMPLEX SHUFFLING
  !     mings = MINIMUM NUMBER OF COMPLEXES REQUIRED, IF THE NUMBER OF
  !         COMPLEXES IS ALLOWED TO REDUCE AS THE OPTIMIZATION PROCEEDS
  !     seed = INITIAL RANDOM SEED
  !     IPCNVG = FLAG INDICATING WHETHER PARAMETER CONVERGENCE IS REACHED
  !         (I.E., CHECK IF GNRNG IS LESS THAN 0.001)
  !         = 0, PARAMETER CONVERGENCE NOT SATISFIED
  !         = 1, PARAMETER CONVERGENCE SATISFIED
  !     ISCE = UNIT NUMBER FOR SCE OUTPUT (MPC ADDITION)
  !     X(.,.) = COORDINATES OF POINTS IN THE POPULATION
  !     XF(.) = FUNCTION VALUES OF X(.,.)
  !     XX(.) = COORDINATES OF A SINGLE POINT IN X
  !     CX(.,.) = COORDINATES OF POINTS IN A COMPLEX
  !     CF(.) = FUNCTION VALUES OF CX(.,.)
  !     S(.,.) = COORDINATES OF POINTS IN THE CURRENT SIMPLEX
  !     SF(.) = FUNCTION VALUES OF S(.,.)
  !     BESTX(.) = BEST POINT AT CURRENT SHUFFLING LOOP
  !     BESTF = FUNCTION VALUE OF BESTX(.)
  !     WORSTX(.) = WORST POINT AT CURRENT SHUFFLING LOOP
  !     WORSTF = FUNCTION VALUE OF WORSTX(.)
  !     XNSTD(.) = STANDARD DEVIATION OF PARAMETERS IN THE POPULATION
  !     GNRNG = NORMALIZED GEOMETRIC MEAN OF PARAMETER RANGES
  !     LCS(.) = INDICES LOCATING POSITION OF S(.,.) IN X(.,.)
  !     BOUND(.) = BOUND ON ITH VARIABLE BEING OPTIMIZED
  !     ngs1 = NUMBER OF COMPLEXES IN CURRENT POPULATION
  !     ngs2 = NUMBER OF COMPLEXES IN LAST POPULATION
  !     CRITER(.) = VECTOR CONTAINING THE BEST CRITERION VALUES OF THE LAST
  
    implicit none
    ! input
    interface 
      function obj_func(pp)
        use nrtype
        implicit none
        real(dp), intent(in) :: pp(:)
        real(dp)             :: obj_func
      end function obj_func
    end interface
    real(dp),              intent(in)    :: pini(:) 
    real(dp),              intent(in)    :: prange(:,:) 
    integer(i8b),          intent(in)    :: maxiter 
    integer(i4b),          intent(in)    :: kstop 
    real(dp),              intent(in)    :: pcento 
    integer(i8b),          intent(in)    :: seed 
    integer(i4b),          intent(in)    :: ngs
    character(*),          intent(in)    :: tmp_file 
    logical,     optional, intent(in)    :: mask(:)        ! parameter to be optimized (true or false)
    ! local
    integer(i4b), parameter              :: INIFLG=1
    integer(i4b), parameter              :: IPRINT=1     
    integer(i4b), parameter              :: ISCE=999
    integer(i4b)                         :: pnum
    integer(i4b)                         :: npg 
    integer(i4b)                         :: nps 
    integer(i4b)                         :: nspl 
    integer(i4b)                         :: mings 
    integer(i4b)                         :: LCS(50)
    integer(i4b)                         :: LPOS 
    integer(i4b)                         :: IPCNVG 
    integer(i4b)                         :: ICALL
    integer(i4b)                         :: NLOOP, LOOP, IGS ! loop index
    integer(i4b)                         :: I, J, K, LDX     ! loop index
    integer(i4b)                         :: K1, K2           ! loop index
    integer(i4b)                         :: npt 
    integer(i4b)                         :: npt1
    integer(i4b)                         :: ngs1 
    integer(i4b)                         :: ngs2 
    real(dp), dimension(size(pini))      :: pval             ! initial value of parameter - either input or restart
    real(dp), dimension(2000,size(pini)) :: X
    real(dp), dimension(2000,size(pini)) :: CX
    real(dp), dimension(50,size(pini))   :: S
    real(dp), dimension(size(pini))      :: XX
    real(dp), dimension(size(pini))      :: BESTX
    real(dp), dimension(size(pini))      :: WORSTX
    real(dp)                             :: XF(2000)
    real(dp)                             :: CF(2000)
    real(dp)                             :: SF(50)
    real(dp)                             :: DIST(2000)
    real(dp), dimension(size(pini))      :: XI
    real(dp), dimension(size(pini))      :: XNSTD
    real(dp), dimension(size(pini))      :: BOUND
    real(dp)                             :: CRITER(10)
    real(dp)                             :: GNRNG
    real(dp)                             :: DENOMI
    real(dp)                             :: TIMEOU
    real(dp)                             :: BESTF
    real(dp)                             :: WORSTF 
    real(dp)                             :: ranval 
    real(dp)                             :: FA
    character(len=4)                     :: XNAME(100)   ! parameter names - maximum 100 parameters are allowed
    logical, dimension(size(pini))       :: maske        ! parameter to be optimized (true or false)
  
    XNAME=['  X1','  X2','  X3','  X4','  X5','  X6','  X7','  X8','  X9',' X10',&
           ' X11',' X12',' X13',' X14',' X15',' X16',' X17',' X18',' X19',' X20',&
           ' X21',' X22',' X23',' X24',' X25',' X26',' X27',' X28',' X29',' X30',&
           ' X31',' X32',' X33',' X34',' X35',' X36',' X37',' X38',' X39',' X40',&
           ' X41',' X42',' X43',' X44',' X45',' X46',' X47',' X48',' X49',' X50',&
           ' X51',' X52',' X53',' X54',' X55',' X56',' X57',' X58',' X59',' X60',&
           ' X61',' X62',' X63',' X64',' X65',' X66',' X67',' X68',' X69',' X70',&
           ' X71',' X72',' X73',' X74',' X75',' X76',' X77',' X78',' X79',' X80',&
           ' X81',' X82',' X83',' X84',' X85',' X86',' X87',' X88',' X89',' X90',&
           ' X91',' X92',' X93',' X94',' X95',' X96',' X97',' X98',' X99','X100']                                          
     
    if (present(mask)) then
       if (count(mask) .eq. 0_i4b) then
          stop 'Input argument mask: At least one element has to be true'
       else
          maske = mask
       end if
    else
       maske = .true.
    endif

    ! setup SCE variables 
    ! set pval to value from input argument (pini)
    pval  = pini
    pnum  = size(pini)               ! number of parameter values 
    npg   = 2_i4b*count(maske)+1_i4b ! number of point per complex
    nps   = count(maske)+1_i4b       ! number of point per subcomplex (points used complex evolution)
    nspl  = 2_i4b*count(maske)+1_i4b ! number of evovlution steps per complex before shuffling 
    mings = ngs
    npt = ngs * npg                  ! number of total sample points (# of complexes time # of point per complex)

    if (size(prange,1) /= pnum) stop 'Error sceua: size(prange,1) /= size(pini)'
    if (size(prange,2) /= 2)    stop 'Error sceua: size(prange,2) /= 2'
    if (maxiter .le. 1)         stop 'Error sceua: maxiter must be greater than 1'
  
    ! INITIALIZE VARIABLES
    open(ISCE,file=trim(adjustl(tmp_file)), action='write', status='unknown')
    write(ISCE,400)
    
    ICALL = 0
    NLOOP = 0
    LOOP = 0
    IGS = 0

    ! Initialize uniform and gaussian andom numbers with seed - see mo_xor4096.f90 for more
    call xor4096(seed,ranval)
    call xor4096g(seed,ranval)
  
    ! COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
    ngs1 = ngs
    npt1 = npt
  
    ! COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
    do J = 1, pnum
      BOUND(J) = prange(J,2)-prange(J,1)
      XI(J)    = pval(J)
    end do

    ! COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT
    FA = obj_func(pval) 
    ICALL = ICALL + 1
  
    ! PRINT THE INITIAL POINT AND ITS CRITERION VALUE
    write(ISCE,500)
    write(ISCE,510) (XNAME(J),J=1,pnum)
    write(ISCE,520) FA,(pval(J),J=1,pnum)
    close(ISCE)
    if (maxiter .eq. ICALL) return 

    ! STEP 1.1  Generate 1st point 
    ! GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
    ! IF INIFLG IS EQUAL TO 1, SET X(1,.) TO INITIAL POINT pval(:)==pini(:)
    if (INIFLG .eq. 1) then
      do J = 1, pnum
        X(1,J) = pval(J)
      end do
      XF(1) = FA
    ! ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO X(1,.)
    else
      do J = 1, pnum
        if ( maske(J) ) then
          call xor4096(0_i8b,ranval) 
          X(1,J) = prange(J,1) + bound(J) * ranval 
        else
          X(1,J) = pval(J)
        end if
        XX(J) = X(1,J)
      end do
      XF(1) = obj_func(XX) 
    end if
    ICALL = ICALL + 1
    open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
    write(ISCE,645) NLOOP,ICALL,XF(1),(X(1,J),J=1,pnum)
    if (ICALL .ge. maxiter) then  !if maxiter is 2 
      if (XF(1) .lt. FA) then
        BESTF=XF(1)
        do J = 1, pnum
          BESTX(J) = XX(J)
        end do 
      else
        BESTF=FA
        do J = 1, pnum
          BESTX(J) = XI(J)
        end do 
      end if
      !  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      write(ISCE,800) maxiter,LOOP,IGS,NLOOP
      write(ISCE,830)
      write(ISCE,510) (XNAME(J),J=1,pnum)
      write(ISCE,520) BESTF,(BESTX(J),J=1,pnum)
      do J = 1, pnum
        pval(J) = BESTX(J)
      end do
      return 
    end if
    close(ISCE)
  
    ! STEP1.2 GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
    ! SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES
    do I = 2, npt1
      do J = 1, pnum
        if (maske(J)) then
          call xor4096(0_i8b,ranval) 
          X(I,J) = prange(J,1)  + BOUND(J) * ranval 
        else
          X(I,J) = pval(J)
        endif
        XX(J) = X(I,J)
      end do
      XF(I) = obj_func(XX)
      ! PRINT THE RESULTS FOR CURRENT POPULATION
      open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
      write(ISCE,645) NLOOP,ICALL,XF(I),(XX(J),J=1,pnum)
      close(ISCE)
      ICALL = ICALL + 1
      if (ICALL .ge. maxiter) then
        npt1 = I
        exit
      end if
    enddo
  
    ! STEP2- ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
    call SORT(npt1,pnum,X,XF)
    ! RECORD THE BEST AND WORST POINTS
    do J = 1, pnum
      BESTX(J) = X(1,J)
      WORSTX(J) = X(npt1,J)
    end do 
    BESTF = XF(1)
    WORSTF = XF(npt1)
  
    ! COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
    call PARSTT(npt1, pnum, X, maske, BOUND, XNSTD, GNRNG, IPCNVG)
  
  !  COMPUTE THE PARAMETER DISTANCE FROM THE INITIAL POPULATION
    call NORMDIST(npt,pnum,X,XI,DIST,BOUND)
  
  !  PRINT THE RESULTS FOR THE INITIAL POPULATION
    open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
    write(ISCE,600)
    write(ISCE,610) (XNAME(J),J=1,pnum)
    write(ISCE,630) NLOOP,ICALL,ngs1,BESTF,WORSTF,DIST(1),(BESTX(J),J=1,pnum)
    if (IPRINT .EQ. 1) then 
      write(ISCE,650) NLOOP
      write(ISCE,615) (XNAME(J),J=1,pnum)
      do I = 1, npt1
        write(ISCE,620) XF(I),(X(I,J),J=1,pnum)
      end do
    end if
    close(ISCE)

    if (ICALL .ge. maxiter) then 
      !  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
      write(ISCE,800) maxiter,LOOP,IGS,NLOOP
      write(ISCE,830)
      write(ISCE,510) (XNAME(J),J=1,pnum)
      write(ISCE,520) BESTF,(BESTX(J),J=1,pnum)
      close(ISCE)
      do J = 1, pnum
        pval(J) = BESTX(J)
      end do
      return 
    end if
  
    if (IPCNVG .eq. 1) then 
      !  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
      write(ISCE,820) GNRNG*100.
      write(ISCE,830)
      write(ISCE,510) (XNAME(J),J=1,pnum)
      write(ISCE,520) BESTF,(BESTX(J),J=1,pnum)
      close(ISCE)
      do J = 1, pnum
        pval(J) = BESTX(J)
      end do
      return 
    end if

  !  BEGIN THE MAIN LOOP ----------------
    main:do
      NLOOP = NLOOP + 1
      ! BEGIN LOOP ON COMPLEXES ----------------
      complexes:do IGS = 1, ngs1
        ! STEP3- ASSIGN POINTS INTO COMPLEXES
        do K1 = 1, npg
          K2 = (K1-1) * ngs1 + IGS
          do J = 1, pnum
            CX(K1,J) = X(K2,J)
          end do
          CF(K1) = XF(K2)
        end do
        ! BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
        subcomplex:do LOOP = 1, nspl
          ! STEP 4. COMPLEX EVOLUTION STEP
          ! CHOOSE A SUB-COMPLEX (nps POINTS) ACCORDING TO A LINEAR
          ! PROBABILITY DISTRIBUTION
          if (nps .EQ. npg) then
            DO K = 1, nps
              LCS(K) = K
            end do
            ! CREATE THE SUB-COMPLEX ARRAYS
            do K = 1, nps
              do J = 1, pnum
                S(K,J) = CX(LCS(K),J)
              end do
              SF(K) = CF(LCS(K))
            end do
          else
            call xor4096(0_i8b,ranval) 
            LCS(1) = 1 + int(npg + 0.5 - SQRT( (npg+.5)**2 - npg * (npg+1) * ranval ))
            do K = 2, nps
              do 
                call xor4096(0_i8b,ranval) 
                LPOS = 1 + int(npg + 0.5 - SQRT((npg+.5)**2 - npg * (npg+1) * ranval ))
                do K1 = 1, K-1
                  if (LPOS .EQ. LCS(K1)) cycle 
                end do
                exit
              end do
              LCS(K) = LPOS
            end do
            !  ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
            call SORT1(nps,LCS)
            !  CREATE THE SUB-COMPLEX ARRAYS
            do K = 1, nps
              do J = 1, pnum
                S(K,J) = CX(LCS(K),J)
              end do
              SF(K) = CF(LCS(K))
            end do
          endif

          ! USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
          call CCE(obj_func,pnum,nps,S,SF,prange,XNSTD,ICALL,maxiter,maske)
  
          ! IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
          ! INTO THE COMPLEX
          do K = 1, nps
            do J = 1, pnum
              CX(LCS(K),J) = S(K,J)
            end do
            CF(LCS(K)) = SF(K)
          end do
  
          ! SORT THE POINTS
          call SORT(npg,pnum,CX,CF)
  
         ! ! RECORD THE BEST AND WORST POINTS
         ! do J = 1, pnum
         !   BESTX(J) = CX(1,J)
         !   WORSTX(J) = CX(npt1,J)
         ! end do
         ! BESTF = CF(1)
         ! WORSTF = CF(npt1)
         ! ! PRINT THE RESULTS FOR CURRENT POPULATION
         ! write(ISCE,640) NLOOP,ICALL,ngs1,BESTF,WORSTF,(BESTX(J),J=1,pnum)
          ! IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
          if (ICALL .GE. maxiter) exit
  
        end do subcomplex 

        ! STEP 5. Suffle complexes
        ! REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY X(.,.)
        do K1 = 1, npg
          K2 = (K1-1) * ngs1 + IGS
          do J = 1, pnum
            X(K2,J) = CX(K1,J)
          end do
          XF(K2) = CF(K1)
        end do 
        if (ICALL .GE. maxiter) exit
  
      end do complexes ! end loop on complexes

      ! RE-SORT THE POINTS
      call SORT(npt1,pnum,X,XF)
  
      ! RECORD THE BEST AND WORST POINTS
      do J = 1, pnum
        BESTX(J) = X(1,J)
        WORSTX(J) = X(npt1,J)
      end do
      BESTF = XF(1)
      WORSTF = XF(npt1)
  
      ! TEST THE POPULATION FOR PARAMETER CONVERGENCE
      call PARSTT(npt1, pnum, X, maske, BOUND, XNSTD, GNRNG, IPCNVG)
  
      ! COMPUTE THE PARAMETER DISTANCE FROM THE INITIAL POPULATION
      call NORMDIST(npt,pnum,X,XI,DIST,BOUND)
  
      ! PRINT THE RESULTS FOR CURRENT POPULATION
      open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
      write(ISCE,610) (XNAME(J),J=1,pnum)
      write(ISCE,630) NLOOP,ICALL,ngs1,BESTF,WORSTF,DIST(1),(BESTX(J),J=1,pnum)
      if (IPRINT .EQ. 1) THEN
        write(ISCE,650) NLOOP
        write(ISCE,615) (XNAME(J),J=1,pnum)
        do I = 1, npt1
          write(ISCE,620) XF(I),(X(I,J),J=1,pnum)
        end do
      end if 
  
      ! TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
      if (ICALL .GE. maxiter) then
        ! PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
        open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
        write(ISCE,800) maxiter,LOOP,IGS,NLOOP
        write(ISCE,830)
        write(ISCE,510) (XNAME(J),J=1,pnum)
        write(ISCE,520) BESTF,(BESTX(J),J=1,pnum)
        close(ISCE) 
        do J = 1, pnum
          pval(J) = BESTX(J)
        end do
        return
      end if
  
      ! COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
      CRITER(10) = BESTF
      if (NLOOP .lt. (kstop+1)) then 
        do LDX = 1, 9
          CRITER(LDX) = CRITER(LDX+1)
        end do
      else
        DENOMI = abs(CRITER(10-kstop) + CRITER(10)) / 2.
        TIMEOU = abs(CRITER(10-kstop) - CRITER(10)) / DENOMI
        if (TIMEOU .lt. pcento) then 
          !  PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
          open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
          write(ISCE,810) pcento*100.,kstop
          write(ISCE,830)
          write(ISCE,510) (XNAME(J),J=1,pnum)
          write(ISCE,520) BESTF,(BESTX(J),J=1,pnum)
          close(ISCE) 
          do J = 1, pnum
            pval(J) = BESTX(J)
          end do
          return
        endif
        do LDX = 1, 9
          CRITER(LDX) = CRITER(LDX+1)
        end do
      end if
  
      ! IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
      if (IPCNVG .eq. 1) then 
        ! PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
        open(unit=ISCE,file=trim(adjustl(tmp_file)), action='write', position='append')
        write(ISCE,820) GNRNG*100.
        write(ISCE,830)
        write(ISCE,510) (XNAME(J),J=1,pnum)
        write(ISCE,520) BESTF,(BESTX(J),J=1,pnum)
        close(ISCE) 
        do J = 1, pnum
          pval(J) = BESTX(J)
        end do
        return
      end if
  
      ! NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH
      ! CHECK FOR COMPLEX NUMBER REDUCTION
      if (ngs1 .gt. mings) then
        ngs2 = ngs1
        ngs1 = ngs1 - 1
        npt1 = ngs1 * npg
        call COMP(pnum,npt1,ngs1,ngs2,npg,X,XF,CX,CF)
      end if
  
    end do main  ! end of main loop

    return
  
    400 format(//,2X,50(1H=),/,2X,'ENTER THE SHUFFLED COMPLEX EVOLUTION GLOBAL SEARCH',/,2X,50(1H=))
    500 format(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION VALUE ***')
    510 format(/,' CRITERION',100(2X,A4,2X),/1X,80(1H-))
    520 format(F10.3,100F8.3)
    600 format(//,1X,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
    610 format(/,1X,'LOOP',2X,'TRIALS',2X,'COMPLXS',5X,'BEST-F',4X,'WORST-F',4X,'PAR-RNG',4X,100(A4,2X))
    615 format(9X,'F',4X,100(A4,2X))
    620 format(F10.3,100(F8.3))
    630 format(I5,3X,I5,3X,I6,1X,3(F10.3,1X),100(F8.3))
    640 format(I5,1X,I5,3X,I5,1X,2F10.3,8x,100(F8.3))
    645 format(I5,1X,I5,8X,F10.3,16X,100(F8.3))
    650 format(/,1X,'POPULATION AT LOOP ',I3,/,1X,22(1H-))
    800 format(//,1X,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',  &
               ' LIMIT ON THE MAXIMUM',/,5X,'NUMBER OF TRIALS ',I5,     &
               ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5X,'SUB-COMPLEX ', &
               I3,' OF COMPLEX ',I3,' IN SHUFFLING LOOP ',I3,' ***')
    810 format(//,1X,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION', &
               ' VALUE HAS NOT CHANGED ',/,5X,F5.3,' PERCENT IN',I3,' SHUFFLING LOOPS ***')
    820 format(//,1X,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',&
               ' HAS CONVERGED INTO ',/,4X,F5.3,' PERCENT OF THE FEASIBLE SPACE ***')
    830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS CRITERION VALUE ***')
  
  end subroutine

!====================================================================
  subroutine CCE(obj_func,pnum,nps,S,SF,prange,XNSTD,ICALL,maxiter,maske)
    use mo_xor4096, only: xor4096g
  ! ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
  ! SUB-COMPLEX VARIABLES
    implicit none
    interface 
      function obj_func(pp)
        use nrtype
        implicit none
        real(dp), intent(in) :: pp(:)
        real(dp)             :: obj_func
      end function obj_func
    end interface
    integer(i4b),          intent(in)    :: pnum
    integer(i4b),          intent(in)    :: nps 
    real(dp),              intent(in)    :: prange(:,:) 
    real(dp),              intent(in)    :: XNSTD(:)
    integer(i8b),          intent(in)    :: maxiter 
    real(dp),              intent(inout) :: S(:,:)
    real(dp),              intent(inout) :: SF(:)
    integer(i4b),          intent(inout) :: ICALL
    logical,               intent(in)    :: maske(:)
    ! local
    integer(i4b)                         :: N
    integer(i4b)                         :: I,J
    integer(i4b)                         :: ibound 
    real(dp)                             :: zvalue 
    real(dp)                             :: FW
    real(dp)                             :: FNEW       ! FUNCTION VALUE OF THE WORST POINT
    real(dp)                             :: WO(pnum)  ! THE WORST POINT OF THE SIMPLEX
    real(dp)                             :: CE(pnum)  ! THE CENTROID OF THE SIMPLEX EXCLUDING WO
    real(dp)                             :: SNEW(pnum)! NEW POINT GENERATED FROM THE SIMPLEX
    real(dp)                             :: STEP(pnum)! VECTOR FROM WO TO CE

    ! EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
    N = nps

    ! IDENTIFY THE WORST POINT WO OF THE SUB-COMPLEX S
    ! COMPUTE THE CENTROID CE OF THE REMAINING POINTS
    ! COMPUTE STEP, THE VECTOR BETWEEN WO AND CE
    ! IDENTIFY THE WORST FUNCTION VALUE FW
    do J = 1, pnum 
      WO(J) = S(N,J)
      CE(J) = 0.0
      do I = 1, N-1
        CE(J) = CE(J) + S(I,J)
      end do
      CE(J) = CE(J)/dble(N-1)
      STEP(J) = CE(J) - WO(J)
    end do
    FW = SF(N)

   ! COMPUTE THE NEW POINT SNEW

   ! FIRST TRY A REFLECTION STEP
    do J = 1, pnum 
      SNEW(J) = WO(J) + 2. * STEP(J)
    end do
    ! CHECK IF SNEW IS WITHIN BOUND OR NOT
    ibound = 0
    do J = 1, pnum 
      if (SNEW(J) .GT. prange(J,2) .or. SNEW(J) .LT. prange(J,1)) then
        ibound = 1
        exit 
      end if
    end do
    ! if SNEW is outside the bound,
    ! CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
    ! A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
    ! AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
    if (ibound .eq. 1) then
      do J = 1, pnum 
        if (maske(J)) then
          do
            call xor4096g(0_i8b,zvalue)
            SNEW(J) = S(1,J) + XNSTD(J)*zvalue*(prange(J,2)-prange(J,1))
            if (SNEW(J) .le. prange(J,2) .and. SNEW(J) .ge. prange(J,1)) exit
          end do
        end if
      end do 
    end if 
    ! COMPUTE THE FUNCTION VALUE AT SNEW
    FNEW = obj_func(SNEW)
    ICALL = ICALL + 1
    !COMPARE FNEW WITH THE WORST FUNCTION VALUE FW
    ! FNEW IS LESS THAN FW, ACCEPT THE NEW POINT SNEW AND RETURN
    if (FNEW .le. FW) then 
      do J = 1, pnum 
        S(N,J) = SNEW(J)
      end do
      SF(N) = FNEW
      return
    end if
    if (ICALL .ge. maxiter) return

   ! FNEW IS GREATER THAN FW, SO TRY A CONTRACTION STEP
    do J = 1, pnum 
      SNEW(J) = WO(J) + 0.5 * STEP(J)
    end do
    ! COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
    FNEW = obj_func(SNEW)
    ICALL = ICALL + 1
    ! COMPARE FNEW TO THE WORST VALUE FW
    ! IF FNEW IS LESS THAN OR EQUAL TO FW, THEN ACCEPT THE POINT AND RETURN
    if (FNEW .le. FW) then
      do J = 1, pnum 
        S(N,J) = SNEW(J)
      end do
      SF(N) = FNEW
      return
    end if
    if (ICALL .ge. maxiter) return

    ! IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
    ! ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
    ! AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
    do J = 1, pnum 
      if (maske(J)) then
        do
          call xor4096g(0_i8b,zvalue)
          SNEW(J) = S(1,J) + XNSTD(J)*zvalue*(prange(J,2)-prange(J,1))
          if (SNEW(J) .le. prange(J,2) .and. SNEW(J) .ge. prange(J,1)) exit
        end do
      end if
    end do
    ! COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
    FNEW = obj_func(SNEW)
    ICALL = ICALL + 1
    ! REPLACE THE WORST POINT BY THE NEW POINT
    DO J = 1, pnum 
      S(N,J) = SNEW(J)
    END DO
    SF(N) = FNEW

    return
  end subroutine

!===================================================================
  subroutine PARSTT(npt, pnum, X, maske, BOUND, XNSTD, GNRNG, IPCNVG)
!  SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
    implicit none
    !input
    integer(i4b),  intent(in)  :: npt
    integer(i4b),  intent(in)  :: pnum
    real(dp),      intent(in)  :: X(:,:)
    logical,       intent(in)  :: maske(:)
    real(dp),      intent(in)  :: BOUND(:)
    !output
    real(dp),      intent(out) :: XNSTD(:)
    real(dp),      intent(out) :: GNRNG
    integer(i4b),  intent(out) :: IPCNVG
    !local
    real(dp)                   :: XMAX(pnum)
    real(dp)                   :: XMIN(pnum)
    real(dp)                   :: XMEAN(pnum)
    real(dp)                   :: XSUM1
    real(dp)                   :: XSUM2
    real(dp)                   :: gsum
    real(dp),parameter         :: delta=1.e-20_dp
    real(dp),parameter         :: peps=1.e-5_dp
    integer(i4b)               :: IDX,KDX ! loop indices 

    ! Initinalize XNSTD
    XNSTD=0.0_dp
    ! COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
    gsum = 0.0_dp
    do KDX = 1, pnum
      ! Account for only adjusting parameter to compute gsum
      if (maske(KDX)) then
        XMAX(KDX) = -1.0e+20_dp
        XMIN(KDX) = 1.0e+20_dp
        XSUM1 = 0.0_dp
        XSUM2 = 0.0_dp
        do IDX = 1, npt
          XMAX(KDX) = max(X(IDX,KDX), XMAX(KDX))
          XMIN(KDX) = min(X(IDX,KDX), XMIN(KDX))
          XSUM1 = XSUM1 + X(IDX,KDX)
          XSUM2 = XSUM2 + X(IDX,KDX)*X(IDX,KDX)
        end do
        ! mean of X(KDX) across sample points
        XMEAN(KDX) = XSUM1 / dble(npt)
        ! std of X(KDX) across sample points
        XNSTD(KDX) = (XSUM2 / dble(npt) - XMEAN(KDX)*XMEAN(KDX))
        if (XNSTD(KDX) .le. delta) XNSTD(KDX) = delta
        XNSTD(KDX) = sqrt(XNSTD(KDX))
        ! Fraction of X(KDX) std to Xrange(KDX)
        XNSTD(KDX) = XNSTD(KDX) / BOUND(KDX)
        gsum = gsum + log( delta + (XMAX(KDX)-XMIN(KDX))/BOUND(KDX) )
      end if
    end do
    GNRNG = exp( gsum/dble(count(maske)) )
!  CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= EPS
    IPCNVG = 0
    if (GNRNG .le. peps) then
        IPCNVG = 1
    end if
    return
  end subroutine

!===================================================================
  subroutine NORMDIST(npt,pnum,X,XI,DIST,BOUND)
!  SUBROUTINE COMPUTING NORMAILZIED DISTANCE FROM INITIAL POINT
!     X(.,.)  - POPULATION
!     XI(.)   - INITIAL POINT
!     DIST(.) - NORMALIZED DISTANCE FROM INITIAL POINT
    implicit none
    integer(i4b),  intent(in)  :: npt
    integer(i4b),  intent(in)  :: pnum
    real(dp),      intent(in)  :: X(:,:)
    real(dp),      intent(in)  :: XI(:)
    real(dp),      intent(in)  :: BOUND(:)
    real(dp),      intent(out) :: DIST(:)
    ! local variables
    integer(i4b)               :: IDX, KDX ! loop indices 

!  COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
    do KDX = 1, npt
       DIST(KDX) = 0.
       do IDX = 1, pnum
          DIST(KDX) = DIST(KDX) + ABS(X(KDX,IDX) - XI(IDX))/BOUND(IDX)
       end do
       DIST(KDX) = DIST(KDX) / pnum
    end do 
    return
  end subroutine

!====================================================================
  subroutine COMP(N,npt,ngs1,ngs2,npg,A,AF,B,BF)
!  THIS SUBROUTINE REDUCE INPUT MATRIX A(N,ngs2*npg) TO MATRIX
!  B(N,ngs1*npg) AND VECTOR AF(ngs2*npg) TO VECTOR BF(ngs1*npg)
    implicit none
    integer(i4b),  intent(in)    :: N
    integer(i4b),  intent(in)    :: npt
    integer(i4b),  intent(in)    :: ngs1
    integer(i4b),  intent(in)    :: ngs2
    integer(i4b),  intent(in)    :: npg
    real(dp),      intent(inout) :: A(:,:)
    real(dp),      intent(inout) :: AF(:)
    real(dp),      intent(out)   :: B(:,:)
    real(dp),      intent(out)   :: BF(:)
    !local
    integer(i4b)                 :: IGS,IPG  !loop indices
    integer(i4b)                 :: K1,K2
    integer(i4b)                 :: IDX,JDX

    do IGS=1, ngs1
      do IPG=1, npg
        K1=(IPG-1)*ngs2 + IGS
        K2=(IPG-1)*ngs1 + IGS
        do IDX=1, N
          B(K2,IDX) = A(K1,IDX)
        end do
        BF(K2) = AF(K1)
      end do
    end do
    do JDX=1, npt
      do IDX=1, N
        A(JDX,IDX) = B(JDX,IDX)
      end do
      AF(JDX) = BF(JDX)
    end do
    return
  end subroutine

!===================================================================
  subroutine SORT(N,M,RB,RA)
!  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!  BY W.H. PRESS ET AL., PP. 233-234
!  LIST OF VARIABLES
!     RA(.) = ARRAY TO BE SORTED
!     RB(.,.) = ARRAYS ORDERED CORRESPONDING TO REARRANGEMENT OF RA(.)
!     WK(.,.), IWK(.) = LOCAL VARIBLES
    implicit none
    !input
    integer(i4b),intent(in)    :: N,M 
    real(dp),    intent(inout) :: RA(:) 
    !output
    real(dp),    intent(out)   :: RB(:,:) 
    !local variable
    real(dp)    ,allocatable   :: WK(:,:) 
    integer(i4b),allocatable   :: INDX(:) 
    integer(i4b),allocatable   :: IWK(:) 
    integer(i4b)               :: I,J 
    
    allocate( WK(size(RB,1),size(RB,2)) ) 
    allocate( INDX(size(RA)) )
    allocate( IWK(size(RA)) )

    call INDEXX(N, RA, IWK)
    do I = 1, N
      WK(I,1) = RA(I)
    END DO
    do I = 1, N
      RA(I) = WK(IWK(I),1)
    end do
    do J = 1, M
      do I = 1, N
        WK(I,J) = RB(I,J)
      end do
    end do
    do J = 1, M
      do I = 1, N
        RB(I,J) = WK(IWK(I),J)
      end do
    end do

    RETURN
  end subroutine

!===========================================================
  subroutine SORT1(N,RA)
!  SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
!  BY W.H. PRESS ET AL., PP. 231
    implicit none
    !input
    integer(i4b), intent(in)     :: N
    integer(i4b), intent(inout)  :: RA(N) 
    !local
    integer(i4b)                 :: I, J     ! loop index
    integer(i4b)                 :: L,IR
    integer(i4b)                 :: RRA 

    L = (N / 2) + 1
    IR = N
    do
      if (L .gt. 1) then
        L = L - 1
        RRA = RA(L)
      else
        RRA = RA(IR)
        RA(IR) = RA(1)
        IR = IR - 1
        IF (IR .EQ. 1) then
          RA(1) = RRA
          return
        end if
      end if
      I = L
      J = L + L
      do
        if (J .LE. IR) then
          if (J .LT. IR) then
            if (RA(J) .lt. RA(J + 1)) J = J + 1
          end if
          if (RRA .LT. RA(J)) then
            RA(I) = RA(J)
            I = J
            J = J + J
          else
            J = IR + 1
          end if
        else
          exit
        end if
      end do
      RA(I) = RRA
    end do 

  end subroutine

!=======================================================
  subroutine INDEXX(N, ARRIN, INDX)
  ! THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
    implicit none
    !input
    integer(i4b), intent(in)  :: N
    real(dp),     intent(in)  :: ARRIN(N) 
    !output
    integer(i4b), intent(out) :: INDX(N) 
    !local 
    integer(i4b)              :: I, J     ! loop index
    integer(i4b)              :: L,IR
    real(i4b)                 :: Q
    integer(i4b)              :: INDXT 
     
    do J = 1, N
      INDX(J) = J
    end do
    L = (N / 2) + 1
    IR = N
    do
      if (L .gt. 1) then
        L = L - 1
        INDXT = INDX(L)
        Q = ARRIN(INDXT)
      else
        INDXT = INDX(IR)
        Q = ARRIN(INDXT)
        INDX(IR) = INDX(1)
        IR = IR - 1
        if (IR .eq. 1) then
          INDX(1) = INDXT
          return
        end if
      end if
      I = L
      J = L + L
      do
        if (J .LE. IR) then
          if (J .LT. IR) then
            if (ARRIN(INDX(J)) .LT. ARRIN(INDX(J + 1))) J = J + 1
          end if
          if (Q .lt. ARRIN(INDX(J))) then
            INDX(I) = INDX(J)
            I = J
            J = J + J
          else
            J = IR + 1
          end if
        else 
          exit
        end if
      end do
      INDX(I) = INDXT
    end do 

  end subroutine

end module sce
