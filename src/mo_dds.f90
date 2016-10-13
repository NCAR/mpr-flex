!> \file mo_dds.f90

!> \brief Dynamically Dimensioned Search (DDS)

!> \details This module provides routines for Dynamically Dimensioned Search (DDS)
!> of Tolson and Shoemaker (2007). It searches the minimum or maximum of a user-specified function,
!> using an n-dimensional continuous global optimization algorithm (DDS).

!> \authors Original by Bryan Tolson and later modified by Rohini Kumar.
!> Matthias Cuntz and Juliane Mai for the module, MDDS, etc.
!> \date Jul 2012

module mo_dds

  ! This module contains routines for Dynamically Dimensioned Search (DDS)

  ! Written Jul 2012, Matthias Cuntz - module version of modified DDS v1.1 of R. Kumar to
  ! original DDS of Bryan Tolson

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012 Matthias Cuntz

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: DDS    ! Dynamically Dimensioned Search (DDS)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !>        \brief DDS

  !     PURPOSE
  !>        \details Searches Minimum or Maximum of a user-specified function using
  !>        Dynamically Dimensioned Search (DDS).
  !>
  !>        DDS is an n-dimensional continuous global optimization algorithm.
  !>        It is coded as a minimizer but one can give maxit=True in a maximization problem,
  !>        so that the algorithm minimizes the negative of the objective function F=(-1*F).

  !>        The function to be minimized is the first argument of DDS and must be defined as \n
  !>        \code
  !>            function func(p)
  !>              use nrtype 
  !>              implicit none
  !>              real(dp), dimension(:), intent(in) :: p
  !>              real(dp) :: func
  !>            end function func
  !>        \endcode

  !     CALLING SEQUENCE
  !         popt = DDS(obj_func, pini, prange, r=r, seed=seed, maxiter=maxiter, maxit=maxit, funcbest=funcbest)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: obj_func(p)"             Function on which to search the minimum
  !>        \param[in] "real(dp) :: pini(:)"                 inital value of decision variables
  !>        \param[in] "real(dp) :: prange(size(pini),2)"    Min/max range of decision variables
  !>        \logical,            :: restart                  logical to read state file and initialize, .false. -> start from begining

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp), optional           :: r"                 DDS perturbation parameter\n
  !>                                                                       (default: 0.2)
  !>        \param[in] "integer(i8b), optional        :: seed"              User seed to initialise the random number generator
  !>                                                                       (default: None)
  !>        \param[in] "integer(i8b), optional        :: maxiter"           Maximum number of iteration or function evaluation
  !>                                                                       (default: 1000)
  !>        \param[in] "logical, optional            :: maxit"             Maximization (.True.) or 
  !>                                                                       minimization (.False.) of function
  !>                                                                       (default: .False.)
  !>        \param[in] "logical, optional            :: mask(size(pini))"  parameter to be optimized (true or false)
  !>                                                                       (default: .True.)
  !>        \param[in]  "character(len=*) , optional :: tmp_file"          file with temporal output

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), optional              :: funcbest"    the best value of the function.
  !>        \param[out] "real(dp), optional, allocatable :: history(:)"  the history of best function values,
  !>                                                                     history(maxiter)=funcbest\n
  !>                                          allocatable only to be in correspondance with other optimization routines

  !     RETURN
  !>        \return real(dp) :: DDS  &mdash;  The parameters of the point which is estimated to minimize the function.

  !     RESTRICTIONS
  !         None.

  !     EXAMPLE
  !         dv_range(:,1) = (/ -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0 /)
  !         dv_range(:,2) = (/ 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0 /)
  !         dv_ini        = (/ -.226265E+01, -.130187E+01, -.151219E+01, 0.133983E+00, 0.988159E+00, &
  !                            -.495074E+01, -.126574E+02, 0.572684E+00, 0.303864E+01, 0.343031E+01 /)
  !         dv_opt = DDS(griewank, dv_ini, dv_range)
  !         -> see also example in test directory

  !     LITERATURE
  !         Tolson, B. A., and C. A. Shoemaker (2007)
  !             Dynamically dimensioned search algorithm for computationally efficient watershed
  !             model calibration, Water Resour. Res., 43, W01413, doi:10.1029/2005WR004723.

  !     HISTORY
  !>        \author Written original Bryan Tolson - DDS v1.1\n
  !>        Modified Rohini Kumar, Matthias Cuntz, Juliane Mai
  !         \date Feb 2007
  !         Modified, Rohini Kumar, Feb 2008
  !         Modified, Matthias Cuntz & Juliane Mai, Jul 2012 - module
  !                   Juliane Mai,                  Aug 2012 - optional argument funcbest added
  !                   Juliane Mai,                  Nov 2012 - masked parameter
  !                   Juliane Mai,                  Dec 2012 - history output

  subroutine DDS(obj_func, pini, prange, restart, restartFile, r, seed, maxiter, maxit, mask, tmp_file, history)

    use nrtype
    use public_var
    use mo_xor4096, only: xor4096, xor4096g

    implicit none

    INTERFACE
       function obj_func(pp)
         use nrtype 
         implicit none
         real(dp), dimension(:), intent(in) :: pp
         real(dp) :: obj_func
       end function obj_func
    END INTERFACE
    real(dp),    dimension(:),             intent(in)  :: pini        ! inital value of decision (parameter) variables
    real(dp),    dimension(:,:),           intent(in)  :: prange      ! Min/max values of decision variables
    logical(lgc),                          intent(in)  :: restart     ! .true.  read state file and initialize, .false. -> start from begining
    character(len=strLen),                 intent(in)  :: restartFile ! name of restart file including iteration, the most recent parameter values 
    real(dp),                    optional, intent(in)  :: r           ! DDS perturbation parameter (-> 0.2 by default)
    integer(i8b),                optional, intent(in)  :: seed        ! User seed to initialise the random number generator
    integer(i8b),                optional, intent(in)  :: maxiter     ! Maximum number of iteration or function evaluation
    logical,                     optional, intent(in)  :: maxit       ! .false. -> minimization of function, or .true. -> Maximization
    logical,     dimension(:),   optional, intent(in)  :: mask        ! parameter to be optimized (true or false)
    character(len=strLen),       optional, intent(in)  :: tmp_file    ! file for temporal output
    real(dp),    dimension(:),   optional, intent(out), &
         allocatable                                   :: history     ! History of objective function values

    ! Local variables
    integer(i4b)                            :: pnum                   ! Total number of decision variables
    real(dp), dimension(size(pini))         :: pval                   ! initial value of decision variables - either input or restart
    integer(i8b)                            :: iseed                  ! User given seed
    integer(i8b)                            :: imaxiter               ! Maximum number of iteration or function evaluation
    integer(i8b)                            :: iStart                 ! starting index of objective function evaluation 
    real(dp)                                :: ir                     ! DDS perturbation parameter
    real(dp)                                :: imaxit                 ! Maximization or minimization of function
    real(dp)                                :: of_new, of_best        ! intermediate results for objective function values
    real(dp)                                :: Pn, new_value          ! intermediate results for parameter value  
    real(dp), dimension(size(pini))         :: pnew                   ! Test value of decision (parameter) variables
    real(dp), dimension(size(pini))         :: pbest                  ! Best value of decision (parameter) variables 
    real(dp)                                :: ranval                 ! random value
    integer(i8b)                            :: i                      ! maxiter=i8b
    integer(i4b)                            :: j, dvn_count, dv       ! pnum=i4b
    integer(i4b)                            :: idummy                 ! dummy vaiable
    integer, dimension(8)                   :: sdate                  ! date_and_time return
    logical, dimension(size(pini))          :: maske                  ! parameter to be optimized (true or false)
    integer(i4b), dimension(:), allocatable :: truepara               ! parameter to be optimized (their indexes)
    logical                                 :: isExistFile            ! logical to check if the file exist or not
    character(50)                           :: rowfmt1                ! string specifying write format for real value
    character(50)                           :: rowfmt2                ! string specifying write format for real value
    
    ! set starting iteration to 1 and pval to value from input argument (pini)
    iStart=1
    pval=pini

    ! Check input
    pnum = size(pini)
    if (size(prange,1) /= pnum) stop 'Error DDS: size(prange,1) /= size(pini)'
    if (size(prange,2) /= 2)    stop 'Error DDS: size(prange,2) /= 2'
    ! restart option?
    if (restart) then
      print*, 'read restart file'
      inquire(file=trim(adjustl(restartFile)), exist=isExistFile)
      if ( isExistFile ) then !  if state file exists, update iStart and pval, otherwise iteration start with very beginning
        open(unit=70,file=trim(adjustl(restartFile)), action='read', status = 'unknown')
        read(70,*) iStart 
        read(70,*) (pval(i),i=1,pnum)    
        close(70)
        print*, iStart
       endif
    endif 
    ! r Perturbation parameter
    ir = 0.2_dp
    if (present(r)) ir = r
    if (ir <= 0.0_dp .or. ir > 1.0_dp) stop 'Error DDS: DDS perturbation parameter (0.0, 1.0]'
    ! max. iteration
    imaxiter = 1000
    if (present(maxiter)) imaxiter = maxiter
    if (imaxiter < 6) stop 'Error DDS: max function evals must be minimum 6'
    ! history output
    if (present(history)) then
       allocate(history(imaxiter))
    end if
    ! Min or max objective function
    imaxit = 1.0_dp
    if (present(maxit)) then
       if (maxit) imaxit = -1.0_dp
    endif
    ! Given seed
    iseed = 0
    if (present(seed)) iseed = seed
    iseed = max(iseed, 0_i8b)

    if (present(mask)) then
       if (count(mask) .eq. 0_i4b) then
          stop 'Input argument mask: At least one element has to be true'
       else
          maske = mask
       end if
    else
       maske = .true.
    endif

    allocate ( truepara(count(maske)) )
    idummy = 0_i4b
    do j=1,size(pval,1)
       if ( maske(j) ) then
          idummy = idummy+1_i4b
          truepara(idummy) = j
       end if
    end do

    ! Seed random numbers
    if (iseed == 0) then
       call date_and_time(values=sdate)
       iseed = sdate(1)*31536000000_i8b + sdate(2)*2592000000_i8b + sdate(3)*86400000_i8b + &
            sdate(5)*3600000_i8b + sdate(6)*60000_i8b + sdate(7)*1000_i8b + sdate(8)
       call xor4096(iseed,ranval)
       call xor4096g(iseed,ranval)
    else
       call xor4096(iseed,ranval)
       call xor4096g(iseed,ranval)
    endif

    ! Temporal file writing
    if(present(tmp_file) .and. (.not.restart) ) then
       open(unit=999,file=trim(adjustl(tmp_file)), action='write', status = 'unknown')
       write(999,*) '# settings :: general'
       write(999,*) '# nIterations    iseed'
       write(999,*) imaxiter, iseed
       write(999,*) '# settings :: dds specific'
       write(999,*) '# dds_r'
       write(999,*) ir
       write(999,*) '# iter   bestf   (bestx(j),j=1,nopt)'
       close(999)
    end if

    ! Evaluate initial solution and return objective function value
    ! and Initialise the other variables (e.g. of_best)
    ! imaxit is 1.0 for MIN problems, -1 for MAX problems
    pbest      = pval
    of_new     = imaxit * obj_func(pval)
    of_best    = of_new
    if (present(history)) history(1) = of_new

    file_write: if (present(tmp_file)) then
       open(unit=999,file=trim(adjustl(tmp_file)), action='write', position='append', recl=(pnum+2)*30)
       write(rowfmt1,'(A,I5,A)') '(I6,1X',(pnum+1),'(1X,ES17.10))'
       if (imaxit .lt. 0.0_dp) then  ! Maximize
          write(999,rowfmt1) iStart-1, -of_best, pval
       else                          ! Minimize
          write(999,rowfmt1) iStart-1, of_best, pval
       end if
       close(999)
    end if file_write

    ! Code below is now the DDS algorithm as presented in Figure 1 of Tolson and Shoemaker (2007)
    do i=iStart, imaxiter-1
       ! Determine Decision Variable (DV) selected for perturbation:
       Pn        = 1.0_dp - log(real(i,dp))/log(real(imaxiter-1,dp)) ! probability each DV selected
       dvn_count = 0                                                 ! counter for how many DVs selected for perturbation
       pnew      = pbest                                               ! define pnew initially as best current solution

       ! Step 3 of Fig 1 of Tolson and Shoemaker (2007)
       do j=1, size(truepara) !pnum
          call xor4096(0_i8b,ranval)                           ! selects next uniform random number in sequence
          ! Step 4 of Fig 1 of Tolson and Shoemaker (2007)
          if (ranval < Pn) then                               ! jth DV selected for perturbation
             dvn_count = dvn_count + 1
             ! call 1-D perturbation function to get new DV value (new_value)
             call neigh_value(pbest(truepara(j)), prange(truepara(j),1), prange(truepara(j),2), ir, new_value)
             pnew(truepara(j)) = new_value
          end if
       end do

       ! Step 3 of Fig 1 of Tolson and Shoemaker (2007) in case {N} empty
       if (dvn_count == 0) then                               ! no DVs selected at random, so select one
          call xor4096(0_i8b,ranval)                           ! selects next uniform random number in sequence
          dv = truepara(int(( ranval * real(size(truepara),dp))  + 1.0_dp, i4b ))  ! index for one DV
          ! call 1-D perturbation function to get new DV value (new_value):
          call neigh_value(pbest(dv), prange(dv,1), prange(dv,2), ir, new_value)
          pnew(dv) = new_value                                ! change relevant DV value in stest
       end if

       ! Step 5 of Fig 1 of Tolson and Shoemaker (2007)
       ! Evaluate obj function value for test
       of_new = imaxit * obj_func(pnew)                       ! imaxit handles min(=1) and max(=-1) problems
       ! update current best solution
       if (of_new <= of_best) then
          of_best = of_new
          pbest     = pnew
       end if
       if (present(history)) history(i+1) = of_best
       
       ! Write out in temp file if exist
       file_write2: if (present(tmp_file)) then
          open(unit=999,file=trim(adjustl(tmp_file)), action='write', position='append')
          if (imaxit .lt. 0.0_dp) then       ! Maximize
             write(999,rowfmt1) i, -of_best, pbest 
          else                               ! Minimize
             write(999,rowfmt1) i, of_best, pbest 
          end if
          close(999)
       end if file_write2

       ! update out in restart file 
       open(unit=80,file=trim(adjustl(restartFile)), action='write', status='replace')
       write(rowfmt2,'(A,I5,A)') '(',pnum,'(1X,ES17.10))'
       write(80,*) i
       write(80,rowfmt2) pbest
       close(80)

    end do

  end subroutine DDS

  ! ------------------------------------------------------------------
  ! Purpose is to generate a neighboring decision variable value for a single
  !  decision variable value being perturbed by the DDS optimization algorithm.
  !  New DV value respects the upper and lower DV bounds.
  !  Coded by Bryan Tolson, Nov 2005.

  ! I/O variable definitions:
  !  x_cur     - current decision variable (DV) value
  !  x_min     - min DV value
  !  x_max     - max DV value
  !  r         - the neighborhood perturbation factor
  !  new_value - new DV variable value (within specified min and max)

  subroutine neigh_value(x_cur, x_min, x_max, r, new_value)

    use nrtype
    use mo_xor4096, only: xor4096g

    implicit none

    real(dp), intent(in)   :: x_cur, x_min, x_max, r
    real(dp), intent(out)  :: new_value
    real(dp)               :: x_range
    real(dp) :: zvalue

    x_range = x_max - x_min

    ! generate a standard normal random variate (zvalue)
    call xor4096g(0_i8b,zvalue)

    ! calculate new decision variable value:
    new_value = x_cur + zvalue*r*x_range

    !  check new value is within bounds. If not, bounds are reflecting.
    if (new_value < x_min) then
       new_value = x_min + (x_min - new_value)
       if (new_value > x_max) then
          ! if reflection goes past x_max then value should be x_min since
          ! without reflection the approach goes way past lower bound.
          ! This keeps x close to lower bound when x_cur is close to lower bound
          ! Practically speaking, this should never happen with r values <0.3.
          new_value = x_min
       end if
    else if (new_value > x_max) then
       new_value = x_max - (new_value - x_max)
       if (new_value < x_min) then
          ! if reflection goes past x_min then value should be x_max for same reasons as above.
          new_value = x_max
       endif
    endif

  end subroutine neigh_value

  ! ------------------------------------------------------------------

end module mo_dds
