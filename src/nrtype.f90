module nrtype
    use, intrinsic :: iso_c_binding,   only: &
             c_short, c_int, c_long, c_float, c_double, c_float_complex, c_double_complex, c_bool

    !> 8byte integer
    !integer, parameter :: i8b = SELECTED_INT_KIND(18)  ! 8byte integer
    integer, parameter :: i8b = c_long 
    !> 2byte integer
    !integer, parameter :: I4B = SELECTED_INT_KIND(9)   
    integer, parameter :: I4B = c_int 
    !> 2byte integer
    !integer, parameter :: I2B = SELECTED_INT_KIND(4)
    integer, parameter :: I2B = c_short 
    !> 1byte integer
    integer, parameter :: I1B = SELECTED_INT_KIND(2)
    !> single precision  real
    !integer, parameter :: SP = KIND(1.0)
    integer, parameter :: SP = c_float 
    !> double precision real
    !integer, parameter :: DP = KIND(1.0D0)
    integer, parameter :: DP = c_double 
    !> single precision  complex 
    !integer, parameter :: SPC = KIND((1.0,1.0))
    integer, parameter :: SPC = c_float_complex
    !> double precision  complex 
    !integer, parameter :: DPC = KIND((1.0D0,1.0D0))
    integer, parameter :: DPC = c_double_complex 
    !> logical
    integer, parameter :: LGT = KIND(.true.)
    ! Real kinds
    integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
    integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real
    ! Integer kinds
    integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
    integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer
    !Complex kinds
    integer, parameter :: kc4 = kr4                            ! single precision complex
    integer, parameter :: kc8 = kr8                            ! double precision complex
    real(SP), parameter :: PI=3.141592653589793238462643383279502884197_sp
    real(SP), parameter :: PIO2=1.57079632679489661923132169163975144209858_sp
    real(SP), parameter :: TWOPI=6.283185307179586476925286766559005768394_sp
    real(SP), parameter :: SQRT2=1.41421356237309504880168872420969807856967_sp
    real(SP), parameter :: EULER=0.5772156649015328606065120900824024310422_sp
    real(DP), parameter :: PI_D=3.141592653589793238462643383279502884197_dp
    real(DP), parameter :: PIO2_D=1.57079632679489661923132169163975144209858_dp
    real(DP), parameter :: TWOPI_D=6.283185307179586476925286766559005768394_dp
    TYPE sprs2_sp
        integer(I4B) :: n,len
        real(SP), dimension(:), pointer :: val
        integer(I4B), dimension(:), pointer :: irow
        integer(I4B), dimension(:), pointer :: jcol
    END TYPE sprs2_sp
    TYPE sprs2_dp
        integer(I4B) :: n,len
        real(DP), dimension(:), pointer :: val
        integer(I4B), dimension(:), pointer :: irow
        integer(I4B), dimension(:), pointer :: jcol
    END TYPE sprs2_dp
end module nrtype
