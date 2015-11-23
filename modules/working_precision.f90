!-----------------------------------------------------------------------
!> @file    working_precision.f90
!! @brief   Selects (real and complex) working precision.
!! @details This module specifies a kind parameter value for single,
!!          double or quadruple precision real and complex arithmetics.
!! @author  bgeneto
!! @date    20150623
!-----------------------------------------------------------------------
module working_precision
    use iso_fortran_env, only: wp => real64    ! real32, real64, real128
end module working_precision

