!> @file tern_kinds.f90
!! @brief Numeric kind parameters shared by all TERN modules.
module tern_kinds
  use iso_fortran_env, only: real64, int32, int64
  implicit none
  private
  public :: dp, i32, i64

  integer, parameter :: dp  = real64
  integer, parameter :: i32 = int32
  integer, parameter :: i64 = int64
end module tern_kinds
