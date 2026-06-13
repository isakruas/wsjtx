!> @file tern_rng.f90
!! @brief Deterministic random number generation for reproducible experiments.
!!
!! @details
!! xorshift64* core with a splitmix64 seeding stage. Every Monte Carlo result
!! in this project is reproducible from (code revision, seed); compiler and
!! platform independence of the stream is part of the experimental contract.
!!
!! Signed-integer overflow is undefined behavior in Fortran, so the 64-bit
!! modular products required by splitmix64 and the xorshift64* output
!! scrambler are evaluated by mul64, which decomposes the operands into
!! 16-bit limbs and reassembles the low 64 bits of the product with bit-model
!! operations only.
module tern_rng
  use tern_kinds, only: dp, i64
  implicit none
  private
  public :: rng_state, rng_init, rng_uniform, rng_normal, rng_cnormal

  type :: rng_state
    integer(i64) :: s = 88172645463325252_i64
  end type rng_state

  integer(i64), parameter :: mask16 = 65535_i64

contains

  !> Low 64 bits of a*b without signed overflow.
  pure function mul64(a, b) result(c)
    integer(i64), intent(in) :: a, b
    integer(i64) :: c
    integer(i64) :: al(4), bl(4), col, carry, limb(4)
    integer :: i, j, n
    do i = 1, 4
      al(i) = iand(ishft(a, -16 * (i - 1)), mask16)
      bl(i) = iand(ishft(b, -16 * (i - 1)), mask16)
    end do
    carry = 0_i64
    do n = 1, 4
      col = carry
      do i = 1, n
        j = n + 1 - i
        col = col + al(i) * bl(j)
      end do
      limb(n) = iand(col, mask16)
      carry = ishft(col, -16)
    end do
    c = ior(ior(limb(1), ishft(limb(2), 16)), &
            ior(ishft(limb(3), 32), ishft(limb(4), 48)))
  end function mul64

  !> Seed the generator. Distinct seeds yield decorrelated streams via a
  !! splitmix64 scrambling pass.
  subroutine rng_init(rng, seed)
    type(rng_state), intent(out) :: rng
    integer, intent(in) :: seed
    integer(i64) :: z
    z = int(seed, i64) + 7046029254386353131_i64
    z = mul64(ieor(z, ishft(z, -30)), -4658895280553007687_i64)
    z = mul64(ieor(z, ishft(z, -27)), -7723592293110705685_i64)
    z = ieor(z, ishft(z, -31))
    if (z == 0_i64) z = 88172645463325252_i64
    rng%s = z
  end subroutine rng_init

  !> Uniform deviate in (0, 1).
  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(i64) :: x
    x = rng%s
    x = ieor(x, ishft(x, 13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x, 17))
    rng%s = x
    x = mul64(x, 2685821657736338717_i64)
    u = (real(ishft(x, -11), dp) + 0.5_dp) * 2.0_dp**(-53)
    if (u <= 0.0_dp) u = 2.0_dp**(-53)
    if (u >= 1.0_dp) u = 1.0_dp - 2.0_dp**(-53)
  end function rng_uniform

  !> Standard normal deviate (Box-Muller).
  function rng_normal(rng) result(g)
    type(rng_state), intent(inout) :: rng
    real(dp) :: g
    real(dp) :: u1, u2
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    g = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
  end function rng_normal

  !> Circularly symmetric complex normal deviate with unit total variance.
  function rng_cnormal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    complex(dp) :: z
    real(dp) :: u1, u2, r
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    r = sqrt(-log(u1))
    z = cmplx(r * cos(2.0_dp * acos(-1.0_dp) * u2), &
              r * sin(2.0_dp * acos(-1.0_dp) * u2), dp)
  end function rng_cnormal
end module tern_rng
