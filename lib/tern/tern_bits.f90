!> @file tern_bits.f90
!! @brief Bit-vector utilities and the frame CRC.
module tern_bits
  use tern_kinds, only: dp
  use tern_params, only: crc14_poly, k_msg, k_crc
  implicit none
  private
  public :: crc14, append_crc, check_crc, parse_bit_string, bits_to_string

contains

  !> CRC-14 over a bit vector, MSB-first, zero initial register, generator
  !! polynomial tern_params::crc14_poly with implicit x^14 term.
  pure function crc14(bits) result(reg)
    integer, intent(in) :: bits(:)
    integer :: reg
    integer :: i, fb
    reg = 0
    do i = 1, size(bits)
      fb = ieor(iand(bits(i), 1), iand(ishft(reg, -(k_crc - 1)), 1))
      reg = iand(ishft(reg, 1), 2**k_crc - 1)
      if (fb == 1) reg = ieor(reg, crc14_poly)
    end do
  end function crc14

  !> Append the CRC of msg(1:k_msg) producing the k_info polar payload.
  pure subroutine append_crc(msg, info)
    integer, intent(in)  :: msg(k_msg)
    integer, intent(out) :: info(k_msg + k_crc)
    integer :: c, i
    info(1:k_msg) = msg
    c = crc14(msg)
    do i = 1, k_crc
      info(k_msg + i) = iand(ishft(c, -(k_crc - i)), 1)
    end do
  end subroutine append_crc

  !> Verify the CRC of a decoded k_info block.
  pure function check_crc(info) result(ok)
    integer, intent(in) :: info(k_msg + k_crc)
    logical :: ok
    integer :: c, i, ref
    c = crc14(info(1:k_msg))
    ref = 0
    do i = 1, k_crc
      ref = ior(ishft(ref, 1), iand(info(k_msg + i), 1))
    end do
    ok = (c == ref)
  end function check_crc

  !> Parse a '0'/'1' string into a fixed-length bit vector, zero-padded on the
  !! right. Returns .false. on invalid characters or overlength input.
  subroutine parse_bit_string(text, bits, ok)
    character(len=*), intent(in) :: text
    integer, intent(out) :: bits(:)
    logical, intent(out) :: ok
    integer :: i, n
    bits = 0
    n = len_trim(text)
    ok = (n <= size(bits))
    if (.not. ok) return
    do i = 1, n
      select case (text(i:i))
      case ('0')
        bits(i) = 0
      case ('1')
        bits(i) = 1
      case default
        ok = .false.
        return
      end select
    end do
  end subroutine parse_bit_string

  !> Render a bit vector as a '0'/'1' string.
  pure function bits_to_string(bits) result(s)
    integer, intent(in) :: bits(:)
    character(len=size(bits)) :: s
    integer :: i
    do i = 1, size(bits)
      s(i:i) = merge('1', '0', iand(bits(i), 1) == 1)
    end do
  end function bits_to_string
end module tern_bits
