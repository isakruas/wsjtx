!> @file tern_frame.f90
!! @brief Frame assembly: codeword-to-symbol mapping and pilot insertion.
!!
!! @details
!! The 256 coded bits fill a 4 x 64 block interleaver: bit position b of data
!! symbol d carries coded bit (b-1)*64 + d. Consecutive coded bits therefore
!! land on consecutive symbols, so a fade of several symbol durations damages
!! bits that are far apart in the codeword, preserving the time diversity the
!! single long codeword is designed to harvest.
module tern_frame
  use tern_kinds, only: dp
  use tern_params, only: n_code, n_symbols, n_data_symbols, bits_per_symbol, &
                         m_tones, k_msg, k_info, pilot_positions, &
                         data_positions, pilot_tone, n_pilot_symbols
  use tern_bits, only: append_crc
  use tern_polar, only: polar_encode
  implicit none
  private
  public :: encode_frame_tones, coded_bit_index, tones_to_llr_order

contains

  !> Coded-bit index carried by bit position b (1..4, MSB first) of data
  !! symbol d (1..64).
  pure integer function coded_bit_index(d, b)
    integer, intent(in) :: d, b
    coded_bit_index = (b - 1) * n_data_symbols + d
  end function coded_bit_index

  !> Encode a 77-bit message into the 85-symbol tone sequence.
  subroutine encode_frame_tones(info_mask, msg, tones)
    integer, intent(in)  :: info_mask(:)
    integer, intent(in)  :: msg(k_msg)
    integer, intent(out) :: tones(n_symbols)
    integer :: info(k_info), coded(n_code)
    integer :: ppos(n_pilot_symbols), dpos(n_data_symbols)
    integer :: d, b, t, k

    call append_crc(msg, info)
    call polar_encode(info_mask, info, coded)

    ppos = pilot_positions()
    dpos = data_positions()

    tones = 0
    do k = 1, n_pilot_symbols
      tones(ppos(k)) = pilot_tone(k)
    end do
    do d = 1, n_data_symbols
      t = 0
      do b = 1, bits_per_symbol
        t = ior(ishft(t, 1), iand(coded(coded_bit_index(d, b)), 1))
      end do
      tones(dpos(d)) = t
    end do
  end subroutine encode_frame_tones

  !> Scatter per-(symbol, bit) LLRs into codeword order for the polar decoder.
  subroutine tones_to_llr_order(llr_sym_bit, llr_coded)
    real(dp), intent(in)  :: llr_sym_bit(:, :)   !< (bits_per_symbol, n_data_symbols)
    real(dp), intent(out) :: llr_coded(n_code)
    integer :: d, b
    do d = 1, n_data_symbols
      do b = 1, bits_per_symbol
        llr_coded(coded_bit_index(d, b)) = llr_sym_bit(b, d)
      end do
    end do
  end subroutine tones_to_llr_order
end module tern_frame
