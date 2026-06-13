!> @file tern_glue.f90
!! @brief WSJT-X integration glue for the vendored TERN modem.
!!
!! @details
!! Centralizes the two pieces of modem state that both the transmitter
!! (gentern/gen_ternwave) and the receiver (tern_decode) must share: the
!! frozen polar information set and the mode geometry. The frozen set is
!! constructed once and cached, since polar_info_set runs the
!! Gaussian-approximation density evolution on every call.
!!
!! WSJT-X submode convention used throughout the integration:
!!     nsubmode = 0  ->  TERN Mode A (0.32 s symbols, 50 Hz, 30 s slot)
!!     nsubmode = 1  ->  TERN Mode B (0.64 s symbols, 25 Hz, 60 s slot)
module tern_glue
  use tern_kinds,  only: dp
  use tern_params, only: tern_mode, make_mode
  use tern_polar,  only: polar_info_set
  implicit none
  private
  public :: tern_info_mask, tern_mode_from_n, tern_audio_synth

  !> Polar design point: LLR mean 4*R*Eb/N0 at Eb/N0 = 2 dB, R = 91/256.
  !! Must match src/api/tern_api.f90 so the frozen set is identical to the
  !! one the reference TERN implementation validates against.
  real(dp), parameter :: design_llr_mean = 4.0_dp * (91.0_dp / 256.0_dp) * &
                                           10.0_dp**(0.2_dp)

contains

  !> Cached frozen information set for the (256,91) polar code.
  subroutine tern_info_mask(mask)
    integer, allocatable, intent(out) :: mask(:)
    integer, allocatable, save :: cached(:)
    logical, save :: have = .false.
    if (.not. have) then
      call polar_info_set(design_llr_mean, cached)
      have = .true.
    end if
    mask = cached
  end subroutine tern_info_mask

  !> Map a WSJT-X submode number to a TERN mode geometry.
  function tern_mode_from_n(nsubmode) result(m)
    integer, intent(in) :: nsubmode
    type(tern_mode) :: m
    if (nsubmode == 1) then
      m = make_mode('B')
    else
      m = make_mode('A')
    end if
  end function tern_mode_from_n

  !> Interpolate a 200 Hz complex baseband frame to a real audio passband
  !! waveform at fsample, centered on offset_hz. Faithful generalization of
  !! tern_audio::baseband_to_audio to an arbitrary integer interpolation
  !! factor (240 for 48 kHz), reusing the same 80 Hz Blackman-windowed sinc
  !! anti-imaging kernel scaled to the output rate. Output is peak-normalized
  !! to 0.9 like the reference transmitter.
  subroutine tern_audio_synth(x, offset_hz, fsample, pcm)
    complex(dp), intent(in) :: x(:)
    real(dp),    intent(in) :: offset_hz, fsample
    real(dp), allocatable, intent(out) :: pcm(:)

    real(dp), parameter :: fs_base = 200.0_dp     !< TERN baseband rate.
    real(dp), parameter :: cutoff_hz = 80.0_dp    !< Anti-imaging cutoff.
    real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer  :: decim, fir_half
    real(dp), allocatable :: h(:)
    integer  :: n_in, n_out, m, n, k, lo, hi
    real(dp) :: ph, peak, xv, w, hv
    complex(dp) :: acc

    decim = nint(fsample / fs_base)               ! 240 for 48 kHz
    fir_half = 600 * decim / 60                    ! same time-support as 12 kHz/600
    allocate(h(-fir_half:fir_half))
    do k = -fir_half, fir_half
      if (k == 0) then
        hv = 2.0_dp * cutoff_hz / fsample
      else
        xv = two_pi * cutoff_hz * real(k, dp) / fsample
        hv = sin(xv) / (pi * real(k, dp))
      end if
      w = 0.42_dp + 0.5_dp * cos(pi * real(k, dp) / real(fir_half, dp)) + &
          0.08_dp * cos(two_pi * real(k, dp) / real(fir_half, dp))
      h(k) = hv * w
    end do

    n_in = size(x)
    n_out = n_in * decim
    allocate(pcm(n_out))
    do m = 1, n_out
      acc = (0.0_dp, 0.0_dp)
      lo = max(1, (m - 1 - fir_half) / decim)
      hi = min(n_in, (m - 1 + fir_half) / decim + 2)
      do n = lo, hi
        k = (m - 1) - (n - 1) * decim
        if (abs(k) > fir_half) cycle
        acc = acc + h(k) * x(n)
      end do
      ph = two_pi * offset_hz * real(m - 1, dp) / fsample
      pcm(m) = real(acc * cmplx(cos(ph), sin(ph), dp), dp) * real(decim, dp)
    end do

    peak = maxval(abs(pcm))
    if (peak > 0.0_dp) pcm = 0.9_dp * pcm / peak
    deallocate(h)
  end subroutine tern_audio_synth

end module tern_glue
