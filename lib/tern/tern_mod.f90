!> @file tern_mod.f90
!! @brief Continuous-phase frequency-shift keying waveform synthesis.
!!
!! @details
!! The transmitted signal is constant-envelope CPFSK over the orthogonal tone
!! grid. The per-sample instantaneous frequency trajectory is smoothed by a
!! Gaussian pulse (BT = 2) before phase integration, bounding the occupied
!! bandwidth without measurably perturbing the integer-cycle phase property
!! of the tone grid. Constant envelope keeps the waveform insensitive to PA
!! nonlinearity and ALC behavior on SSB transmitters.
module tern_mod
  use tern_kinds, only: dp
  use tern_params, only: tern_mode, two_pi, fs_hz, n_symbols
  implicit none
  private
  public :: synthesize_tones

contains

  !> Synthesize the complex baseband frame for a tone sequence.
  subroutine synthesize_tones(mode, tones, iq)
    type(tern_mode), intent(in) :: mode
    integer, intent(in) :: tones(:)
    complex(dp), allocatable, intent(out) :: iq(:)

    real(dp), allocatable :: f_inst(:), f_smooth(:), kernel(:)
    real(dp) :: sigma, ksum, phase
    integer :: n, s, k, i, j, hw, idx

    n = size(tones) * mode%ns
    allocate(f_inst(n), f_smooth(n))
    do s = 1, size(tones)
      f_inst((s - 1) * mode%ns + 1 : s * mode%ns) = mode%tone_freq(tones(s))
    end do

    ! Gaussian frequency pulse, BT = 2: sigma = sqrt(ln 2) * Ns / (4 pi).
    sigma = sqrt(log(2.0_dp)) * real(mode%ns, dp) / (4.0_dp * acos(-1.0_dp))
    hw = max(1, nint(4.0_dp * sigma))
    allocate(kernel(-hw:hw))
    ksum = 0.0_dp
    do k = -hw, hw
      kernel(k) = exp(-0.5_dp * (real(k, dp) / sigma)**2)
      ksum = ksum + kernel(k)
    end do
    kernel = kernel / ksum

    do i = 1, n
      f_smooth(i) = 0.0_dp
      do k = -hw, hw
        idx = min(max(i + k, 1), n)
        f_smooth(i) = f_smooth(i) + kernel(k) * f_inst(idx)
      end do
    end do

    allocate(iq(n))
    phase = 0.0_dp
    do i = 1, n
      iq(i) = cmplx(cos(phase), sin(phase), dp)
      phase = phase + two_pi * f_smooth(i) / fs_hz
      if (phase > two_pi)  phase = phase - two_pi
      if (phase < -two_pi) phase = phase + two_pi
    end do

    deallocate(f_inst, f_smooth, kernel)
  end subroutine synthesize_tones
end module tern_mod
