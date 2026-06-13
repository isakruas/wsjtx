!> @file tern_channel.f90
!! @brief HF channel emulation: ITU-R F.1487 Watterson model, AWGN referenced
!!        to 2500 Hz, Middleton Class A impulsive noise, CFO, and clock error.
!!
!! @details
!! Watterson model implementation. The channel is a two-tap delay line. Each
!! tap gain is an independent zero-mean circular complex Gaussian process
!! (Rayleigh envelope) whose power spectrum is Gaussian with standard
!! deviation sigma_f = spread/2, where spread is the ITU frequency-spread
!! parameter (2 sigma). The fading processes are produced by FIR-filtering
!! complex white Gaussian noise at a decimated rate with the spectral square
!! root of the target spectrum, i.e. an impulse response
!!
!!     h(t) proportional to exp(-4 pi^2 sigma_f^2 t^2),
!!
!! then cubic-interpolated to the signal rate. The second tap is delayed by
!! the profile's differential delay using a windowed-sinc fractional-delay
!! filter. Tap powers are 1/2 each so the average channel power gain is
!! unity, leaving the SNR definition unambiguous.
!!
!! All routines accept the sample rate as an argument, so the same channel
!! instrument applies both to the 200 Hz TERN baseband and to full-rate
!! analytic audio signals when calibrating external modems.
!!
!! Statistical compliance (Rayleigh envelope, Gaussian Doppler spectrum width,
!! tap independence, unit power gain) is asserted by test/test_watterson.f90;
!! see docs/VALIDATION.md for the validation methodology.
module tern_channel
  use tern_kinds, only: dp
  use tern_params, only: two_pi, pi, eps, snr_ref_bw_hz, itu_profile
  use tern_rng, only: rng_state, rng_normal, rng_cnormal, rng_uniform
  implicit none
  private
  public :: watterson_fading, watterson_apply, apply_awgn, apply_cfo
  public :: apply_clock_error, apply_impulsive_noise, pad_signal, signal_power

contains

  !> Mean complex-sample power.
  pure function signal_power(x) result(p)
    complex(dp), intent(in) :: x(:)
    real(dp) :: p
    if (size(x) == 0) then
      p = 0.0_dp
    else
      p = sum(abs(x)**2) / real(size(x), dp)
    end if
  end function signal_power

  !> Generate a unit-power complex Gaussian fading process with Gaussian
  !! Doppler spectrum of standard deviation sigma_f = spread_hz / 2.
  !!
  !! @param[in]  spread_hz Frequency spread (2 sigma) in hertz.
  !! @param[in]  n         Number of output samples at rate fs.
  !! @param[in]  fs        Output sample rate in hertz.
  !! @param[inout] rng     Generator state.
  !! @param[out] g         Allocated fading samples, E|g|^2 = 1.
  subroutine watterson_fading(spread_hz, n, fs, rng, g)
    real(dp), intent(in) :: spread_hz, fs
    integer, intent(in) :: n
    type(rng_state), intent(inout) :: rng
    complex(dp), allocatable, intent(out) :: g(:)

    real(dp) :: sigma_f, sigma_t, fs_fade, dt_fade, t, ksum
    real(dp), allocatable :: h(:)
    complex(dp), allocatable :: w(:), gf(:)
    integer :: hw, n_fade, i, k, i0
    real(dp) :: x, frac
    complex(dp) :: p0, p1, p2, p3

    sigma_f = max(spread_hz, 1.0e-3_dp) * 0.5_dp
    ! Time-domain standard deviation of the filter impulse response.
    sigma_t = 1.0_dp / (2.0_dp * sqrt(2.0_dp) * pi * sigma_f)

    ! Fading process rate: comfortably above the Doppler support.
    fs_fade = max(16.0_dp * spread_hz, 1.0_dp)
    fs_fade = min(fs_fade, fs)
    dt_fade = 1.0_dp / fs_fade

    hw = max(2, nint(4.0_dp * sigma_t / dt_fade))
    allocate(h(-hw:hw))
    ksum = 0.0_dp
    do k = -hw, hw
      t = real(k, dp) * dt_fade
      h(k) = exp(-4.0_dp * pi * pi * sigma_f * sigma_f * t * t)
      ksum = ksum + h(k) * h(k)
    end do
    h = h / sqrt(ksum)                 ! Unit noise gain: output variance 1.

    n_fade = int(real(n, dp) * fs_fade / fs) + 2 * hw + 8
    allocate(w(n_fade + 2 * hw), gf(n_fade))
    do i = 1, size(w)
      w(i) = rng_cnormal(rng)
    end do
    do i = 1, n_fade
      gf(i) = (0.0_dp, 0.0_dp)
      do k = -hw, hw
        gf(i) = gf(i) + h(k) * w(i + hw + k)
      end do
    end do

    ! Catmull-Rom interpolation to the signal rate.
    allocate(g(n))
    do i = 1, n
      x = real(i - 1, dp) * fs_fade / fs
      i0 = int(x) + 1
      frac = x - real(i0 - 1, dp)
      i0 = min(max(i0, 2), n_fade - 2)
      p0 = gf(i0 - 1); p1 = gf(i0); p2 = gf(i0 + 1); p3 = gf(i0 + 2)
      g(i) = p1 + 0.5_dp * frac * (p2 - p0 + frac * (2.0_dp * p0 - 5.0_dp * p1 + &
             4.0_dp * p2 - p3 + frac * (3.0_dp * (p1 - p2) + p3 - p0)))
    end do
    deallocate(h, w, gf)
  end subroutine watterson_fading

  !> Apply the two-path Watterson channel. Average power gain is unity.
  subroutine watterson_apply(profile, x, fs, rng, y)
    type(itu_profile), intent(in) :: profile
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: fs
    type(rng_state), intent(inout) :: rng
    complex(dp), allocatable, intent(out) :: y(:)

    complex(dp), allocatable :: g1(:), g2(:), xd(:)
    real(dp) :: inv_sqrt2
    integer :: i, n

    n = size(x)
    if (.not. profile%active) then
      allocate(y(n))
      y = x
      return
    end if

    call watterson_fading(profile%spread_hz, n, fs, rng, g1)
    call watterson_fading(profile%spread_hz, n, fs, rng, g2)
    call fractional_delay(x, profile%delay_ms * 1.0e-3_dp * fs, xd)

    inv_sqrt2 = 1.0_dp / sqrt(2.0_dp)
    allocate(y(n))
    do i = 1, n
      y(i) = inv_sqrt2 * (g1(i) * x(i) + g2(i) * xd(i))
    end do
    deallocate(g1, g2, xd)
  end subroutine watterson_apply

  !> Windowed-sinc fractional delay (Hann window, 65 taps).
  subroutine fractional_delay(x, delay_samples, y)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: delay_samples
    complex(dp), allocatable, intent(out) :: y(:)
    integer, parameter :: hw = 32
    real(dp) :: c, arg, win
    integer :: i, k, n, idx

    n = size(x)
    allocate(y(n))
    y = (0.0_dp, 0.0_dp)
    do i = 1, n
      do k = -hw, hw
        idx = i - k
        if (idx < 1 .or. idx > n) cycle
        arg = real(k, dp) - delay_samples
        if (abs(arg) < eps) then
          c = 1.0_dp
        else
          c = sin(pi * arg) / (pi * arg)
        end if
        win = 0.5_dp + 0.5_dp * cos(pi * arg / real(hw + 1, dp))
        if (abs(arg) > real(hw, dp)) win = 0.0_dp
        y(i) = y(i) + c * win * x(idx)
      end do
    end do
  end subroutine fractional_delay

  !> Additive white Gaussian noise at a target SNR referenced to 2500 Hz.
  !!
  !! With signal power P and N0 = P / (snr * 2500), the per-sample complex
  !! noise variance is N0 * fs. ref_power allows the caller to fix the SNR to
  !! the active-frame power when the buffer contains guard padding, or to the
  !! real-signal power when x is an analytic representation.
  subroutine apply_awgn(x, snr_db_2500, fs, rng, y, noise_var, ref_power)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: snr_db_2500, fs
    type(rng_state), intent(inout) :: rng
    complex(dp), allocatable, intent(out) :: y(:)
    real(dp), intent(out) :: noise_var
    real(dp), intent(in), optional :: ref_power
    real(dp) :: p, n0, sigma
    integer :: i

    if (present(ref_power)) then
      p = max(ref_power, eps)
    else
      p = max(signal_power(x), eps)
    end if
    n0 = p / (max(10.0_dp**(snr_db_2500 / 10.0_dp), eps) * snr_ref_bw_hz)
    noise_var = n0 * fs
    sigma = sqrt(noise_var)
    allocate(y(size(x)))
    do i = 1, size(x)
      y(i) = x(i) + sigma * rng_cnormal(rng)
    end do
  end subroutine apply_awgn

  !> Middleton Class A impulsive noise added on top of an existing buffer.
  !!
  !! Per-sample variance is gauss_var * (m/A + Gp) / (1 + Gp) with m Poisson
  !! distributed with mean A. Gp is the Gaussian-to-impulsive power ratio.
  !! The process reduces to pure AWGN as A grows large.
  subroutine apply_impulsive_noise(x, gauss_var, a_overlap, gamma_ratio, rng, y)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: gauss_var, a_overlap, gamma_ratio
    type(rng_state), intent(inout) :: rng
    complex(dp), allocatable, intent(out) :: y(:)
    real(dp) :: u, p_acc, p_m, var_m, lam
    integer :: i, m

    allocate(y(size(x)))
    lam = max(a_overlap, 1.0e-6_dp)
    do i = 1, size(x)
      u = rng_uniform(rng)
      m = 0
      p_m = exp(-lam)
      p_acc = p_m
      do while (u > p_acc .and. m < 50)
        m = m + 1
        p_m = p_m * lam / real(m, dp)
        p_acc = p_acc + p_m
      end do
      var_m = gauss_var * (real(m, dp) / lam + gamma_ratio) / (1.0_dp + gamma_ratio)
      y(i) = x(i) + sqrt(var_m) * rng_cnormal(rng)
    end do
  end subroutine apply_impulsive_noise

  !> Constant carrier-frequency offset.
  subroutine apply_cfo(x, cfo_hz, fs, y)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: cfo_hz, fs
    complex(dp), allocatable, intent(out) :: y(:)
    real(dp) :: ph
    integer :: i
    allocate(y(size(x)))
    do i = 1, size(x)
      ph = two_pi * cfo_hz * real(i - 1, dp) / fs
      y(i) = x(i) * cmplx(cos(ph), sin(ph), dp)
    end do
  end subroutine apply_cfo

  !> Sample-clock error by Catmull-Rom resampling.
  subroutine apply_clock_error(x, ppm, y)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: ppm
    complex(dp), allocatable, intent(out) :: y(:)
    real(dp) :: scale, src, frac
    integer :: n_in, n_out, i, i0
    complex(dp) :: p0, p1, p2, p3

    n_in = size(x)
    if (abs(ppm) < 1.0e-9_dp) then
      allocate(y(n_in))
      y = x
      return
    end if
    scale = 1.0_dp + ppm * 1.0e-6_dp
    n_out = int(real(n_in - 4, dp) / scale)
    allocate(y(n_out))
    do i = 1, n_out
      src = real(i - 1, dp) * scale
      i0 = int(src) + 1
      frac = src - real(i0 - 1, dp)
      i0 = min(max(i0, 2), n_in - 2)
      p0 = x(i0 - 1); p1 = x(i0); p2 = x(i0 + 1); p3 = x(i0 + 2)
      y(i) = p1 + 0.5_dp * frac * (p2 - p0 + frac * (2.0_dp * p0 - 5.0_dp * p1 + &
             4.0_dp * p2 - p3 + frac * (3.0_dp * (p1 - p2) + p3 - p0)))
    end do
  end subroutine apply_clock_error

  !> Zero-pad before and after the frame.
  subroutine pad_signal(x, prefix_s, suffix_s, fs, y)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: prefix_s, suffix_s, fs
    complex(dp), allocatable, intent(out) :: y(:)
    integer :: pre, post, i
    pre  = max(0, nint(prefix_s * fs))
    post = max(0, nint(suffix_s * fs))
    allocate(y(pre + size(x) + post))
    y = (0.0_dp, 0.0_dp)
    do i = 1, size(x)
      y(pre + i) = x(i)
    end do
  end subroutine pad_signal
end module tern_channel
