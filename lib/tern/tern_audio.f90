!> @file tern_audio.f90
!! @brief Audio passband interface: WAV file I/O and conversion between the
!!        12 kHz real audio passband and the 200 Hz complex baseband.
!!
!! @details
!! TERN operates on SSB transceivers as a narrow tone group placed at a
!! configurable audio offset (1500 Hz by default) inside the transmit
!! passband. This module converts between that representation and the
!! internal complex baseband:
!!
!!   - Downconversion: the real audio signal is mixed by exp(-i 2 pi f0 t),
!!     low-pass filtered, and decimated by the integer factor 60. The mixing
!!     image at -2 f0 lies far outside the filter passband and is rejected
!!     by the same filter.
!!   - Upconversion: the complex baseband is interpolated by 60 with the
!!     same low-pass kernel, mixed to the audio offset, and the real part
!!     is scaled to 16-bit PCM.
!!
!! The anti-aliasing kernel is a Blackman-windowed sinc with 80 Hz cutoff
!! at the 12 kHz rate, giving better than 70 dB image and alias rejection
!! for the 50 Hz occupied bandwidth.
!!
!! The WAV interface accepts canonical RIFF files: mono, 16-bit PCM,
!! 12000 Hz (the WSJT-X audio convention), with unknown chunks skipped.
module tern_audio
  use tern_kinds, only: dp
  use tern_params, only: fs_hz, two_pi, pi
  implicit none
  private
  public :: audio_rate_hz, decim_factor
  public :: load_wav, save_wav, audio_to_baseband, baseband_to_audio
  public :: analytic_signal

  real(dp), parameter :: audio_rate_hz = 12000.0_dp
  integer,  parameter :: decim_factor  = 60
  integer,  parameter :: fir_half      = 600   !< Kernel half-length in audio samples.
  real(dp), parameter :: fir_cutoff_hz = 80.0_dp

contains

  !> Blackman-windowed sinc kernel value at integer audio-sample offset k,
  !! normalized for unit passband gain at the audio rate.
  pure function fir_tap(k) result(h)
    integer, intent(in) :: k
    real(dp) :: h
    real(dp) :: x, w
    if (abs(k) > fir_half) then
      h = 0.0_dp
      return
    end if
    x = two_pi * fir_cutoff_hz * real(k, dp) / audio_rate_hz
    if (k == 0) then
      h = 2.0_dp * fir_cutoff_hz / audio_rate_hz
    else
      h = sin(x) / (pi * real(k, dp))
    end if
    w = 0.42_dp + 0.5_dp * cos(pi * real(k, dp) / real(fir_half, dp)) + &
        0.08_dp * cos(two_pi * real(k, dp) / real(fir_half, dp))
    h = h * w
  end function fir_tap

  !> Convert real audio samples to complex baseband at fs_hz.
  !!
  !! @param[in]  pcm       Audio samples at audio_rate_hz, any real scale.
  !! @param[in]  offset_hz Audio frequency translated to baseband zero.
  !! @param[out] x         Allocated baseband samples at fs_hz.
  subroutine audio_to_baseband(pcm, offset_hz, x)
    real(dp), intent(in) :: pcm(:)
    real(dp), intent(in) :: offset_hz
    complex(dp), allocatable, intent(out) :: x(:)
    complex(dp), allocatable :: mixed(:)
    real(dp) :: h(-fir_half:fir_half)
    integer :: n_out, n_in, m, k, idx
    real(dp) :: ph
    complex(dp) :: acc

    n_in = size(pcm)
    n_out = n_in / decim_factor
    if (n_out < 1) error stop "audio_to_baseband: input shorter than one output sample"

    do k = -fir_half, fir_half
      h(k) = fir_tap(k)
    end do
    allocate(mixed(n_in))
    do idx = 1, n_in
      ph = -two_pi * offset_hz * real(idx - 1, dp) / audio_rate_hz
      mixed(idx) = pcm(idx) * cmplx(cos(ph), sin(ph), dp)
    end do

    allocate(x(n_out))
    do m = 1, n_out
      acc = (0.0_dp, 0.0_dp)
      do k = -fir_half, fir_half
        idx = (m - 1) * decim_factor + 1 - k
        if (idx < 1 .or. idx > n_in) cycle
        acc = acc + h(k) * mixed(idx)
      end do
      ! Factor 2 restores the analytic-signal amplitude lost when mixing a
      ! real passband signal (half the power sits in the rejected image).
      x(m) = 2.0_dp * acc
    end do
    deallocate(mixed)
  end subroutine audio_to_baseband

  !> Convert complex baseband at fs_hz to real audio samples.
  !!
  !! @param[in]  x         Baseband samples at fs_hz.
  !! @param[in]  offset_hz Audio frequency of baseband zero.
  !! @param[out] pcm       Allocated audio samples, peak-normalized to 0.9.
  subroutine baseband_to_audio(x, offset_hz, pcm)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: offset_hz
    real(dp), allocatable, intent(out) :: pcm(:)
    integer :: n_in, n_out, m, n, k
    real(dp) :: ph, peak
    complex(dp) :: acc

    real(dp) :: h(-fir_half:fir_half)

    do k = -fir_half, fir_half
      h(k) = fir_tap(k)
    end do
    n_in = size(x)
    n_out = n_in * decim_factor
    allocate(pcm(n_out))
    do m = 1, n_out
      acc = (0.0_dp, 0.0_dp)
      ! Only baseband samples within the kernel support contribute.
      do n = max(1, (m - 1 - fir_half) / decim_factor), &
             min(n_in, (m - 1 + fir_half) / decim_factor + 2)
        k = (m - 1) - (n - 1) * decim_factor
        if (abs(k) > fir_half) cycle
        acc = acc + h(k) * x(n)
      end do
      ph = two_pi * offset_hz * real(m - 1, dp) / audio_rate_hz
      pcm(m) = real(acc * cmplx(cos(ph), sin(ph), dp), dp) * real(decim_factor, dp)
    end do
    peak = maxval(abs(pcm))
    if (peak > 0.0_dp) pcm = 0.9_dp * pcm / peak
  end subroutine baseband_to_audio

  !> Write mono 16-bit PCM WAV at audio_rate_hz. Input is clipped to [-1, 1].
  subroutine save_wav(path, pcm)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: pcm(:)
    integer :: unit, ios, i
    integer(kind=2), allocatable :: samples(:)
    integer :: n, data_bytes
    real(dp) :: v

    n = size(pcm)
    data_bytes = 2 * n
    allocate(samples(n))
    do i = 1, n
      v = max(-1.0_dp, min(1.0_dp, pcm(i)))
      samples(i) = int(nint(v * 32767.0_dp), kind=2)
    end do

    open(newunit=unit, file=path, access="stream", form="unformatted", &
         status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "save_wav: cannot open output file"
    write(unit) "RIFF", int32_le(36 + data_bytes), "WAVE"
    write(unit) "fmt ", int32_le(16), int16_le(1), int16_le(1)
    write(unit) int32_le(nint(audio_rate_hz)), int32_le(nint(audio_rate_hz) * 2)
    write(unit) int16_le(2), int16_le(16)
    write(unit) "data", int32_le(data_bytes)
    write(unit) samples
    close(unit)
    deallocate(samples)
  end subroutine save_wav

  !> Read mono 16-bit PCM WAV at audio_rate_hz, skipping non-essential chunks.
  subroutine load_wav(path, pcm)
    character(len=*), intent(in) :: path
    real(dp), allocatable, intent(out) :: pcm(:)
    integer :: unit, ios, i
    character(len=4) :: tag
    integer :: chunk_size, fmt_code, n_channels, sample_rate, bits
    integer :: data_bytes, n
    integer(kind=2), allocatable :: samples(:)
    logical :: have_fmt
    integer(kind=8) :: pos

    open(newunit=unit, file=path, access="stream", form="unformatted", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_wav: cannot open input file"
    read(unit) tag
    if (tag /= "RIFF") error stop "load_wav: not a RIFF file"
    read(unit) chunk_size
    read(unit) tag
    if (tag /= "WAVE") error stop "load_wav: not a WAVE file"

    have_fmt = .false.
    data_bytes = -1
    do
      read(unit, iostat=ios) tag, chunk_size
      if (ios /= 0) exit
      select case (tag)
      case ("fmt ")
        block
          integer(kind=2) :: w16(2), w16b(2)
          integer :: w32(2)
          read(unit) w16, w32, w16b
          fmt_code = int(w16(1)); n_channels = int(w16(2))
          sample_rate = w32(1); bits = int(w16b(2))
          if (chunk_size > 16) then
            inquire(unit=unit, pos=pos)
            read(unit, pos=pos + int(chunk_size - 16, kind=8)) ! reposition
          end if
          have_fmt = .true.
        end block
      case ("data")
        data_bytes = chunk_size
        exit
      case default
        inquire(unit=unit, pos=pos)
        if (mod(chunk_size, 2) == 1) chunk_size = chunk_size + 1
        read(unit, pos=pos + int(chunk_size, kind=8), iostat=ios)
        if (ios /= 0) exit
      end select
    end do

    if (.not. have_fmt) error stop "load_wav: missing fmt chunk"
    if (data_bytes < 0) error stop "load_wav: missing data chunk"
    if (fmt_code /= 1) error stop "load_wav: PCM format required"
    if (n_channels /= 1) error stop "load_wav: mono required"
    if (bits /= 16) error stop "load_wav: 16-bit samples required"
    if (sample_rate /= nint(audio_rate_hz)) &
         error stop "load_wav: 12000 Hz sample rate required"

    n = data_bytes / 2
    allocate(samples(n), pcm(n))
    read(unit) samples
    close(unit)
    do i = 1, n
      pcm(i) = real(samples(i), dp) / 32767.0_dp
    end do
    deallocate(samples)
  end subroutine load_wav

  !> Analytic-signal construction by a windowed FIR Hilbert transformer.
  !!
  !! @param[in]  x Real samples at any rate.
  !! @param[out] z Allocated analytic samples z = x + i*H(x). The transformer
  !!              is applied noncausally, centered on each sample, so Re(z)
  !!              equals the input exactly and the quadrature component is
  !!              time-aligned with it, apart from edge transients within
  !!              half a kernel length of the record boundaries.
  !!
  !! The transformer is an odd-length type III FIR (801 taps, Blackman
  !! window): impulse response 2/(pi*k) for odd lags, zero for even lags.
  !! Passband ripple is below 0.1 dB for 100 Hz to 5900 Hz at the 12 kHz
  !! audio rate, covering the entire SSB passband.
  subroutine analytic_signal(x, z)
    real(dp), intent(in) :: x(:)
    complex(dp), allocatable, intent(out) :: z(:)
    integer, parameter :: hh = 400
    real(dp) :: h(-hh:hh), w, q
    integer :: i, k, n, idx

    do k = -hh, hh
      if (mod(abs(k), 2) == 1) then
        w = 0.42_dp + 0.5_dp * cos(pi * real(k, dp) / real(hh + 1, dp)) + &
            0.08_dp * cos(two_pi * real(k, dp) / real(hh + 1, dp))
        h(k) = 2.0_dp / (pi * real(k, dp)) * w
      else
        h(k) = 0.0_dp
      end if
    end do

    n = size(x)
    allocate(z(n))
    do i = 1, n
      q = 0.0_dp
      do k = -hh + 1, hh - 1, 2
        idx = i - k
        if (idx < 1 .or. idx > n) cycle
        q = q + h(k) * x(idx)
      end do
      z(i) = cmplx(x(i), q, dp)
    end do
  end subroutine analytic_signal

  !> 32-bit little-endian integer for stream output.
  pure function int32_le(v) result(b)
    integer, intent(in) :: v
    character(len=4) :: b
    b = achar(iand(v, 255)) // achar(iand(ishft(v, -8), 255)) // &
        achar(iand(ishft(v, -16), 255)) // achar(iand(ishft(v, -24), 255))
  end function int32_le

  !> 16-bit little-endian integer for stream output.
  pure function int16_le(v) result(b)
    integer, intent(in) :: v
    character(len=2) :: b
    b = achar(iand(v, 255)) // achar(iand(ishft(v, -8), 255))
  end function int16_le
end module tern_audio
