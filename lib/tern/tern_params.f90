!> @file tern_params.f90
!! @brief Mode geometry, frame structure, code parameters, and ITU-R F.1487
!!        channel profiles.
!!
!! @details
!! TERN transmits K_MSG = 77 message bits protected by a 14-bit CRC inside a
!! single polar codeword of length N_CODE = 256. The codeword is mapped onto
!! N_DATA = 64 symbols of an M_TONES = 16 orthogonal alphabet (4 coded bits
!! per symbol). Three length-7 Costas arrays (start, middle, end of frame)
!! provide time-frequency synchronization and act as coherent anchors for the
!! Gauss-Markov channel tracker. Tone spacing equals 1/T_s, so the alphabet is
!! orthogonal over one symbol and, because every tone frequency is an integer
!! multiple of 1/T_s, continuous-phase keying accumulates an integer number of
!! cycles per symbol: the channel phase observed by the demodulator is
!! data-independent modulo 2*pi.
!!
!! Two operating modes share the frame structure and differ only in symbol
!! duration (hence bandwidth, airtime, and threshold):
!!
!!   Mode A: T_s = 0.32 s, tone spacing 3.125 Hz, occupied width 50 Hz,
!!           frame 27.2 s (30 s slot). Measured AWGN threshold -22.7 dB.
!!   Mode B: T_s = 0.64 s, tone spacing 1.5625 Hz, occupied width 25 Hz,
!!           frame 54.4 s (60 s slot). Measured AWGN threshold -25.4 dB.
!!
!! All SNR figures in this project are referenced to a 2500 Hz noise
!! bandwidth (the WSJT-X convention) so thresholds are directly comparable
!! with published FT8/FST4 numbers.
module tern_params
  use tern_kinds, only: dp
  implicit none
  private
  public :: pi, two_pi, eps
  public :: fs_hz, m_tones, n_symbols, n_data_symbols, n_pilot_symbols
  public :: bits_per_symbol, k_msg, k_crc, k_info, n_code, scl_list_size
  public :: crc14_poly, snr_ref_bw_hz, costas7
  public :: pilot_positions, data_positions, pilot_tone
  public :: tern_mode, make_mode
  public :: itu_profile, get_itu_profile, itu_profile_names

  real(dp), parameter :: pi     = acos(-1.0_dp)
  real(dp), parameter :: two_pi = 2.0_dp * pi
  real(dp), parameter :: eps    = 1.0e-12_dp

  !> Internal complex-baseband sample rate. Narrowband by construction: at
  !! 200 Hz a +-50 ppm sound-card clock error accumulates < 0.3 sample over a
  !! full Mode A frame, eliminating the need for timing tracking loops.
  real(dp), parameter :: fs_hz = 200.0_dp

  integer, parameter :: m_tones         = 16   !< Orthogonal alphabet size.
  integer, parameter :: bits_per_symbol = 4    !< log2(m_tones).
  integer, parameter :: n_symbols       = 85   !< Symbols per frame.
  integer, parameter :: n_data_symbols  = 64   !< Payload-bearing symbols.
  integer, parameter :: n_pilot_symbols = 21   !< Three Costas-7 arrays.

  integer, parameter :: k_msg  = 77            !< Message bits per frame.
  integer, parameter :: k_crc  = 14            !< CRC bits appended to message.
  integer, parameter :: k_info = k_msg + k_crc !< Polar information length.
  integer, parameter :: n_code = 256           !< Polar codeword length.
  integer, parameter :: scl_list_size = 8      !< SCL decoder list size.

  !> CRC-14 generator polynomial (x^14 implicit), MSB-first bitwise update.
  integer, parameter :: crc14_poly = int(z'2757')

  !> Reference noise bandwidth for all reported SNR values.
  real(dp), parameter :: snr_ref_bw_hz = 2500.0_dp

  !> Order-7 Costas array. Cyclic auto-ambiguity sidelobes are at most 1 hit,
  !! giving an unambiguous time-frequency acquisition surface.
  integer, parameter :: costas7(7) = [3, 1, 4, 0, 6, 5, 2]

  !> Mode-dependent waveform geometry.
  type :: tern_mode
    character(len=1) :: id = 'A'
    integer  :: ns = 64                !< Samples per symbol at fs_hz.
    real(dp) :: t_sym = 0.32_dp        !< Symbol duration in seconds.
    real(dp) :: df = 3.125_dp          !< Tone spacing = 1/t_sym in hertz.
  contains
    procedure :: frame_samples
    procedure :: tone_freq
    procedure :: occupied_bw
  end type tern_mode

  !> ITU-R F.1487 Watterson profile: two independent Rayleigh paths with a
  !! Gaussian Doppler spectrum. Frequency spread is defined as twice the
  !! standard deviation of the Doppler spectrum (2*sigma).
  type :: itu_profile
    character(len=24) :: name = "none"
    real(dp) :: delay_ms  = 0.0_dp     !< Differential path delay.
    real(dp) :: spread_hz = 0.0_dp     !< Frequency spread (2*sigma).
    logical  :: active    = .false.
  end type itu_profile

contains

  !> Construct the waveform geometry for mode 'A' or 'B'.
  function make_mode(id) result(m)
    character(len=1), intent(in) :: id
    type(tern_mode) :: m
    select case (id)
    case ('A', 'a')
      m%id = 'A'; m%ns = 64
    case ('B', 'b')
      m%id = 'B'; m%ns = 128
    case default
      error stop "make_mode: mode must be 'A' or 'B'"
    end select
    m%t_sym = real(m%ns, dp) / fs_hz
    m%df    = fs_hz / real(m%ns, dp)
  end function make_mode

  pure integer function frame_samples(self)
    class(tern_mode), intent(in) :: self
    frame_samples = n_symbols * self%ns
  end function frame_samples

  !> Baseband frequency of tone index m in [0, m_tones-1]. Integer multiples
  !! of the spacing keep continuous-phase symbols at integer cycle counts.
  pure real(dp) function tone_freq(self, m)
    class(tern_mode), intent(in) :: self
    integer, intent(in) :: m
    tone_freq = real(m - m_tones / 2, dp) * self%df
  end function tone_freq

  pure real(dp) function occupied_bw(self)
    class(tern_mode), intent(in) :: self
    occupied_bw = real(m_tones, dp) * self%df
  end function occupied_bw

  !> One-based symbol positions of the three Costas pilot blocks.
  pure function pilot_positions() result(pos)
    integer :: pos(n_pilot_symbols)
    integer :: i
    do i = 1, 7
      pos(i)      = i
      pos(7 + i)  = 39 + i
      pos(14 + i) = 78 + i
    end do
  end function pilot_positions

  !> One-based symbol positions of the 64 data symbols.
  pure function data_positions() result(pos)
    integer :: pos(n_data_symbols)
    integer :: i
    do i = 1, 32
      pos(i)      = 7 + i
      pos(32 + i) = 46 + i
    end do
  end function data_positions

  !> Tone index used by the k-th pilot symbol (k in 1..21). The Costas values
  !! are dilated across the full 16-tone grid to sharpen the frequency
  !! discrimination of the acquisition metric.
  pure integer function pilot_tone(k)
    integer, intent(in) :: k
    pilot_tone = 2 * costas7(mod(k - 1, 7) + 1) + 1
  end function pilot_tone

  !> Resolve a named ITU-R F.1487 profile (Table 1 channel parameters).
  function get_itu_profile(name) result(p)
    character(len=*), intent(in) :: name
    type(itu_profile) :: p
    p%name = name
    p%active = .true.
    select case (trim(name))
    case ("none")
      p%active = .false.
    case ("low_quiet")
      p%delay_ms = 0.5_dp;  p%spread_hz = 0.5_dp
    case ("low_moderate")
      p%delay_ms = 2.0_dp;  p%spread_hz = 1.5_dp
    case ("low_disturbed")
      p%delay_ms = 6.0_dp;  p%spread_hz = 10.0_dp
    case ("mid_quiet")
      p%delay_ms = 0.5_dp;  p%spread_hz = 0.1_dp
    case ("mid_moderate")
      p%delay_ms = 1.0_dp;  p%spread_hz = 0.5_dp
    case ("mid_disturbed")
      p%delay_ms = 2.0_dp;  p%spread_hz = 1.0_dp
    case ("high_quiet")
      p%delay_ms = 1.0_dp;  p%spread_hz = 0.5_dp
    case ("high_moderate")
      p%delay_ms = 3.0_dp;  p%spread_hz = 10.0_dp
    case ("high_disturbed")
      p%delay_ms = 7.0_dp;  p%spread_hz = 30.0_dp
    case default
      error stop "get_itu_profile: unknown profile name"
    end select
  end function get_itu_profile

  !> Comma-separated list of supported profile names for CLI help output.
  pure function itu_profile_names() result(s)
    character(len=160) :: s
    s = "none, low_quiet, low_moderate, low_disturbed, mid_quiet, " // &
        "mid_moderate, mid_disturbed, high_quiet, high_moderate, high_disturbed"
  end function itu_profile_names
end module tern_params
