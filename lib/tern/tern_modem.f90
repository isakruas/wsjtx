!> @file tern_modem.f90
!! @brief Modem facade and passband scanner, vendored from the TERN reference
!!        implementation (src/api/tern_api.f90), trimmed to the receive path.
!!
!! @details
!! This is the validated tern_modem type and scan_audio procedure, kept
!! verbatim so WSJT-X decodes through exactly the algorithm the TERN campaign
!! qualified. The Monte-Carlo-only entry points (simulate_frame,
!! simulate_repeats, random_message) and their tern_channel/tern_rng
!! dependencies are omitted. The frozen polar set comes from tern_glue so it is
!! computed once for both transmit and receive. scan_hit additionally carries
!! t0_s so the caller can report a WSJT-X time offset.
module tern_modem
  use tern_kinds,  only: dp
  use tern_params, only: tern_mode, make_mode, k_msg
  use tern_frame,  only: encode_frame_tones
  use tern_mod,    only: synthesize_tones
  use tern_demod,  only: demod_opts, tern_result, demodulate
  use tern_audio,  only: audio_to_baseband
  use tern_glue,   only: tern_info_mask
  implicit none
  private
  public :: tern_rx_modem, scan_hit, scan_audio

  type :: tern_rx_modem
    type(tern_mode) :: mode
    integer, allocatable :: info_mask(:)
  contains
    procedure :: init => modem_init
    procedure :: demodulate => modem_demodulate
  end type tern_rx_modem

  type :: scan_hit
    integer  :: msg_bits(k_msg) = 0
    real(dp) :: audio_hz = 0.0_dp
    real(dp) :: es_n0_db = 0.0_dp
    real(dp) :: sync_score = 0.0_dp
    real(dp) :: t0_s = 0.0_dp        !< Frame-start estimate (for WSJT-X dt).
  end type scan_hit

contains

  subroutine modem_init(self, mode_id)
    class(tern_rx_modem), intent(out) :: self
    character(len=1), intent(in) :: mode_id
    self%mode = make_mode(mode_id)
    call tern_info_mask(self%info_mask)
  end subroutine modem_init

  subroutine modem_demodulate(self, rx, opts, res)
    class(tern_rx_modem), intent(in) :: self
    complex(dp), intent(in) :: rx(:)
    type(demod_opts), intent(in) :: opts
    type(tern_result), intent(out) :: res
    call demodulate(self%mode, self%info_mask, rx, opts, res)
  end subroutine modem_demodulate

  !> Decode every TERN signal in an audio passband capture. Verbatim from
  !! tern_api::scan_audio (validated), with t0_s carried into scan_hit.
  subroutine scan_audio(modem, pcm, opts, min_hz, max_hz, hits)
    type(tern_rx_modem), intent(in) :: modem
    real(dp), intent(in) :: pcm(:)
    type(demod_opts), intent(in) :: opts
    real(dp), intent(in) :: min_hz, max_hz
    type(scan_hit), allocatable, intent(out) :: hits(:)

    real(dp), parameter :: chunk_step_hz = 100.0_dp
    real(dp), parameter :: chunk_search_hz = 55.0_dp
    type(scan_hit), allocatable :: found(:)
    type(demod_opts) :: chunk_opts
    type(tern_result) :: res
    complex(dp), allocatable :: x(:)
    real(dp) :: center
    integer :: n_found, i
    logical :: duplicate

    chunk_opts = opts
    chunk_opts%sync_search = .true.
    chunk_opts%f_search_hz = chunk_search_hz

    allocate(found(int((max_hz - min_hz) / chunk_step_hz) + 2))
    n_found = 0
    center = min_hz + 0.5_dp * chunk_step_hz
    do while (center <= max_hz - 0.5_dp * chunk_step_hz + 1.0e-9_dp)
      call audio_to_baseband(pcm, center, x)
      call modem%demodulate(x, chunk_opts, res)
      deallocate(x)
      if (res%crc_ok) then
        duplicate = .false.
        do i = 1, n_found
          if (all(found(i)%msg_bits == res%msg_bits)) then
            duplicate = .true.
            if (res%sync_score > found(i)%sync_score) then
              found(i)%audio_hz = center + res%f0_hz
              found(i)%es_n0_db = res%es_n0_db
              found(i)%sync_score = res%sync_score
              found(i)%t0_s = res%t0_s
            end if
            exit
          end if
        end do
        if (.not. duplicate) then
          n_found = n_found + 1
          found(n_found)%msg_bits = res%msg_bits
          found(n_found)%audio_hz = center + res%f0_hz
          found(n_found)%es_n0_db = res%es_n0_db
          found(n_found)%sync_score = res%sync_score
          found(n_found)%t0_s = res%t0_s
        end if
      end if
      center = center + chunk_step_hz
    end do

    allocate(hits(n_found))
    hits = found(1:n_found)
    deallocate(found)
  end subroutine scan_audio

end module tern_modem
