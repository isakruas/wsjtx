!> @file tern_decode.f90
!! @brief TERN receiver as a WSJT-X decoder class, mirroring fst4_decoder.
!!
!! @details
!! The jt9 dispatcher (lib/decoder.f90) hands us the slot's 12 kHz int16 audio
!! (id2) plus the search band [nfa,nfb]. Decoding is delegated to the validated
!! tern_modem::scan_audio (the same routine the TERN campaign qualified): it
!! scans the band, rejects false candidates by the frame CRC, and merges
!! adjacent-chunk duplicates. Each surviving decode is reported through the
!! callback, whose formatted stdout line is parsed by
!! MainWindow::readFromStdout, exactly as for the native modes.
!!
!! SNR is converted from the receiver's per-symbol Es/N0 estimate to the
!! WSJT-X / TERN 2500 Hz reference: SNR_2500 = Es/N0 + 10*log10(Rs / 2500),
!! with symbol rate Rs = 1 / t_sym. This is the convention under which every
!! threshold in docs/VALIDATION.md is stated.
module tern_decode

  use tern_kinds,  only: dp
  use tern_params, only: k_msg, snr_ref_bw_hz
  use tern_demod,  only: demod_opts, tern_result, demodulate, decode_combined
  use tern_audio,  only: audio_rate_hz, audio_to_baseband
  use tern_modem,  only: tern_rx_modem, scan_hit, scan_audio
  implicit none
  private
  public :: tern_decoder

  ! Cross-frame combining state (TERN repeat protocol). Persisted across slots:
  ! weak non-decoding frames at the Rx frequency have their pilot-only tone
  ! log-likelihoods accumulated, and decode_combined is attempted each slot.
  ! Reset whenever a frame decodes on its own (the QSO moved to a new message)
  ! or the accumulation goes stale, which keeps mismatched messages from being
  ! summed.
  real(dp), allocatable, save :: logp_acc(:,:)
  integer,  save :: nacc = 0
  integer,  save :: acc_utc = -1
  real(dp), save :: acc_freq = -1.0_dp
  integer,  parameter :: max_combine = 8

  type :: tern_decoder
    procedure(tern_decode_callback), pointer :: callback => null()
  contains
    procedure :: decode
  end type tern_decoder

  abstract interface
    subroutine tern_decode_callback(this, nutc, sync, nsnr, dt, freq,    &
        decoded, nap, qual, ntrperiod, fmid, w50)
      import tern_decoder
      implicit none
      class(tern_decoder), intent(inout) :: this
      integer, intent(in) :: nutc
      real,    intent(in) :: sync
      integer, intent(in) :: nsnr
      real,    intent(in) :: dt
      real,    intent(in) :: freq
      character(len=37), intent(in) :: decoded
      integer, intent(in) :: nap
      real,    intent(in) :: qual
      integer, intent(in) :: ntrperiod
      real,    intent(in) :: fmid
      real,    intent(in) :: w50
    end subroutine tern_decode_callback
  end interface

  !> Nominal Tx lead-in: the modulator starts audio delay_ms=1000 ms into the
  !! slot, so a perfectly-timed frame is estimated at t0 = 1.0 s.
  real(dp), parameter :: nominal_lead_s = 1.0_dp

contains

  subroutine decode(this, callback, iwave, nutc, nfa, nfb, nfqso, ntol, &
      nsubmode, ndepth, ntrperiod, nlen)
    use packjt77
    class(tern_decoder), intent(inout) :: this
    procedure(tern_decode_callback) :: callback
    integer*2, intent(in) :: iwave(*)
    integer,   intent(in) :: nutc, nfa, nfb, nfqso, ntol
    integer,   intent(in) :: nsubmode, ndepth, ntrperiod, nlen

    real(dp), parameter :: chunk_step_hz = 100.0_dp

    type(tern_rx_modem) :: modem
    type(demod_opts) :: opts
    type(scan_hit), allocatable :: hits(:)
    real(dp), allocatable :: pcm(:)
    character(len=1) :: mode_id
    real(dp) :: lo_hz, hi_hz, snr2500, snr_offset
    integer :: i, n, nsnr
    logical :: unpk_ok
    character*77 c77
    character*37 msgsent

    this%callback => callback

    if (nsubmode == 1) then
      mode_id = 'B'
    else
      mode_id = 'A'
    end if
    call modem%init(mode_id)

    n = max(0, nlen)
    if (n < 2 * nint(audio_rate_hz)) return       ! need a couple of seconds
    allocate(pcm(n))
    do i = 1, n
      pcm(i) = real(iwave(i), dp) / 32767.0_dp
    enddo

    opts%sync_search  = .true.
    opts%rx_spread_hz = 0.5_dp
    opts%max_iters    = max(1, iand(ndepth, 3) + 1)

    lo_hz = real(max(200, nfa), dp)
    hi_hz = real(min(4000, nfb), dp)
    if (hi_hz - lo_hz < chunk_step_hz) then        ! fall back to nfqso +- tol
      lo_hz = real(nfqso - max(ntol, 50), dp)
      hi_hz = real(nfqso + max(ntol, 50), dp)
    endif

    call scan_audio(modem, pcm, opts, lo_hz, hi_hz, hits)

    ! Es/N0 (per symbol) -> SNR in 2500 Hz, the TERN/WSJT-X reference.
    snr_offset = 10.0_dp * log10((1.0_dp / modem%mode%t_sym) / snr_ref_bw_hz)

    do i = 1, size(hits)
      write(c77, '(77i1)') hits(i)%msg_bits
      call unpack77(c77, 0, msgsent, unpk_ok)
      if (.not. unpk_ok) cycle
      snr2500 = hits(i)%es_n0_db + snr_offset
      nsnr = nint(snr2500)
      if (nsnr < -35) nsnr = -35
      if (nsnr >  49) nsnr =  49
      call this%callback(nutc, real(hits(i)%sync_score), nsnr,           &
           real(hits(i)%t0_s - nominal_lead_s), real(hits(i)%audio_hz),  &
           msgsent, 0, 1.0, ntrperiod, -999.0, -999.0)
    enddo

    call combine_at_rx(this, pcm, modem, opts, nfqso, ntol, nutc,        &
                       ntrperiod, snr_offset, hits)

    deallocate(pcm)
  end subroutine decode

  !> Cross-frame combining at the Rx frequency. Accumulates pilot-only tone
  !! LLRs of weak frames across slots and reports a decode_combined success.
  subroutine combine_at_rx(this, pcm, modem, opts, nfqso, ntol, nutc,    &
                           ntrperiod, snr_offset, hits)
    use packjt77
    class(tern_decoder), intent(inout) :: this
    real(dp), intent(in) :: pcm(:)
    type(tern_rx_modem), intent(in) :: modem
    type(demod_opts), intent(in) :: opts
    integer, intent(in) :: nfqso, ntol, nutc, ntrperiod
    real(dp), intent(in) :: snr_offset
    type(scan_hit), allocatable, intent(in) :: hits(:)

    type(demod_opts) :: qopts
    type(tern_result) :: rq
    complex(dp), allocatable :: xq(:)
    integer :: cmsg(k_msg), nsnr, i
    logical :: ccrc, unpk_ok, already
    real(dp) :: fq, snr2500
    character*77 c77
    character*37 msgsent

    if (nutc == acc_utc) return                ! one accumulation per slot
    fq = real(nfqso, dp)

    qopts = opts
    qopts%sync_search = .true.
    qopts%f_search_hz = real(max(ntol, 20), dp)
    call audio_to_baseband(pcm, fq, xq)
    call demodulate(modem%mode, modem%info_mask, xq, qopts, rq)
    if (allocated(xq)) deallocate(xq)

    if (rq%crc_ok) then                        ! decoded alone -> fresh start
      nacc = 0
      if (allocated(logp_acc)) deallocate(logp_acc)
      return
    end if
    if (.not. allocated(rq%tone_logp)) return

    if (nacc == 0 .or. abs(fq - acc_freq) > 50.0_dp) then
      logp_acc = rq%tone_logp
      nacc = 1
    else
      logp_acc = logp_acc + rq%tone_logp
      nacc = nacc + 1
    end if
    acc_utc = nutc
    acc_freq = fq

    if (nacc >= 2) then
      call decode_combined(modem%info_mask, logp_acc, cmsg, ccrc)
      if (ccrc) then
        ! Suppress if the scan already reported this exact message this slot.
        already = .false.
        do i = 1, size(hits)
          if (all(hits(i)%msg_bits == cmsg)) already = .true.
        enddo
        if (.not. already) then
          write(c77, '(77i1)') cmsg
          call unpack77(c77, 0, msgsent, unpk_ok)
          if (unpk_ok) then
            snr2500 = rq%es_n0_db + snr_offset + 10.0_dp * log10(real(nacc, dp))
            nsnr = nint(snr2500)
            if (nsnr < -35) nsnr = -35
            if (nsnr >  49) nsnr =  49
            ! nap = combine depth, shown as the "a<n>" annotation.
            call this%callback(nutc, real(rq%sync_score), nsnr,          &
                 real(rq%t0_s - nominal_lead_s), real(fq), msgsent,      &
                 nacc, 1.0, ntrperiod, -999.0, -999.0)
          endif
        endif
        nacc = 0
        if (allocated(logp_acc)) deallocate(logp_acc)
      endif
    endif
    if (nacc > max_combine) then
      nacc = 0
      if (allocated(logp_acc)) deallocate(logp_acc)
    endif
  end subroutine combine_at_rx

end module tern_decode
