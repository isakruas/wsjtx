!> @file ternsim.f90
!! @brief TERN validation driver for the WSJT-X build: TX -> ITU-R F.1487
!!        Watterson channel -> AWGN -> RX, sweeping SNR, mirroring the TERN
!!        reference campaign (tern_api::simulate_frame). Reproduces the decode
!!        thresholds of docs/VALIDATION.md and checks that the reported SNR
!!        (Es/N0 estimate converted to the 2500 Hz reference, as tern_decode
!!        does) tracks the input SNR.
!!
!! Usage: ternsim [mode A|B] [profile] [snr_lo snr_hi snr_step] [ntrials]
!! Example: ternsim A none -28 -16 1 200
program ternsim
  use tern_kinds,   only: dp
  use tern_params,  only: k_msg, make_mode, tern_mode, fs_hz, snr_ref_bw_hz, &
                          itu_profile, get_itu_profile
  use tern_frame,   only: encode_frame_tones
  use tern_mod,     only: synthesize_tones
  use tern_demod,   only: demod_opts, tern_result, demodulate
  use tern_channel, only: pad_signal, watterson_apply, apply_awgn, signal_power
  use tern_rng,     only: rng_state, rng_init, rng_uniform
  use tern_glue,    only: tern_info_mask
  implicit none

  character(len=8)  :: mode_arg
  character(len=24) :: prof_arg
  character(len=32) :: a
  type(tern_mode)   :: mode
  type(itu_profile) :: profile
  type(demod_opts)  :: opts
  type(tern_result) :: res
  type(rng_state)   :: rng
  integer, allocatable :: mask(:)
  integer :: msg(k_msg), tones(85)
  complex(dp), allocatable :: tx(:), padded(:), faded(:), noisy(:)
  real(dp) :: snr_lo, snr_hi, snr_step, snr, p_active, nvar, snr_offset
  real(dp) :: sum_rep, snr2500
  integer :: ntrials, ndec, t, i, nargs
  character(len=1) :: mid

  ! defaults
  mode_arg = 'A'; prof_arg = 'none'
  snr_lo = -28.0_dp; snr_hi = -16.0_dp; snr_step = 1.0_dp; ntrials = 200
  nargs = command_argument_count()
  if (nargs >= 1) call get_command_argument(1, mode_arg)
  if (nargs >= 2) call get_command_argument(2, prof_arg)
  if (nargs >= 5) then
    call get_command_argument(3, a); read(a,*) snr_lo
    call get_command_argument(4, a); read(a,*) snr_hi
    call get_command_argument(5, a); read(a,*) snr_step
  end if
  if (nargs >= 6) then
    call get_command_argument(6, a); read(a,*) ntrials
  end if

  mid = mode_arg(1:1)
  mode = make_mode(mid)
  profile = get_itu_profile(trim(prof_arg))
  call tern_info_mask(mask)
  snr_offset = 10.0_dp * log10((1.0_dp / mode%t_sym) / snr_ref_bw_hz)

  opts%sync_search = .true.
  opts%f_search_hz = 20.0_dp
  opts%rx_spread_hz = 0.5_dp
  opts%max_iters = 3

  print '(a,a,a,a,a,f6.2,a)', '# TERN-', mid, '  profile=', trim(prof_arg), &
        '  SNR offset(Es/N0->2500Hz)=', snr_offset, ' dB'
  print '(a)', '#  SNR_in   decode%   mean_reported_SNR   bias'
  snr = snr_lo
  do while (snr <= snr_hi + 1.0e-9_dp)
    ndec = 0; sum_rep = 0.0_dp
    do t = 1, ntrials
      call rng_init(rng, 100003 * t + nint(1000.0_dp * snr))
      do i = 1, k_msg
        msg(i) = merge(1, 0, rng_uniform(rng) > 0.5_dp)
      end do
      call encode_frame_tones(mask, msg, tones)
      call synthesize_tones(mode, tones, tx)
      p_active = signal_power(tx)
      call pad_signal(tx, 1.0_dp, 0.5_dp, fs_hz, padded)
      call watterson_apply(profile, padded, fs_hz, rng, faded)
      call apply_awgn(faded, snr, fs_hz, rng, noisy, nvar, ref_power=p_active)
      call demodulate(mode, mask, noisy, opts, res)
      if (res%crc_ok) then
        if (all(res%msg_bits == msg)) then
          ndec = ndec + 1
          snr2500 = res%es_n0_db + snr_offset
          sum_rep = sum_rep + snr2500
        end if
      end if
    end do
    if (ndec > 0) then
      print '(f8.1,f10.1,f15.2,f12.2)', snr, 100.0_dp*ndec/ntrials, &
            sum_rep/ndec, sum_rep/ndec - snr
    else
      print '(f8.1,f10.1,a)', snr, 100.0_dp*ndec/ntrials, '             --          --'
    end if
    snr = snr + snr_step
  end do
end program ternsim
