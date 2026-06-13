!> @file gentern.f90
!! @brief TERN transmit shims callable from WSJT-X (mainwindow.cpp), mirroring
!!        the genfst4 / gen_fst4wave pair.
!!
!!   gentern       message text -> pack77 -> 77 bits -> 85-symbol tone sequence
!!   gen_ternwave  tone sequence -> real audio passband waveform at fsample
!!
!! TERN performs its own 14-bit CRC and (256,91) polar coding inside
!! encode_frame_tones, so the raw 77 source-encoded bits from pack77 are fed
!! straight in (no FST4-style rvec scrambling): the TERN encode/decode pair is
!! self-consistent and round-trips the WSJT-X 77-bit payload.

subroutine gentern(msg0, ichk, msgsent, i4tone, nsubmode)
  use, intrinsic :: iso_c_binding
  use packjt77
  use tern_kinds,  only: dp
  use tern_params, only: k_msg, n_symbols
  use tern_frame,  only: encode_frame_tones
  use tern_glue,   only: tern_info_mask
  implicit none

  character*37 msg0, msgsent
  integer ichk
  integer i4tone(n_symbols)
  integer nsubmode

  character*37 message
  character*77 c77
  integer i3, n3, i
  integer msg(k_msg)
  integer, allocatable :: mask(:)
  logical unpk77_success

  message = msg0
  do i = 1, 37
    if (ichar(message(i:i)) .eq. 0) then
      message(i:37) = ' '
      exit
    endif
  enddo
  do i = 1, 37                              !Strip leading blanks
    if (message(1:1) .ne. ' ') exit
    message = message(i+1:)
  enddo

  i3 = -1
  n3 = -1
  call pack77(message, i3, n3, c77)
  call unpack77(c77, 0, msgsent, unpk77_success)  !Unpack to get msgsent

  i4tone = 0
  if (ichk .eq. 1) return
  if (.not. unpk77_success) then
    msgsent = '*** bad message ***                  '
    return
  endif

  read(c77, '(77i1)') msg(1:77)
  call tern_info_mask(mask)
  call encode_frame_tones(mask, msg, i4tone)
end subroutine gentern


subroutine gen_ternwave(i4tone, nsym, nsps, fsample, f0, cwave, nsubmode)
  use tern_kinds,  only: dp
  use tern_params, only: tern_mode
  use tern_mod,    only: synthesize_tones
  use tern_glue,   only: tern_mode_from_n, tern_audio_synth
  implicit none

  integer nsym, nsps, nsubmode
  integer i4tone(nsym)
  real fsample, f0
  real cwave(*)                 ! foxcom_.wave; written up to (nsym+2)*nsps

  type(tern_mode) :: mode
  complex(dp), allocatable :: iq(:)
  real(dp), allocatable :: pcm(:)
  integer :: i, ntot, ncopy

  mode = tern_mode_from_n(nsubmode)
  call synthesize_tones(mode, i4tone(1:nsym), iq)
  call tern_audio_synth(iq, real(f0, dp), real(fsample, dp), pcm)

  ntot  = nsym * nsps
  ncopy = min(ntot, size(pcm))
  do i = 1, ncopy
    cwave(i) = real(pcm(i))
  enddo
  ! Zero the remainder plus a two-symbol guard: the modulator's fade-out reads
  ! one sample past nsym*nsps (i1 = symbolsLength*4*nsps_12k).
  do i = ncopy + 1, ntot + 2 * nsps
    cwave(i) = 0.0
  enddo
end subroutine gen_ternwave
