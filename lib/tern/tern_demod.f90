!> @file tern_demod.f90
!! @brief Block-noncoherent receiver with Gauss-Markov channel smoothing.
!!
!! @details
!! Detection model. After tone-bank projection, symbol s yields complex bin
!! outputs y_{s,j}. Under hypothesis "tone j sent", y_{s,j} = sqrt(Es) h_s +
!! n with h_s the flat-fading channel gain and n circular Gaussian of
!! variance sigma^2; all other bins are noise only. The Watterson process
!! has Gaussian autocorrelation R(tau) = exp(-2 pi^2 sigma_f^2 tau^2); at the
!! symbol lag it is approximated by an AR(1) (Gauss-Markov) process
!!
!!     h_{s+1} = a h_s + w_s,   a = R(T_s),   var(w) = (1 - a^2) P_h,
!!
!! and the posterior p(h_s | pilots [, decisions]) is computed exactly for
!! that model by a Kalman forward pass and RTS backward smoothing. The
!! resulting per-symbol Gaussian posterior CN(hhat_s, P_s) gives marginal
!! tone likelihoods
!!
!!     log p(y | tone j) = |y_j|^2/sigma^2
!!                       - |y_j - sqrt(Es) hhat_s|^2 / (sigma^2 + Es P_s)
!!                       - log(1 + Es P_s / sigma^2) + const,
!!
!! which interpolate continuously between coherent detection (P_s -> 0 near
!! pilots) and classical noncoherent square-law detection (P_s -> P_h far
!! from any anchor). Bit LLRs follow by log-sum-exp over the tone alphabet.
!!
!! After the first polar decode attempt, decided data tones are fed back as
!! channel observations (decision-directed turbo iteration), tightening the
!! posterior and re-deriving LLRs. Iteration stops on CRC success.
!!
!! Impulsive-noise robustness. Atmospheric HF noise is heavy-tailed
!! (ITU-R P.372); a Gaussian-derived metric is fragile against it because a
!! single impulse can dominate a symbol correlation. The receiver therefore
!! applies a memoryless envelope soft limiter before the tone bank: samples
!! whose magnitude exceeds three times the median envelope are scaled down
!! to that bound. Under circular Gaussian noise the median envelope is
!! sigma*sqrt(ln 4), the threshold sits at 3.53 sigma, and fewer than
!! 5e-6 of the samples are touched, so the limiter is statistically inert
!! on the Gaussian channel while bounding the influence of any impulse —
!! the classical locally-optimal detector structure for Middleton Class A
!! noise reduced to its robust-statistics core.
module tern_demod
  use tern_kinds, only: dp
  use tern_params, only: tern_mode, fs_hz, two_pi, pi, eps, m_tones, &
                         n_symbols, n_data_symbols, n_pilot_symbols, &
                         bits_per_symbol, k_msg, k_info, n_code, &
                         pilot_positions, data_positions, pilot_tone
  use tern_frame, only: tones_to_llr_order, coded_bit_index
  use tern_polar, only: polar_decode, polar_encode
  implicit none
  private
  public :: demod_opts, tern_result, acquire, demodulate, decode_combined

  !> Receiver configuration.
  type :: demod_opts
    logical  :: sync_search   = .true.   !< Run time-frequency acquisition.
    real(dp) :: f_search_hz   = 20.0_dp  !< One-sided CFO search range.
    real(dp) :: rx_spread_hz  = 0.5_dp   !< Assumed Doppler spread (2 sigma).
    integer  :: max_iters     = 3        !< Total decode attempts (1 + DD).
  end type demod_opts

  !> Decoding outcome and diagnostics.
  type :: tern_result
    integer  :: msg_bits(k_msg) = 0
    logical  :: crc_ok = .false.
    integer  :: iterations = 0
    real(dp) :: t0_s = 0.0_dp        !< Estimated frame start time.
    real(dp) :: f0_hz = 0.0_dp       !< Estimated carrier offset.
    real(dp) :: sync_score = 0.0_dp  !< Normalized pilot-energy metric.
    real(dp) :: es_n0_db = 0.0_dp    !< Estimated per-symbol Es over noise PSD bin.
    real(dp) :: noise_var = 0.0_dp   !< Estimated per-bin noise variance.
    !> Per-symbol tone log-likelihoods from the pilot-only channel posterior,
    !! retained for cross-frame combining. The pilot-only pass is used, not
    !! the decision-directed iterations, so that a frame that failed its CRC
    !! cannot contaminate the combining statistic with wrong decisions.
    real(dp), allocatable :: tone_logp(:, :)
  end type tern_result

contains

  !> Tone-bank projection of one aligned frame.
  !!
  !! Bin amplitudes are normalized so the per-bin noise variance equals the
  !! per-sample variance. Tone phase reference is the symbol start, which the
  !! integer-cycle grid makes consistent across symbols; the CFO term uses
  !! global time so that the implied channel gain remains continuous in s.
  subroutine tone_bank(mode, rx, start0, f0, y)
    type(tern_mode), intent(in) :: mode
    complex(dp), intent(in) :: rx(:)
    integer, intent(in) :: start0          !< Zero-based frame start sample.
    real(dp), intent(in) :: f0
    complex(dp), intent(out) :: y(m_tones, n_symbols)

    complex(dp) :: acc, rot
    real(dp) :: w, ph
    integer :: s, j, k, base

    do s = 1, n_symbols
      base = start0 + (s - 1) * mode%ns
      do j = 1, m_tones
        w = two_pi * (mode%tone_freq(j - 1) + f0) / fs_hz
        acc = (0.0_dp, 0.0_dp)
        do k = 0, mode%ns - 1
          ph = w * real(k, dp)
          acc = acc + rx(base + k + 1) * cmplx(cos(ph), -sin(ph), dp)
        end do
        ph = two_pi * f0 * real(base, dp) / fs_hz
        rot = cmplx(cos(ph), -sin(ph), dp)
        y(j, s) = acc * rot / sqrt(real(mode%ns, dp))
      end do
    end do
  end subroutine tone_bank

  !> Time-frequency acquisition by noncoherent pilot-energy maximization.
  !!
  !! Coarse grid: T_s/8 in time, df/8 in frequency, refined to one sample and
  !! df/64. The 21-symbol Costas pattern yields a single dominant peak on the
  !! (t, f) surface; the metric is invariant to the fading phase.
  subroutine acquire(mode, rx, f_range_hz, t0, f0, score)
    type(tern_mode), intent(in) :: mode
    complex(dp), intent(in) :: rx(:)
    real(dp), intent(in) :: f_range_hz
    integer, intent(out) :: t0
    real(dp), intent(out) :: f0, score

    integer :: ppos(n_pilot_symbols)
    integer :: n_search, t_step, t, p, k, fi, nf, best_t
    real(dp) :: f, f_step, best_f, best_m, m, mean_e
    complex(dp), allocatable :: tab(:, :)
    complex(dp) :: acc
    integer :: base

    ppos = pilot_positions()
    n_search = size(rx) - mode%frame_samples()
    if (n_search < 0) then
      t0 = 0; f0 = 0.0_dp; score = 0.0_dp
      return
    end if

    allocate(tab(mode%ns, n_pilot_symbols))
    t_step = max(1, mode%ns / 8)
    f_step = mode%df / 8.0_dp
    nf = nint(f_range_hz / f_step)

    best_m = -1.0_dp
    best_t = 0
    best_f = 0.0_dp
    do fi = -nf, nf
      f = real(fi, dp) * f_step
      call build_pilot_tables(mode, f, tab)
      do t = 0, n_search, t_step
        m = pilot_metric(mode, rx, t, tab, ppos)
        if (m > best_m) then
          best_m = m; best_t = t; best_f = f
        end if
      end do
    end do

    ! Local refinement to sample/fine-frequency resolution.
    block
      integer :: tt
      real(dp) :: ff
      do fi = -8, 8
        ff = best_f + real(fi, dp) * f_step / 8.0_dp
        call build_pilot_tables(mode, ff, tab)
        do tt = max(0, best_t - t_step), min(n_search, best_t + t_step)
          m = pilot_metric(mode, rx, tt, tab, ppos)
          if (m > best_m) then
            best_m = m; best_t = tt; best_f = ff
          end if
        end do
      end do
    end block

    t0 = best_t
    f0 = best_f
    mean_e = sum(abs(rx)**2) / max(1, size(rx))
    score = best_m / max(real(n_pilot_symbols * mode%ns, dp) * mean_e, eps)
    deallocate(tab)
  end subroutine acquire

  subroutine build_pilot_tables(mode, f0, tab)
    type(tern_mode), intent(in) :: mode
    real(dp), intent(in) :: f0
    complex(dp), intent(out) :: tab(:, :)
    integer :: p, k
    real(dp) :: w, ph
    do p = 1, n_pilot_symbols
      w = two_pi * (mode%tone_freq(pilot_tone(p)) + f0) / fs_hz
      do k = 0, mode%ns - 1
        ph = w * real(k, dp)
        tab(k + 1, p) = cmplx(cos(ph), -sin(ph), dp)
      end do
    end do
  end subroutine build_pilot_tables

  function pilot_metric(mode, rx, t, tab, ppos) result(m)
    type(tern_mode), intent(in) :: mode
    complex(dp), intent(in) :: rx(:), tab(:, :)
    integer, intent(in) :: t, ppos(:)
    real(dp) :: m
    complex(dp) :: acc
    integer :: p, k, base
    m = 0.0_dp
    do p = 1, n_pilot_symbols
      base = t + (ppos(p) - 1) * mode%ns
      acc = (0.0_dp, 0.0_dp)
      do k = 1, mode%ns
        acc = acc + rx(base + k) * tab(k, p)
      end do
      m = m + abs(acc)**2 / real(mode%ns, dp)
    end do
  end function pilot_metric

  !> Full receiver: acquisition, channel smoothing, LLR computation, and
  !! CRC-aided polar decoding with decision-directed iterations.
  subroutine demodulate(mode, info_mask, rx, opts, res)
    type(tern_mode), intent(in) :: mode
    integer, intent(in) :: info_mask(:)
    complex(dp), intent(in) :: rx(:)
    type(demod_opts), intent(in) :: opts
    type(tern_result), intent(out) :: res

    complex(dp) :: y(m_tones, n_symbols)
    complex(dp) :: z(n_symbols)
    logical :: has_obs(n_symbols)
    complex(dp) :: hhat(n_symbols)
    real(dp) :: pvar(n_symbols)
    real(dp) :: llr_sym(bits_per_symbol, n_data_symbols), llr_coded(n_code)
    real(dp) :: logp(m_tones, n_data_symbols)
    integer :: ppos(n_pilot_symbols), dpos(n_data_symbols)
    integer :: info_bits(k_info), coded(n_code), tones_dd(n_data_symbols)
    integer :: t0, iter, p, d, s, b, j
    real(dp) :: f0, score, sigma2, es, f_fine, a_corr, p_h, r_obs
    logical :: crc

    ppos = pilot_positions()
    dpos = data_positions()

    res = tern_result()
    if (size(rx) < mode%frame_samples()) return

    block
      complex(dp), allocatable :: rx_lim(:)
      call envelope_soft_limit(rx, rx_lim)

      if (opts%sync_search) then
        call acquire(mode, rx_lim, opts%f_search_hz, t0, f0, score)
      else
        t0 = 0; f0 = 0.0_dp; score = 1.0_dp
      end if

      call tone_bank(mode, rx_lim, t0, f0, y)
      call estimate_noise_es(y, ppos, sigma2, es)

      ! Fine CFO from within-block adjacent pilot pairs, then one
      ! re-projection of the limited samples.
      call fine_cfo(mode, y, ppos, sigma2, f_fine)
      if (abs(f_fine) > 1.0e-4_dp) then
        f0 = f0 + f_fine
        call tone_bank(mode, rx_lim, t0, f0, y)
        call estimate_noise_es(y, ppos, sigma2, es)
      end if
      deallocate(rx_lim)
    end block

    res%t0_s = real(t0, dp) / fs_hz
    res%f0_hz = f0
    res%sync_score = score
    res%noise_var = sigma2
    res%es_n0_db = 10.0_dp * log10(max(es / max(sigma2, eps), eps))

    ! Gauss-Markov parameters from the assumed Doppler spread.
    a_corr = exp(-2.0_dp * pi * pi * (0.5_dp * opts%rx_spread_hz)**2 * mode%t_sym**2)
    a_corr = min(max(a_corr, 0.0_dp), 0.999999_dp)
    r_obs = sigma2 / max(es, eps)

    ! Pilot channel observations normalized to unit nominal gain.
    has_obs = .false.
    z = (0.0_dp, 0.0_dp)
    do p = 1, n_pilot_symbols
      z(ppos(p)) = y(pilot_tone(p) + 1, ppos(p)) / sqrt(max(es, eps))
      has_obs(ppos(p)) = .true.
    end do
    p_h = max(0.1_dp, min(4.0_dp, channel_power(z, has_obs) - r_obs))

    allocate(res%tone_logp(m_tones, n_data_symbols))
    do iter = 1, max(1, opts%max_iters)
      call gauss_markov_smooth(z, has_obs, a_corr, p_h, r_obs, hhat, pvar)
      call tone_logp_matrix(y, dpos, hhat, pvar, sigma2, es, logp)
      if (iter == 1) res%tone_logp = logp
      call llrs_from_logp(logp, llr_sym)
      call tones_to_llr_order(llr_sym, llr_coded)
      call polar_decode(info_mask, llr_coded, info_bits, crc)
      res%iterations = iter
      if (crc) exit
      if (iter == opts%max_iters) exit
      ! Decision feedback: decoded tones become channel observations with the
      ! observation noise inflated to hedge residual decision errors.
      call polar_encode(info_mask, info_bits, coded)
      do d = 1, n_data_symbols
        j = 0
        do b = 1, bits_per_symbol
          j = ior(ishft(j, 1), iand(coded(coded_bit_index(d, b)), 1))
        end do
        tones_dd(d) = j
        z(dpos(d)) = y(j + 1, dpos(d)) / sqrt(max(es, eps))
        has_obs(dpos(d)) = .true.
      end do
      r_obs = 1.5_dp * sigma2 / max(es, eps)
    end do

    res%crc_ok = crc
    res%msg_bits = info_bits(1:k_msg)
  end subroutine demodulate

  !> Memoryless envelope soft limiter bounding impulsive-noise influence.
  !!
  !! The clipping bound is three times the median envelope, estimated from a
  !! uniform subsample of at most 4096 points. Samples below the bound pass
  !! unchanged; samples above it keep their phase and are scaled to the
  !! bound. See the module preamble for the statistical rationale.
  subroutine envelope_soft_limit(rx, out)
    complex(dp), intent(in) :: rx(:)
    complex(dp), allocatable, intent(out) :: out(:)
    real(dp), allocatable :: env(:)
    real(dp) :: med, bound, a
    integer :: n, ne, step, i

    n = size(rx)
    step = max(1, n / 4096)
    ne = (n + step - 1) / step
    allocate(env(ne))
    do i = 1, ne
      env(i) = abs(rx((i - 1) * step + 1))
    end do
    med = median_of(env)
    deallocate(env)

    allocate(out(n))
    if (med <= 0.0_dp) then
      out = rx
      return
    end if
    bound = 3.0_dp * med
    do i = 1, n
      a = abs(rx(i))
      if (a > bound) then
        out(i) = rx(i) * (bound / a)
      else
        out(i) = rx(i)
      end if
    end do
  end subroutine envelope_soft_limit

  !> Per-bin noise variance and pilot symbol energy from pilot symbols.
  !! Noise uses the median of off-tone bin energies (exponential median is
  !! sigma^2 ln 2), robust to GFSK spectral leakage and narrowband QRM.
  subroutine estimate_noise_es(y, ppos, sigma2, es)
    complex(dp), intent(in) :: y(:, :)
    integer, intent(in) :: ppos(:)
    real(dp), intent(out) :: sigma2, es
    real(dp) :: e(n_pilot_symbols * (m_tones - 1)), ep
    integer :: p, j, k, n
    n = 0
    ep = 0.0_dp
    do p = 1, n_pilot_symbols
      k = pilot_tone(p)
      do j = 1, m_tones
        if (j - 1 == k) cycle
        n = n + 1
        e(n) = abs(y(j, ppos(p)))**2
      end do
      ep = ep + abs(y(k + 1, ppos(p)))**2
    end do
    sigma2 = max(median_of(e(1:n)) / log(2.0_dp), eps)
    es = max(ep / real(n_pilot_symbols, dp) - sigma2, 10.0_dp * eps)
  end subroutine estimate_noise_es

  !> Residual CFO from the average phase increment between adjacent pilot
  !! symbols inside each Costas block. The fading contribution is zero-mean;
  !! 18 pair products average it down.
  subroutine fine_cfo(mode, y, ppos, sigma2, f_fine)
    type(tern_mode), intent(in) :: mode
    complex(dp), intent(in) :: y(:, :)
    integer, intent(in) :: ppos(:)
    real(dp), intent(in) :: sigma2
    real(dp), intent(out) :: f_fine
    complex(dp) :: acc, z1, z2
    integer :: blk, k, p
    acc = (0.0_dp, 0.0_dp)
    do blk = 0, 2
      do k = 1, 6
        p = blk * 7 + k
        z1 = y(pilot_tone(p) + 1, ppos(p))
        z2 = y(pilot_tone(p + 1) + 1, ppos(p + 1))
        acc = acc + z2 * conjg(z1)
      end do
    end do
    if (abs(acc) < 10.0_dp * sigma2) then
      f_fine = 0.0_dp
    else
      f_fine = atan2(aimag(acc), real(acc, dp)) / (two_pi * mode%t_sym)
    end if
  end subroutine fine_cfo

  !> Mean observed channel power over observation symbols.
  function channel_power(z, has_obs) result(p)
    complex(dp), intent(in) :: z(:)
    logical, intent(in) :: has_obs(:)
    real(dp) :: p
    integer :: s, n
    p = 0.0_dp
    n = 0
    do s = 1, size(z)
      if (.not. has_obs(s)) cycle
      n = n + 1
      p = p + abs(z(s))**2
    end do
    p = p / max(1, n)
  end function channel_power

  !> Exact posterior for the AR(1) channel model: Kalman filter forward,
  !! Rauch-Tung-Striebel smoother backward.
  subroutine gauss_markov_smooth(z, has_obs, a, p_h, r, hhat, pvar)
    complex(dp), intent(in) :: z(:)
    logical, intent(in) :: has_obs(:)
    real(dp), intent(in) :: a, p_h, r
    complex(dp), intent(out) :: hhat(:)
    real(dp), intent(out) :: pvar(:)

    complex(dp) :: mf(n_symbols), mp(n_symbols)
    real(dp) :: pf(n_symbols), pp(n_symbols)
    complex(dp) :: m_prev
    real(dp) :: p_prev, q, kgain, c
    integer :: s

    q = (1.0_dp - a * a) * p_h
    m_prev = (0.0_dp, 0.0_dp)
    p_prev = p_h
    do s = 1, n_symbols
      if (s == 1) then
        mp(s) = (0.0_dp, 0.0_dp)
        pp(s) = p_h
      else
        mp(s) = a * m_prev
        pp(s) = a * a * p_prev + q
      end if
      if (has_obs(s)) then
        kgain = pp(s) / (pp(s) + r)
        mf(s) = mp(s) + kgain * (z(s) - mp(s))
        pf(s) = (1.0_dp - kgain) * pp(s)
      else
        mf(s) = mp(s)
        pf(s) = pp(s)
      end if
      m_prev = mf(s)
      p_prev = pf(s)
    end do

    hhat(n_symbols) = mf(n_symbols)
    pvar(n_symbols) = pf(n_symbols)
    do s = n_symbols - 1, 1, -1
      c = a * pf(s) / max(pp(s + 1), eps)
      hhat(s) = mf(s) + c * (hhat(s + 1) - mp(s + 1))
      pvar(s) = pf(s) + c * c * (pvar(s + 1) - pp(s + 1))
      pvar(s) = max(pvar(s), 0.0_dp)
    end do
  end subroutine gauss_markov_smooth

  !> Per-symbol tone log-likelihoods from the Gaussian channel posterior.
  subroutine tone_logp_matrix(y, dpos, hhat, pvar, sigma2, es, logp)
    complex(dp), intent(in) :: y(:, :), hhat(:)
    integer, intent(in) :: dpos(:)
    real(dp), intent(in) :: pvar(:), sigma2, es
    real(dp), intent(out) :: logp(:, :)

    real(dp) :: v
    complex(dp) :: mu
    integer :: d, s, j

    do d = 1, n_data_symbols
      s = dpos(d)
      v = sigma2 + es * pvar(s)
      mu = sqrt(es) * hhat(s)
      do j = 1, m_tones
        logp(j, d) = abs(y(j, s))**2 / sigma2 - abs(y(j, s) - mu)**2 / v - &
                     log(v / sigma2)
      end do
    end do
  end subroutine tone_logp_matrix

  !> Marginal bit LLRs from per-symbol tone log-likelihoods.
  subroutine llrs_from_logp(logp, llr_sym)
    real(dp), intent(in) :: logp(:, :)
    real(dp), intent(out) :: llr_sym(:, :)
    real(dp) :: l0, l1
    integer :: d, j, b

    do d = 1, n_data_symbols
      do b = 1, bits_per_symbol
        l0 = -huge(1.0_dp)
        l1 = -huge(1.0_dp)
        do j = 1, m_tones
          if (ibits(j - 1, bits_per_symbol - b, 1) == 0) then
            l0 = logsumexp2(l0, logp(j, d))
          else
            l1 = logsumexp2(l1, logp(j, d))
          end if
        end do
        llr_sym(b, d) = max(-40.0_dp, min(40.0_dp, l0 - l1))
      end do
    end do
  end subroutine llrs_from_logp

  !> Decode from accumulated tone log-likelihoods.
  !!
  !! Cross-frame combining: when the same message is transmitted in several
  !! frames over independent channel realizations, the per-symbol tone
  !! log-likelihoods are additive, so summing the tone_logp matrices of the
  !! individual receptions and decoding the sum realizes full diversity
  !! combining — every failed reception still contributes its energy,
  !! in contrast with discard-and-retry protocols.
  !!
  !! @param[in]  info_mask Polar information-set mask.
  !! @param[in]  logp_sum  Accumulated (m_tones, n_data_symbols) matrix.
  !! @param[out] msg_bits  Decoded message bits.
  !! @param[out] crc_ok    True when the decode satisfies the frame CRC.
  subroutine decode_combined(info_mask, logp_sum, msg_bits, crc_ok)
    integer, intent(in) :: info_mask(:)
    real(dp), intent(in) :: logp_sum(:, :)
    integer, intent(out) :: msg_bits(k_msg)
    logical, intent(out) :: crc_ok

    real(dp) :: llr_sym(bits_per_symbol, n_data_symbols), llr_coded(n_code)
    integer :: info_bits(k_info)

    call llrs_from_logp(logp_sum, llr_sym)
    call tones_to_llr_order(llr_sym, llr_coded)
    call polar_decode(info_mask, llr_coded, info_bits, crc_ok)
    msg_bits = info_bits(1:k_msg)
  end subroutine decode_combined

  pure function logsumexp2(a, b) result(c)
    real(dp), intent(in) :: a, b
    real(dp) :: c, hi, lo
    hi = max(a, b)
    lo = min(a, b)
    if (hi - lo > 40.0_dp) then
      c = hi
    else
      c = hi + log(1.0_dp + exp(lo - hi))
    end if
  end function logsumexp2

  function median_of(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    real(dp), allocatable :: w(:)
    integer :: i, j, n
    real(dp) :: key
    n = size(x)
    allocate(w(n))
    w = x
    do i = 2, n
      key = w(i)
      j = i - 1
      do while (j >= 1)
        if (w(j) <= key) exit
        w(j + 1) = w(j)
        j = j - 1
      end do
      w(j + 1) = key
    end do
    if (mod(n, 2) == 1) then
      m = w((n + 1) / 2)
    else
      m = 0.5_dp * (w(n / 2) + w(n / 2 + 1))
    end if
    deallocate(w)
  end function median_of
end module tern_demod
