!> @file tern_polar.f90
!! @brief Polar encoder and CRC-aided successive-cancellation list decoder.
!!
!! @details
!! Code construction uses Gaussian-approximation density evolution
!! (Trifonov's phi-function approximation) for a BI-AWGN design point, which
!! tracks the LLR mean through the polar transform and is markedly more
!! accurate than the Bhattacharyya/BEC heuristic for short codes.
!!
!! Decoding is LLR-domain SCL. Final path selection is CRC-aided: among the
!! surviving paths ordered by metric, the first whose information bits pass
!! the frame CRC is selected; if none passes, the minimum-metric path is
!! returned. CRC-aided selection converts the list into a near-ML decoder for
!! the concatenated CRC+polar code and is worth roughly 0.5-1 dB over plain
!! SCL at this block length.
!!
!! LLR sign convention: positive favors bit 0.
module tern_polar
  use tern_kinds, only: dp
  use tern_params, only: n_code, k_info, scl_list_size, eps, pi
  use tern_bits, only: check_crc
  implicit none
  private
  public :: polar_info_set, polar_encode, polar_decode

contains

  !> Gaussian-approximation construction. Returns a mask with 1 at the
  !! k_info most reliable synthetic-channel positions.
  !!
  !! @param[in]  design_llr_mean Mean of the channel LLR at the design point
  !!                             (4*Es/N0 per coded bit for BI-AWGN).
  !! @param[out] info_mask       Length n_code; 1 = information, 0 = frozen.
  subroutine polar_info_set(design_llr_mean, info_mask)
    real(dp), intent(in) :: design_llr_mean
    integer, allocatable, intent(out) :: info_mask(:)
    real(dp), allocatable :: m(:), mnext(:)
    integer, allocatable :: order(:)
    integer :: levels, lvl, half, i
    real(dp) :: p

    levels = nint(log(real(n_code, dp)) / log(2.0_dp))
    if (ishft(1, levels) /= n_code) error stop "polar_info_set: N not a power of 2"

    allocate(m(n_code), mnext(n_code))
    m = 0.0_dp
    m(1) = max(design_llr_mean, 1.0e-6_dp)
    do lvl = 1, levels
      half = ishft(1, lvl - 1)
      do i = 1, half
        p = phi_fun(m(i))
        mnext(2 * i - 1) = phi_inv(min(1.0_dp, 2.0_dp * p - p * p))
        mnext(2 * i)     = 2.0_dp * m(i)
      end do
      m(1:2 * half) = mnext(1:2 * half)
    end do

    allocate(order(n_code))
    do i = 1, n_code
      order(i) = i
    end do
    ! Descending reliability = descending LLR mean.
    call sort_indices(order, -m)

    allocate(info_mask(n_code))
    info_mask = 0
    do i = 1, k_info
      info_mask(order(i)) = 1
    end do
    deallocate(m, mnext, order)
  end subroutine polar_info_set

  !> phi(m) = 1 - E[tanh(L/2)] for L ~ N(m, 2m); Trifonov approximation.
  pure function phi_fun(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x <= 1.0e-9_dp) then
      y = 1.0_dp
    else if (x < 10.0_dp) then
      y = exp(-0.4527_dp * x**0.86_dp + 0.0218_dp)
      y = min(y, 1.0_dp)
    else
      y = sqrt(pi / x) * exp(-x / 4.0_dp) * (1.0_dp - 10.0_dp / (7.0_dp * x))
      y = max(y, 0.0_dp)
    end if
  end function phi_fun

  !> Inverse of phi_fun by bisection on its monotone-decreasing domain.
  pure function phi_inv(y) result(x)
    real(dp), intent(in) :: y
    real(dp) :: x
    real(dp) :: lo, hi, mid, val
    integer :: it
    if (y >= 1.0_dp) then
      x = 0.0_dp
      return
    end if
    if (y <= 1.0e-30_dp) then
      x = 300.0_dp
      return
    end if
    lo = 0.0_dp
    hi = 300.0_dp
    do it = 1, 80
      mid = 0.5_dp * (lo + hi)
      val = phi_fun(mid)
      if (val > y) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function phi_inv

  !> Encode k_info information bits into one n_code codeword.
  subroutine polar_encode(info_mask, info_bits, coded)
    integer, intent(in)  :: info_mask(:)
    integer, intent(in)  :: info_bits(:)
    integer, intent(out) :: coded(:)
    integer :: u(n_code)
    integer :: i, j, cnt, half, step

    u = 0
    cnt = 0
    do i = 1, n_code
      if (info_mask(i) == 1) then
        cnt = cnt + 1
        u(i) = iand(info_bits(cnt), 1)
      end if
    end do

    half = 1
    do while (half < n_code)
      step = half * 2
      do i = 1, n_code, step
        do j = 0, half - 1
          u(i + j) = ieor(u(i + j), u(i + j + half))
        end do
      end do
      half = step
    end do
    coded(1:n_code) = u
  end subroutine polar_encode

  !> CRC-aided SCL decoding of one codeword.
  !!
  !! @param[in]  info_mask Information-set mask from polar_info_set.
  !! @param[in]  llrs      n_code channel LLRs, positive favoring bit 0.
  !! @param[out] info_out  k_info decoded information bits.
  !! @param[out] crc_pass  True when the selected path satisfies the CRC.
  subroutine polar_decode(info_mask, llrs, info_out, crc_pass)
    integer, intent(in)  :: info_mask(:)
    real(dp), intent(in) :: llrs(:)
    integer, intent(out) :: info_out(:)
    logical, intent(out) :: crc_pass

    integer :: levels, l_max
    real(dp), allocatable :: llr_tree(:, :, :)
    integer,  allocatable :: u_tree(:, :, :)
    real(dp), allocatable :: path_metric(:)
    logical,  allocatable :: path_active(:)
    integer,  allocatable :: order(:), cand_info(:)
    real(dp), allocatable :: act_metric(:)
    integer,  allocatable :: act_idx(:)
    integer :: phi, l, i, n_act, best

    levels = nint(log(real(n_code, dp)) / log(2.0_dp))
    l_max = scl_list_size

    allocate(llr_tree(0:levels, n_code, l_max))
    allocate(u_tree(0:levels, n_code, l_max))
    allocate(path_metric(l_max), path_active(l_max))
    path_metric = 0.0_dp
    path_active = .false.
    path_active(1) = .true.
    do l = 1, l_max
      llr_tree(0, 1:n_code, l) = llrs(1:n_code)
      u_tree(:, :, l) = 0
    end do

    do phi = 0, n_code - 1
      do l = 1, l_max
        if (.not. path_active(l)) cycle
        call descend_llr(phi, levels, llr_tree(:, :, l), u_tree(:, :, l))
      end do

      if (info_mask(phi + 1) == 0) then
        do l = 1, l_max
          if (.not. path_active(l)) cycle
          path_metric(l) = path_metric(l) + log1p_exp(-llr_tree(levels, phi + 1, l))
          u_tree(levels, phi + 1, l) = 0
        end do
      else
        call split_paths(l_max, levels, phi, llr_tree, u_tree, path_metric, path_active)
      end if

      do l = 1, l_max
        if (.not. path_active(l)) cycle
        call ascend_partial_sum(phi, levels, u_tree(:, :, l))
      end do
    end do

    ! CRC-aided selection over metric-ordered surviving paths.
    n_act = count(path_active)
    allocate(act_metric(n_act), act_idx(n_act), order(n_act), cand_info(k_info))
    i = 0
    do l = 1, l_max
      if (.not. path_active(l)) cycle
      i = i + 1
      act_metric(i) = path_metric(l)
      act_idx(i) = l
    end do
    do i = 1, n_act
      order(i) = i
    end do
    call sort_indices(order, act_metric)

    crc_pass = .false.
    best = act_idx(order(1))
    do i = 1, n_act
      l = act_idx(order(i))
      call extract_info(info_mask, u_tree(:, :, l), levels, cand_info)
      if (check_crc(cand_info)) then
        best = l
        crc_pass = .true.
        exit
      end if
    end do
    call extract_info(info_mask, u_tree(:, :, best), levels, info_out)

    deallocate(llr_tree, u_tree, path_metric, path_active)
    deallocate(act_metric, act_idx, order, cand_info)
  end subroutine polar_decode

  !> Collect the information bits of one decoded path.
  subroutine extract_info(info_mask, u, levels, info_out)
    integer, intent(in) :: info_mask(:)
    integer, intent(in) :: u(0:, :)
    integer, intent(in) :: levels
    integer, intent(out) :: info_out(:)
    integer :: i, cnt
    cnt = 0
    do i = 1, n_code
      if (info_mask(i) == 1) then
        cnt = cnt + 1
        info_out(cnt) = u(levels, i)
      end if
    end do
  end subroutine extract_info

  !> Compute LLR messages from the channel level down to the current leaf.
  subroutine descend_llr(phi, levels, llr, u)
    integer, intent(in) :: phi, levels
    real(dp), intent(inout) :: llr(0:, :)
    integer,  intent(inout) :: u(0:, :)
    integer :: depth, span, half, idx, j, bit
    real(dp) :: la, lb

    do depth = 0, levels - 1
      span = ishft(1, levels - depth)
      half = span / 2
      idx  = (phi / span) * span + 1
      if (mod(phi / half, 2) == 0) then
        do j = 0, half - 1
          la = llr(depth, idx + j)
          lb = llr(depth, idx + j + half)
          llr(depth + 1, idx + j) = sign(1.0_dp, la) * sign(1.0_dp, lb) * &
                                    min(abs(la), abs(lb))
        end do
      else
        do j = 0, half - 1
          la = llr(depth, idx + j)
          lb = llr(depth, idx + j + half)
          bit = u(depth + 1, idx + j)
          if (bit == 0) then
            llr(depth + 1, idx + j + half) = lb + la
          else
            llr(depth + 1, idx + j + half) = lb - la
          end if
        end do
      end if
    end do
  end subroutine descend_llr

  !> Propagate a decided leaf upward through the partial-sum tree.
  subroutine ascend_partial_sum(phi, levels, u)
    integer, intent(in) :: phi, levels
    integer, intent(inout) :: u(0:, :)
    integer :: depth, span, half, idx, j

    do depth = levels - 1, 0, -1
      span = ishft(1, levels - depth)
      half = span / 2
      idx  = (phi / span) * span + 1
      if (mod(phi / half, 2) == 1) then
        do j = 0, half - 1
          u(depth, idx + j)        = ieor(u(depth + 1, idx + j), &
                                          u(depth + 1, idx + j + half))
          u(depth, idx + j + half) = u(depth + 1, idx + j + half)
        end do
      else
        return
      end if
    end do
  end subroutine ascend_partial_sum

  !> Duplicate every active path on an information bit and retain the l_max
  !! lowest-metric candidates.
  subroutine split_paths(l_max, levels, phi, llr_tree, u_tree, path_metric, path_active)
    integer, intent(in) :: l_max, levels, phi
    real(dp), intent(inout) :: llr_tree(0:, :, :)
    integer,  intent(inout) :: u_tree(0:, :, :)
    real(dp), intent(inout) :: path_metric(:)
    logical,  intent(inout) :: path_active(:)

    real(dp), allocatable :: cand_metric(:), new_llr(:, :, :), new_metric(:)
    integer,  allocatable :: cand_path(:), cand_bit(:), order(:), new_u(:, :, :)
    logical,  allocatable :: new_active(:)
    integer :: n_cand, k, l, n_keep, i, parent, bit, n
    real(dp) :: lv

    n_cand = 2 * count(path_active)
    allocate(cand_metric(n_cand), cand_path(n_cand), cand_bit(n_cand))
    k = 0
    do l = 1, l_max
      if (.not. path_active(l)) cycle
      lv = llr_tree(levels, phi + 1, l)
      k = k + 1; cand_metric(k) = path_metric(l) + log1p_exp(-lv)
      cand_path(k) = l; cand_bit(k) = 0
      k = k + 1; cand_metric(k) = path_metric(l) + log1p_exp(lv)
      cand_path(k) = l; cand_bit(k) = 1
    end do

    allocate(order(n_cand))
    do i = 1, n_cand
      order(i) = i
    end do
    call sort_indices(order, cand_metric)
    n_keep = min(l_max, n_cand)

    n = size(llr_tree, 2)
    allocate(new_llr(0:levels, n, l_max), new_u(0:levels, n, l_max))
    allocate(new_metric(l_max), new_active(l_max))
    new_metric = 0.0_dp
    new_active = .false.
    new_llr = 0.0_dp
    new_u = 0

    do i = 1, n_keep
      parent = cand_path(order(i))
      bit    = cand_bit(order(i))
      new_llr(:, :, i) = llr_tree(:, :, parent)
      new_u(:, :, i)   = u_tree(:, :, parent)
      new_u(levels, phi + 1, i) = bit
      new_metric(i) = cand_metric(order(i))
      new_active(i) = .true.
    end do

    llr_tree = new_llr
    u_tree = new_u
    path_metric = new_metric
    path_active = new_active

    deallocate(cand_metric, cand_path, cand_bit, order)
    deallocate(new_llr, new_u, new_metric, new_active)
  end subroutine split_paths

  !> Numerically stable log(1 + exp(x)).
  pure real(dp) function log1p_exp(x) result(y)
    real(dp), intent(in) :: x
    if (x > 30.0_dp) then
      y = x
    else if (x < -30.0_dp) then
      y = 0.0_dp
    else
      y = log(1.0_dp + exp(x))
    end if
  end function log1p_exp

  !> Insertion sort of an index permutation by ascending associated value.
  pure subroutine sort_indices(order, values)
    integer,  intent(inout) :: order(:)
    real(dp), intent(in) :: values(:)
    integer :: i, j, key
    do i = 2, size(order)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (values(order(j)) <= values(key)) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do
  end subroutine sort_indices
end module tern_polar
