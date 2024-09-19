! SPDX-License-Identifier: GPL-3.0-or-later


subroutine overlap_save(x, Nx, h, M, N, y)
  ! Note that H is assumed hermitian: if H_k = DFT(h_n) over N points,
  ! then
  !    for 0 <= k <= N/2:
  !    H_0 = H(0)
  !    H_{N/2} = H(1)
  !    H_k = H(2*k) + j H(2*k+1)
  !    
  ! and for -N/2 + 1 <= k <= -1
  !    H_k = H(-2*k) + j H(-2*k + 1)
  !
  ! WARNING: N must be EVEN
  use iso_c_binding, only: c_double, c_long
  implicit none

  ! signal to be convolved
  integer(kind=c_long), intent(in) :: Nx
  real(kind=c_double),  intent(in) :: x(0:Nx-1)

  ! filter parameters
  integer(kind=c_long), intent(in) :: N    ! length of DFT
  real(kind=c_double),  intent(in) :: h(0:M-1)    ! filter DFT
  integer(kind=c_long), intent(in) :: M    ! filter duration in time

  ! output sequence
  real(kind=c_double),  intent(inout):: y(0:Nx+M-2)


  ! routine auxiliaries
  real(kind=c_double) :: X_FFT(0:N-1), Y_FFT(0:N-1), H_FFT(0:N-1)
  integer(kind=c_long) :: x_ptr, x_cntr, x_left, y_left
  ! integer(kind=c_double) :: i
  integer(kind=c_long) :: N_new, L_tot
  real(kind=c_double) :: x_tmp(0:N-1)

  N_new = N - M + 1
  L_tot = Nx+M-1

  
  H_FFT(0:M-1) = h(:)
  H_FFT(M:)    = 0
  
  call rfft(H_FFT, N, -1)

  ! --- First iteration ---
  x_tmp(0:M-2) = 0
  x_tmp(M-1:N-1) = x(0:N_new-1)

  call overlap_save_iteration(x_tmp, H_FFT, N, M, y, N_new)

  ! update parameters for the next iteration
  x_ptr  = N_new
  x_cntr = N_new + N


  ! --- Other iterations ---
  ! Note that x_cntr leads y_cntr, so it "hits the top" before y_cntr does
  do while (x_cntr <= Nx)
     call overlap_save_iteration(x(x_ptr), H_FFT, N, M, y(x_ptr), N_new)
     x_cntr = x_cntr + N_new
     x_ptr  = x_ptr  + N_new
  end do


  ! --- latest iterations ---
  x_left = Nx - x_ptr
  y_left = L_tot - x_ptr
  do while (y_left > 0)
     ! note that x_ptr is, by chance, also the number of currently valid y samples
     if (x_left > 0) then 
        x_tmp(0:x_left-1) = x(x_ptr:)
        x_tmp(x_left:) = 0
     else
        ! should never happen
        x_tmp(:) = 0
     end if
     call overlap_save_iteration(x_tmp, H_FFT, N, M, y(x_ptr), y_left)
     
     x_ptr  = x_ptr  + N_new
     x_left = x_left - N_new
     y_left = y_left - N_new
  end do

end subroutine overlap_save


subroutine overlap_save_iteration(x, H, N, M, y, N_new)
  use iso_c_binding, only: c_double, c_long
  implicit none

  ! inout variables
  integer(kind=c_long), intent(in) :: N, M
  real(kind=c_double),  intent(in) :: x(0:N-1)

  real(kind=c_double),  intent(in) :: H(0:N-1)

  integer(kind=c_long), intent(in) :: N_new
  real(kind=c_double),  intent(inout) :: y(0:N_new-1)


  ! internal variables
  real(kind=c_double) :: x_fft(0:N-1), y_fft(0:N-1)
  
  x_fft(0:N-1) = x(0:N-1)

  call rfft(x_fft, N, -1)

  call mult_hermitian(x_fft, H, N, y_fft)

  call ifft_hermitian(y_fft, N)

  ! Note that normally N_new = N - M + 1,
  ! yielding N - 1 as the upper end,
  ! but this way we can also deal with the exit transition
  ! so less then N - M + 1 new data are needed
  if (N_new >= N - M + 1) then
     y(0:N-M) = y_fft(M-1 : N - 1)
  else
     y(:) = y_fft(M-1 : M - 2 + N_new)
  end if
  
end subroutine overlap_save_iteration





subroutine mult_hermitian(X, H, N, Y)
  use iso_c_binding, only: c_double, c_long
  implicit none

  integer(kind=c_long), intent(in)    :: N
  real(kind=c_double),  intent(in)    :: X(0:N-1)
  real(kind=c_double),  intent(in)    :: H(0:N-1)
  real(kind=c_double),  intent(inout) :: Y(0:N-1)
  
  integer :: i, ireal, iimag, N_half

  N_half = ishft(N, -1)

  Y(0) = X(0) * H(0)
  Y(1) = X(1) * H(1)

  do i = 1, N_half - 1
     ireal = 2*i
     iimag = ireal + 1
     Y(ireal) = H(ireal)*X(ireal) - H(iimag)*X(iimag)
     Y(iimag) = H(iimag)*X(ireal) + H(ireal)*X(iimag)
  end do
end subroutine mult_hermitian




subroutine ifft_hermitian(Y, N)
  ! assuming Y(0) = Y_0, Y(1) = Y[N/2], Y(2*k) = Re(Y[k]), Y(2*k+1) = Im(y[k])
  use iso_c_binding, only: c_double, c_long
  implicit none

  integer(kind=c_long), intent(in) :: N
  real(kind=c_double),  intent(inout) :: Y(0:N-1)


  real(kind=c_double) :: Y_tmp(0:2*N-1)
  integer(kind=c_long) :: i, ire, iim, conjre, conjim, twoN, N_half

  twoN = ishft(N,1)
  N_half = ishft(N, -1)
  
  Y_tmp(0) = Y(0)/N
  Y_tmp(1) = 0

  Y_tmp(N) = Y(1)/N
  Y_tmp(N + 1) = 0

  do i = 1, N_half-1
     ire = ishft(i, 1)
     iim = ire + 1
     Y_tmp(ire) = Y(ire)/N
     Y_tmp(iim) = Y(iim)/N

     conjre = twoN - ire
     conjim = conjre + 1
     Y_tmp(conjre) = Y_tmp(ire)
     Y_tmp(conjim) = -Y_tmp(iim)
  end do

  call cfft(Y_tmp, N, 1)

  Y(0:N-1) = Y_tmp(0:twoN-1:2) ! output real part only
end subroutine ifft_hermitian
