  ! SPDX-License-Identifier: GPL-3.0-or-later

subroutine filter_butter2_impulseinvariance_iir(w0T, x, Nx, y, y_initial)
  use iso_c_binding, only: c_double, c_long
  implicit none

  ! ----------------------
  ! --- i/o parameters ---
  ! ----------------------

  ! Normalized frequency
  real(kind=c_double), intent(in) :: w0T
  
  ! Input waveform
  integer(kind=c_long), intent(in) :: Nx
  real(kind=c_double),  intent(in) :: x(Nx)

  ! Initial memory value
  real(kind=c_double), optional, intent(inout) :: y_initial(2)

  ! Result
  real(kind=c_double), intent(inout) :: y(Nx) ! will be overwritten

  ! ------------------------
  ! --- other parameters ---
  ! ------------------------

  ! Parameters of Z transform
  real (kind=c_double) :: b(2)
  real (kind=c_double) :: a(3)

  b(1) = 0
  b(2) = sqrt(2.0)*w0T*sin(w0T/sqrt(2.0))*exp(-w0T/sqrt(2.0))

  a(1) = 1
  a(2) = -2*cos(w0T/sqrt(2.0))*exp(-w0T/sqrt(2.0))
  a(3) = exp(-w0T*sqrt(2.0))


  if (present(y_initial)) then
     call filteriir(b, 2, a, 3, x, Nx, y, y_initial)
  else
     y_initial(1) = 0
     y_initial(2) = 0
     call filteriir(b, 2, a, 3, x, Nx, y, y_initial)
  end if

end subroutine filter_butter2_impulseinvariance_iir
