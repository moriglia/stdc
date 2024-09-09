! SPDX-License-Identifier: GPL-3.0-or-later

subroutine filteriir(b, Nb, a, Na, x, Nx, y, y_initial)
  use iso_c_binding, only: c_double, c_long
  implicit none

  ! ----------------------
  ! --- i/o parameters ---
  ! ----------------------

  ! Numerator of Z transform (decreasing order)
  integer(kind=c_long), intent(in) :: Nb
  real(kind=c_double),  intent(in) :: b(Nb)
  
  ! Denominator coefficient of z transform (decreasing order)
  integer(kind=c_long), intent(in) :: Na
  real(kind=c_double),  intent(in) :: a(Na)

  ! Input waveform
  integer(kind=c_long), intent(in) :: Nx
  real(kind=c_double),  intent(in) :: x(Nx)

  ! Initial memory value
  real(kind=c_double), optional, intent(inout) :: y_initial(Na-1)

  ! Result
  real(kind=c_double), intent(inout) :: y(Nx) ! will be overwritten

  ! ------------------------
  ! --- other parameters ---
  ! ------------------------

  ! filter memory
  real(kind=c_double) :: y_memory(Na - 1)
  real(kind=c_double) :: x_memory(Nb)

  ! normalized coefficients
  real(kind=c_double) :: anorm(Na-1)
  real(kind=c_double) :: bnorm(Nb)

  ! iteration index
  integer :: n

  ! ------------------------
  ! --- setup simulation ---
  ! ------------------------

  ! setup memory
  if (.not. present(y_initial)) then
     y_memory(:) = 0
  else
     y_memory = y_initial
  end if

  x_memory(:Nb-1) = 0

  ! normalize coefficients
  if (a(1) .ne. 1.0) then
     bnorm = b/a(1)
     anorm = -a(2:)/a(1)
  else
     bnorm = b
     anorm = -a(2:)
  end if
  
  ! ------------------
  ! --- simulation ---
  ! ------------------
  do n = 1, Nx
     ! shift memories
     x_memory(2:) = x_memory(1:Nb-1)
     x_memory(1) = x(n)

     y_memory(2:) = y_memory(1:Na-2)
     y_memory(1) = y(n)

     ! update output
     y(n) = sum(bnorm*x_memory) + sum(anorm*y_memory)
  end do

  if (present(y_initial)) then
     y_initial = y_memory
  end if
  
end subroutine filteriir
