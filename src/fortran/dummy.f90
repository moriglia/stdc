! SPDX-License-Identifier: GPL-3.0-or-later

subroutine dummy() bind(c)
  print *, 'Hello from fortran code'
end subroutine dummy


function mult_ints(a, b) bind(C) result(res)
  use iso_c_binding
  
  integer(kind=c_int), intent(in) :: a
  integer(kind=c_int), intent(in) :: b
  integer(kind=c_int) :: res
  
  res = a*b
end function mult_ints


function inner_prod(a,b) bind(C)
  use iso_c_binding, only: c_double
  implicit none

  real(kind=c_double), intent(in) :: a(3)
  real(kind=c_double), intent(in) :: b(3)
  real(kind=c_double) :: inner_prod


  inner_prod = sum(a*b)

end function inner_prod

  
subroutine print_array(x, N)
  use iso_c_binding, only: c_double, c_long

  integer(kind=c_long), intent(in) :: N
  real(kind=c_double), intent(in) :: x(N)

  print *, x(1:N)
end subroutine print_array
