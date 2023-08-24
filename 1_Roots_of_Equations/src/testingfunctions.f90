module TestingFunctions
    implicit none
    
contains
    
real function polynomialsExplicit(x, coef, size_coef) result(rv)
    integer, intent(in) :: size_coef
    real, intent(in) :: x, coef(size_coef)
    integer :: i

    rv = 0
    do i = 1, size_coef
        rv = rv + coef(i)*x**(size_coef - i)
    end do
end function polynomialsExplicit

real function polynomialsImplicit(x, coef) result(rv)
    real, intent(in) :: x, coef(:)
    integer :: coef_size
    integer :: i

    coef_size = size(coef)
    rv = 0
    do i = 1, coef_size
        rv = rv + coef(i)*x**(coef_size - i)
    end do
end function polynomialsImplicit

end module TestingFunctions