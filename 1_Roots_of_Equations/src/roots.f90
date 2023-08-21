program root
    !
    ! Purpose:
    !   Roots of the equations are values are values that turn an equation 
    !   that is separated by an equal sign (=) into correct equality. 
    !   The codes for solving these equations is divided into 3 parts 
    !   which described in their own separate modules.
    !       1. Bracketing methods
    !       2. Open methods
    !       3. Roots of polynomials
    !
    ! Record of revisions:
    !   Date            Programmer                  Description of changes
    !   ==========      ======================      ============================
    !   YYYY/MM/DD      A. Y. Sunanhadikusuma       ...
    !
    !
    use Bracketing
    use OpenMethods
    use PolyRoots
    implicit none

    !--- variables dictionary
    ! bracketing values
    real, external  :: poly
    real            :: u_bound = 0.0, l_bound = -3.0
    real            :: tol = 1.0E-8
    integer         :: max_iter = 1000
    logical         :: verbose = .false.
    real            :: found_root

    ! function parameters
    integer         :: order
    real, allocatable, dimension(:) :: coefs

    ! placeholders
    integer :: i

    !--- processes
    write(*, *) "enter the order of the polynomial function:"
    read(*, *)  order
    allocate (coefs(order))
    write(*, *) "enter coefficients in decending order:"
    do i = 1, order + 1
        read(*, *) coefs(i)
    end do

    call bisection(             &
        poly, coefs, order,     &   ! function parameters
        u_bound, l_bound,       &   ! bracketing boundaries
        tol, max_iter, verbose, &   ! iteration options
        found_root              &   ! return value
    )
end program root


real function poly(x, coef, poly_order) result(rv)
    
    integer, intent(in) :: poly_order
    real, intent(in) :: x, coef(poly_order + 1)
    integer :: i

    rv = 0
    do i = 1, poly_order + 1
        rv = rv + coef(i)*x**(poly_order + 1 - i)
    end do    

end function poly