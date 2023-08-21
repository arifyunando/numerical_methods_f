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
    !   Date             Programmer                   Description of changes
    !   ==========       ======================       =================================
    !   YYYY/MM/DD       A. Y. Sunanhadikusuma        ...
    !
    !
    use Bracketing
    use OpenMethods
    use PolyRoots
    implicit none

    !--- variables dictionary
        ! bracketing arguments
        real, external  :: polynomials                                      ! import function
        real            :: u_bound = 0.0, l_bound = -3.0, tol = 1.0E-8      ! bracket parameters
        integer         :: max_iter = 1000                                  ! maximum iteration
        logical         :: verbose = .false.                                ! print out iteration
        real            :: found_root                                       ! results

        ! function parameters
        real, allocatable, dimension(:) :: coefs    ! polynomial coefficients
        integer         :: order                    ! polynomial order

        ! placeholders
        integer :: i

    !--- process: call subroutines for each data input
        write(*, *) "enter the order of the polynomial function:"
        read(*, *)  order
        allocate (coefs(order))
        write(*, *) "enter coefficients in decending order:"
        do i = 1, order + 1
            read(*, *) coefs(i)
        end do

        call bisection(                     &
            polynomials, coefs, order,      &   ! function parameters
            u_bound, l_bound,               &   ! bracketing boundaries
            tol, max_iter, verbose,         &   ! iteration options
            found_root                      &   ! return value
        )
end program root


real function polynomials(x, coef, poly_order) result(rv)
    
    integer, intent(in) :: poly_order
    real, intent(in) :: x, coef(poly_order + 1)
    integer :: i

    rv = 0
    do i = 1, poly_order + 1
        rv = rv + coef(i)*x**(poly_order + 1 - i)
    end do

end function polynomials