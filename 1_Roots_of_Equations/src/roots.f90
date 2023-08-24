! "Numerical Methods Fortran" is a study repository that contains implementations
! of the numerical methods discussed in Chapra S. C. & Canale R. P. (2015). 
! Numerical methods for engineers (7th ed.). McGraw-Hill Higher Education.
! 
! © 2023 Arif Y Sunanhadikusuma


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
    use TestingFunctions
    implicit none

    !--- variables dictionary
        ! bracketing arguments
        real            :: u_bound = 0.0, l_bound = -3.0, tol = 1.0E-8      ! bracket parameters
        integer         :: max_iter = 1000                                  ! maximum iteration
        real            :: found_root                                       ! results
        logical         :: verbose = .true.

        ! function parameters
        real, allocatable, dimension(:) :: coefs    ! polynomial coefficients
        integer         :: order                    ! polynomial order

        ! placeholders
        integer :: i

    !--- process option 1: input data directly
        coefs = [1.0, 0.0, -4.0]
    
    !--- process option 2: input data from CLI
        ! write(*, *) "enter the order of the polynomial function:"
        ! read(*, *)  order
        ! allocate (coefs(order))
        ! write(*, *) "enter coefficients in decending order:"
        ! do i = 1, order + 1
        !     read(*, *) coefs(i)
        ! end do
    
    !--- process: call subroutines for solving roots
        print 300, size(coefs) - 1
        print 310, coefs

        call bisection(                     &
            polynomialsExplicit, coefs,     &   ! function parameters
            u_bound, l_bound,               &   ! bracketing boundaries
            tol, max_iter, verbose,         &   ! iteration options
            found_root                      &   ! return value
        )
        if (.not. verbose) print*, found_root
        
        call regulaFalsi(                   &
            polynomialsImplicit, coefs,     &   ! function parameters
            u_bound, l_bound,               &   ! bracketing boundaries
            tol, max_iter, verbose,         &   ! iteration options
            found_root                      &   ! return value
        )
        if (.not. verbose) print*, found_root

        call regulaFalsiModified(                   &
            polynomialsImplicit, coefs,     &   ! function parameters
            u_bound, l_bound,               &   ! bracketing boundaries
            tol, max_iter, verbose,         &   ! iteration options
            found_root                      &   ! return value
        )
        if (.not. verbose) print*, found_root

    !--- formatters
        300 format (' ', "Polynomial order :", I3)
        310 format (' ', "Coefficients     :", 3F5.1)
end program root


