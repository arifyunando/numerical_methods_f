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
    !   2023/08/20       A. Y. Sunanhadikusuma        Original code
    !   2023/08/24       A. Y. Sunanhadikusuma        Add calls to The Regula Falsi and 
    !                                                 Regula Falsi Modified method 
    !   2023/08/26       A. Y. Sunanhadikusuma        Add calls to the Open Methods
    !   2023/08/29       A. Y. Sunanhadikusuma        Add calls to Newton Raphson and Secant Method
    !
    use Bracketing
    use OpenMethods
    use PolyRoots
    use NonLinearSystem
    use Utilities
    implicit none

    !--- variables dictionary
        ! bracketing arguments
        real                :: u_bound = 0.0, l_bound = -3.0, tol = 1.0E-8      ! bracket parameters
        integer             :: max_iter = 1000, err_code                        ! maximum iteration
        real, allocatable   :: roots(:), errors(:)                              ! results

        ! function parameters
        real, allocatable, dimension(:) :: func_args

        ! placeholders

    !--- process 1 option 1: input data directly
        func_args = [1.0, 1.0, -2.0]
    
    !--- process 1 option 2: input data from CLI

    
    !--- process 2: call subroutines for solving roots
        print 200, size(func_args) - 1
        print 210, func_args

        
        ! Bracketing Methods

        call bisection(                             &
            polynomialsImplicit, func_args,         &   
            u_bound, l_bound,                       &   
            tol, max_iter,                          & 
            roots, errors, err_code                 &   
        )
        call logResults(roots, errors, err_code, "bisection")
        
        call regulaFalsi(                           &
            polynomialsImplicit, func_args,         &   
            u_bound, l_bound,                       &   
            tol, max_iter,                          & 
            roots, errors, err_code                 &   
        )
        call logResults(roots, errors, err_code, "false point method")

        call regulaFalsiModified(                   &
            polynomialsImplicit, func_args,         &   
            u_bound, l_bound,                       &   
            tol, max_iter,                          & 
            roots, errors, err_code                 &   
        )
        call logResults(roots, errors, err_code, "false point modified method")


        ! Open Methods

        call fixedPointInteration(                  &
            polynomialsFixedPoint, func_args,       &   
            u_bound,                                &    
            tol, max_iter,                          &   
            roots, errors, err_code                 &   
        )
        call logResults(roots, errors, err_code, "fixed point iterations")

        call newtonRaphson(                                         &
            polynomialsImplicit, polynomialsDerivative, func_args,  &   
            u_bound,                                                &    
            tol, max_iter,                                          &   
            roots, errors, err_code                                 &   
        )
        call logResults(roots, errors, err_code, "newton raphson method")

        call secantMethod(                          &
            polynomialsImplicit, func_args,         &   
            l_bound, u_bound,                       &    
            tol, max_iter,                          &   
            roots, errors, err_code                 &   
        )
        call logResults(roots, errors, err_code, "secant method")

        call secantMethodModified(                  &
            polynomialsImplicit, func_args,         &   
            l_bound, 0.01,                          &    
            tol, max_iter,                          &   
            roots, errors, err_code                 &   
        )
        call logResults(roots, errors, err_code, "secant method")

    !--- formatters
        200 format (' ', "Polynomial order :", I3)
        210 format (' ', "Coefficients     :", 3F5.1)
end program root
