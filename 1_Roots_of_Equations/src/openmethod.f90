! "Numerical Methods Fortran" is a study repository that contains implementations
! of the numerical methods discussed in Chapra S. C. & Canale R. P. (2015). 
! Numerical methods for engineers (7th ed.). McGraw-Hill Higher Education.
! 
! © 2023 Arif Y Sunanhadikusuma


module OpenMethods
    ! 
    ! Introduction:
    !   
    !   
    !   
    !   
    !   
    !     
    !     
    ! Record of revisions:
    !   Date             Programmer                   Description of changes
    !   ==========       ======================       =================================
    !   2023/08/26       A. Y. Sunanhadikusuma        Add fixed point interations method
    !   2023/08/29       A. Y. Sunanhadikusuma        Add Newton Raphson and Secant Method
    ! 
    ! 
    !     
    implicit none
    
contains

subroutine fixedPointInteration(                                            &                     
        func, coef, x_init, tol, max_iter, roots, errors, err_code          &
    )
    ! 
    ! Introduction:
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    !         
    implicit none

    !--- argument list
        interface
            function func(x, array)
                real, intent(in) :: x
                real, dimension(:), intent(in) :: array
                real :: func
            end function
        end interface

        integer, intent(in) :: max_iter
        real, intent(in)    :: x_init, tol, coef(:)

        real, allocatable, intent(out) :: roots(:), errors(:)
        integer, intent(out) :: err_code

    !--- local variable
        real :: root, error, x_current
        integer :: iter

    !--- process 1: initialization
        if (allocated(roots)) deallocate(roots)
        if (allocated(errors)) deallocate(errors)
        allocate(roots(max_iter), errors(max_iter))
        err_code = 0

        error = 1
        iter = 1
        x_current = x_init

    !--- process 2: iteration
        do while ( (error > tol) .and. (iter < max_iter) )
            root = func(x_current, coef)
            if (root .ne. 0.0) error = abs((root - x_current)/root)
            x_current = root

            roots(iter) = root
            errors(iter) = error
            iter = iter + 1
        end do
    
    !--- process 3: closure
        if (error > tol .or. isnan(error)) err_code = 2
        roots  = roots(:iter - 1)
        errors = errors(:iter - 1)
end subroutine fixedPointInteration

subroutine newtonRaphson(                                                   &                     
        func, dfunc, coef, x_init, tol, max_iter, roots, errors, err_code   &
    )
    ! 
    ! Introduction:
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    !         
    implicit none
    
    !--- argument list
        interface
            function func(x, array)
                real, intent(in) :: x
                real, dimension(:), intent(in) :: array
                real :: func
            end function
            function dfunc(x, array)
                real, intent(in) :: x
                real, dimension(:), intent(in) :: array
                real :: dfunc
            end function
        end interface

        integer, intent(in) :: max_iter
        real, intent(in)    :: x_init, tol, coef(:)

        real, allocatable, intent(out) :: roots(:), errors(:)
        integer, intent(out) :: err_code

    !--- local variable
        real :: x_current, x_new, error
        integer :: iter
    
    !--- process 1: initialization
        if (allocated(roots)) deallocate(roots)
        if (allocated(errors)) deallocate(errors)
        allocate(roots(max_iter), errors(max_iter))
        err_code = 0

        error = 1.0
        iter = 1
        x_current = x_init

    !--- process 2: iteration
        do while ( (error > tol) .and. (iter < max_iter) )
            x_new = x_current - func(x_current, coef)/dfunc(x_current, coef)
            if (x_new .ne. 0.0) error = abs((x_new - x_current)/x_new)
            x_current = x_new

            roots(iter) = x_new
            errors(iter) = error
            iter = iter + 1
        end do

    !--- process 3: closure
        if (error > tol .or. isnan(error)) err_code = 2
        roots  = roots(:iter - 1)
        errors = errors(:iter - 1)
end subroutine newtonRaphson

subroutine secantMethod(                                                    &                     
        func, coef, x_0, x_1, tol, max_iter, roots, errors, err_code        &
    )
    ! 
    ! Introduction:
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    !         
    implicit none
    
    !--- argument list
        interface
            function func(x, array)
                real, intent(in) :: x
                real, dimension(:), intent(in) :: array
                real :: func
            end function
        end interface

        integer, intent(in) :: max_iter
        real, intent(in)    :: x_0, x_1, tol, coef(:)

        real, allocatable, intent(out) :: roots(:), errors(:)
        integer, intent(out) :: err_code

    !--- local variable
        real :: x_curr, x_prev, x_next, error
        integer :: iter

    !--- process 1: initialization
        if (allocated(roots)) deallocate(roots)
        if (allocated(errors)) deallocate(errors)
        allocate(roots(max_iter), errors(max_iter))
        err_code = 0

        error = 1.0
        iter = 1
        x_curr = x_1
        x_prev = x_0
    
    !--- process 2: iteration
        do while ( (error > tol) .and. (iter < max_iter) )
            x_next = x_curr - func(x_curr, coef)*(x_prev - x_curr)/ &
                              (func(x_prev, coef) - func(x_curr, coef))
            if (x_next .ne. 0.0) error = abs((x_next - x_curr)/x_next)
            x_prev = x_curr
            x_curr = x_next

            roots(iter) = x_next
            errors(iter) = error
            iter = iter + 1
        end do

    !--- process 3: closure
        if (error > tol .or. isnan(error)) err_code = 2
        roots  = roots(:iter - 1)
        errors = errors(:iter - 1)
end subroutine secantMethod

subroutine secantMethodModified(                                            &                     
        func, coef, x_init, dx, tol, max_iter, roots, errors, err_code      &
    )
    ! 
    ! Introduction:
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    !         
    implicit none

    !--- argument list
        interface
            function func(x, array)
                real, intent(in) :: x
                real, dimension(:), intent(in) :: array
                real :: func
            end function
        end interface

        integer, intent(in) :: max_iter
        real, intent(in)    :: x_init, dx, tol, coef(:)

        real, allocatable, intent(out) :: roots(:), errors(:)
        integer, intent(out) :: err_code

    !--- local variable
        real :: x_curr, x_next, error
        integer :: iter

    !--- process 1: initialization
        if (allocated(roots)) deallocate(roots)
        if (allocated(errors)) deallocate(errors)
        allocate(roots(max_iter), errors(max_iter))
        err_code = 0

        error = 1.0
        iter = 1
        x_curr = x_init

    !--- process 2: iteration
        do while ( (error > tol) .and. (iter < max_iter) )
            x_next = x_curr - func(x_curr, coef)*dx/ &
                            (func(x_curr + dx, coef) - func(x_curr, coef))
            if (x_next .ne. 0.0) error = abs((x_next - x_curr)/x_next)
            x_curr = x_next

            roots(iter) = x_next
            errors(iter) = error
            iter = iter + 1
        end do

    !--- process 3: closure
        if (error > tol .or. isnan(error)) err_code = 2
        roots  = roots(:iter - 1)
        errors = errors(:iter - 1)
end subroutine secantMethodModified

subroutine brentMethod()
    ! 
    ! Introduction:
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    ! 
    !         
    implicit none
    
    !--- argument list

    !--- local variable

    !--- process

    !--- formatters

    
end subroutine brentMethod
    
end module OpenMethods