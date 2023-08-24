module Bracketing
    implicit none
    
contains

subroutine bisection(                                                       &
        func, coef, u_bound, l_bound, tol, max_iter, verbose, root          &
    )
    !
    !
    !
    !
    !
    !
    !
    !
    !
    implicit none

    !--- arguments list
        real, external      :: func
        real, intent(in)    :: coef(:), u_bound, l_bound, tol
        logical, intent(in) :: verbose
        integer, intent(in) :: max_iter
        real, intent(out)   :: root
        
    !--- local variable list
        integer :: iter
        real    :: curr_u, curr_l            
        real    :: error, midpoint

    !--- process 1: initiatization
        if (verbose) write(*,*) ""
        if (verbose) write(*,*) "finding roots with 'bracketing method'..."
        if (verbose) write(*, 210)

        error = 1
        iter = 0
        curr_l = l_bound
        curr_u = u_bound

        ! check valid roots
        if (func(u_bound, coef, size(coef)) * func(l_bound, coef, size(coef)) >= 0) then
            print*, "ERROR: none or multiple roots exist inside the bounded domain"
            if (verbose) write(*, 210) 
            root = 0
            return
        end if

    !--- process 2: iteration
        do while ( (error > tol) .and. (iter < max_iter) )
            midpoint = (curr_u + curr_l) / 2
            if (func(midpoint, coef, size(coef)) * func(curr_l, coef, size(coef)) < 0) then
                error = abs(abs(midpoint - curr_u)/midpoint)
                curr_u = midpoint
            else 
                error = abs(abs(midpoint - curr_l)/midpoint)
                curr_l = midpoint
            end if
            iter = iter + 1
            if (verbose) write(*,200) iter, midpoint, error*100
        end do

    !--- process 3: closure
        if (verbose) write(*, 210) 
        if (verbose) write(*,*) "results:", midpoint 
        
    !--- formatting
        200 format ('  ', "Iteration ", I5, ':', F9.4, '    error: ', F6.2, ' %')
        210 format (":-----------------------------------------------:")
end subroutine bisection

subroutine regulaFalsi(                                                     &
        func, coef, u_bound, l_bound, tol, max_iter, verbose, root          &
    )
    !
    ! Purpose:
    !   Regula Falsi / False Point method is a bracketing method 
    !   that utilize linear function between two boundary points 
    !   to find the interesecting root. This method is supposedly
    !   better than the bracketing method whose nature is "brute-force" 
    !   and relatively inefficient
    !
    !   The equation to find the intersecting point 
    !   can be solved with this following equation.
    !   x_r = x_u - f(x_u)(x_l - x_u)/(f(x_l) - f(x_u))
    !
    !
    implicit none
    
    !--- argument list
        interface
            function func(x, args)
                real, intent(in) :: x                       ! evaluation point
                real, dimension(:), intent(in) :: args      ! function coefficients / parameters
                real :: func
            end function
        end interface
        real, intent(in)    :: coef(:), u_bound, l_bound, tol
        integer, intent(in) :: max_iter
        logical, intent(in) :: verbose
        real, intent(out)   :: root

    !--- local variable list
        real :: intersect_point
        real :: curr_l, curr_u
        real :: error
        integer :: iter

    !--- process 1: initialization
        if (verbose) write(*,*) ""
        if (verbose) write(*,*) "finding roots with 'false point method'..."
        if (verbose) write(*, 210)

        error = 1
        iter = 0
        curr_l = l_bound
        curr_u = u_bound

        ! check valid roots
        if (func(u_bound, coef) * func(l_bound, coef) > 0 ) then
            print*, "ERROR: none or multiple roots exist inside the bounded domain"
            if (verbose) write(*, 210) 
            root = 0
            return
        end if

    !--- process 2: iteration
        do while( (error > tol) .and. (iter < max_iter) )
            intersect_point = curr_u - func(curr_u, coef)*(curr_l - curr_u)         &
                                    / (func(curr_l, coef) - func(curr_u, coef))
            
            if (func(intersect_point, coef) * func(curr_l, coef) < 0) then
                error = abs(abs(intersect_point - curr_u)/intersect_point)
                curr_u = intersect_point
            else
                error = abs(abs(intersect_point - curr_l)/intersect_point)
                curr_l = intersect_point
            end if
            iter = iter + 1
            if (verbose) write(*,200) iter, intersect_point, error*100
        end do

    !--- process 3: closure
        if (verbose) write(*, 210) 
        if (verbose) write(*,*) "results:", intersect_point

    !--- formatting
        200 format ('  ', "Iteration ", I5, ':', F9.4, '    error: ', F6.2, ' %')
        210 format (":-----------------------------------------------:")
end subroutine regulaFalsi

subroutine regulaFalsiModified(                                             &
        func, coef, u_bound, l_bound, tol, max_iter, verbose, root          &
    )
    !--- argument list
        interface
            function func(x, args)
                real, intent(in) :: x                       ! evaluation point
                real, dimension(:), intent(in) :: args      ! function coefficients / parameters
                real :: func
            end function
        end interface
        real, intent(in)    :: coef(:), u_bound, l_bound, tol
        integer, intent(in) :: max_iter
        logical, intent(in) :: verbose
        real, intent(out)   :: root

    !--- local variable list
        real :: intersect_point
        real :: curr_l, curr_u, eval_l, eval_u
        real :: error
        integer :: iter, counter_u, counter_l

    !--- process 1: initialization
        if (verbose) write(*,*) ""
        if (verbose) write(*,*) "finding roots with 'modified false point method'..."
        if (verbose) write(*, 210)

        error = 1
        iter = 0
        curr_l = l_bound
        curr_u = u_bound
        eval_l = func(l_bound, coef)
        eval_u = func(u_bound, coef)

        ! check valid roots
        if (eval_u * eval_l > 0) then
            print*, "ERROR: none or multiple roots exist inside the bounded domain"
            if (verbose) write(*, 210) 
            root = 0
            return
        end if       

    !--- process 2: iteration
        do while( (error > tol) .and. (iter < max_iter) )
            intersect_point = curr_u - eval_u*(curr_l - curr_u)         &
                                    / (eval_l - eval_u)
            
            if (func(intersect_point, coef) * eval_l < 0) then
                error = abs(abs(intersect_point - curr_u)/intersect_point)
                curr_u = intersect_point
                eval_u = func(curr_u, coef)
                counter_u = 0
                counter_l = counter_l + 1
                if (counter_l > 2) then
                    eval_l = eval_l/2
                end if
            else
                error = abs(abs(intersect_point - curr_l)/intersect_point)
                curr_l = intersect_point
                eval_l = func(curr_l, coef)
                counter_l = 0
                counter_u = counter_u + 1
                if (counter_u > 2) then
                    eval_u = eval_u/2
                end if
            end if
            iter = iter + 1
            if (verbose) write(*,200) iter, intersect_point, error*100
        end do

    !--- process 3: closure
        if (verbose) write(*, 210) 
        if (verbose) write(*,*) "results:", intersect_point


    !--- formatting
        200 format ('  ', "Iteration ", I5, ':', F9.4, '    error: ', F6.2, ' %')
        210 format (":-----------------------------------------------:")
    
end subroutine regulaFalsiModified

end module Bracketing