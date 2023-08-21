module Bracketing
    implicit none
    
contains

    subroutine bisection(                                                       &
        func, coef, p_order, u_bound, l_bound, tol, max_iter, verbose, root     &
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
        real, intent(in)    :: u_bound, l_bound, tol
        logical, intent(in) :: verbose
        integer, intent(in) :: p_order, max_iter
        real, intent(in)    :: coef(:)
        real, intent(out)   :: root
        
        !--- local variable list
        integer :: iter
        real    :: curr_u_bound, curr_l_bound            
        real    :: error, midpoint

        !--- process
        error = 1
        iter = 0
        curr_l_bound = l_bound
        curr_u_bound = u_bound

        ! check valid roots
        if (func(u_bound, coef, p_order) * func(l_bound, coef, p_order) > 0) then
            print*, "ERROR: none or multiple roots exist inside the bounded domain"
            root = 0
            return
        end if
         
        if (verbose) write(*,*) "finding roots with 'bracketing method'..."
        if (verbose) write(*, 210)
        do while ( (error > tol) .and. (iter < max_iter))
            midpoint = (curr_u_bound + curr_l_bound) / 2
            if (func(midpoint, coef, p_order) * func(curr_l_bound, coef, p_order) < 0) then
                error = abs(abs(midpoint - curr_u_bound)/midpoint)
                curr_u_bound = midpoint
            else 
                error = abs(abs(midpoint - curr_l_bound)/midpoint)
                curr_l_bound = midpoint
            end if
            iter = iter + 1
            if (verbose) write(*,200) iter, midpoint, error*100
        end do

        if (verbose) write(*, 210) 
        write(*,*) "results:", midpoint 
        
        !--- formatting
        200 format ('  ', "Iteration ", I5, ':', F9.4, '    error: ', F6.2, ' %')
        210 format (":-----------------------------------------------:")
    end subroutine bisection

    subroutine false_point()
        !
        !
        !
        !
        !
        !
        !
        !
        write (*, *) "===================================="
        write (*, *) "        False Point Method          "
        write (*, *) "===================================="
        

        
    end subroutine false_point
end module Bracketing