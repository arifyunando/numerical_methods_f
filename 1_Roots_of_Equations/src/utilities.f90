! "Numerical Methods Fortran" is a study repository that contains implementations
! of the numerical methods discussed in Chapra S. C. & Canale R. P. (2015). 
! Numerical methods for engineers (7th ed.). McGraw-Hill Higher Education.
! 
! © 2023 Arif Y Sunanhadikusuma


module Utilities
    implicit none
    
contains

!--- Testing Functions

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
    integer :: i, coef_size

    coef_size = size(coef)
    rv = 0
    do i = 1, coef_size
        rv = rv + coef(i)*x**(coef_size - i)
    end do
end function polynomialsImplicit

real function polynomialsFixedPoint(x, coef) result(rv)
    real, intent(in) :: x, coef(:)
    integer :: coef_size
    integer :: i
    
    ! check first order term exists
    coef_size = size(coef)

    if (coef(coef_size - 1) == 0.0) then
        rv = 0
        print*, "ERROR: polynomialsFixedPoint only accepts function with first order term"
        return
    end if

    ! evaluate return values
    rv = 0
    do i = 1, coef_size
        rv = rv + coef(i)*x**(coef_size - i)
    end do
    rv = - (rv - coef(coef_size - 1)*x) / coef(coef_size - 1)
end function polynomialsFixedPoint

real function polynomialsDerivative(x, coef) result(rv)
    real, intent(in) :: x, coef(:)
    integer :: i, coef_size

    coef_size = size(coef)
    rv = 0
    do i = 1, coef_size - 1
        rv = rv + (coef_size - i)*coef(i)*x**(coef_size - i - 1)
    end do
end function polynomialsDerivative

!--- Logging Subroutines

subroutine logResults(values, errors, err_code, method_name)
    !
    !
    !
    !
    !
    implicit none

    !--- argument list
        real, intent(in) :: values(:), errors(:)
        integer, intent(in) :: err_code
        character(*), intent(in) :: method_name
        integer :: i
        
        print*
        print*, "finding roots with ", method_name, "..."
        if ( err_code > 0 ) then
            if ( err_code == 1 ) then
                print*, "ERROR: none or multiple roots exist inside the bounded domain"
                return
            end if

            if ( err_code == 2 ) &
                print*, "WARNING: tolerance not met, showing final iteration value"
        end if
        print 210

        do i = 1, size(values)
           write(*,200) i, values(i), errors(i)*100
        end do
        
        print 210
        write(*,*) "results    :", values(size(values))
        write(*,*) "tolerance  :", errors(size(errors))
        

    !--- formatters
        200 format ('  ', "Iteration ", I5, ':', F9.4, '    error: ', F6.2, ' %')
        210 format (":-----------------------------------------------:")
    
end subroutine logResults

end module Utilities