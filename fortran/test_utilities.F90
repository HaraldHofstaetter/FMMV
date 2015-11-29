
#ifdef USE_SINGLE_PRECISION
module test_utilities_f
#else
module test_utilities
#endif

#ifdef USE_SINGLE_PRECISION
    integer, parameter :: tu_prec = 4
#else		
    integer, parameter :: tu_prec = 8
#endif	

    interface set_from_command_line
        module procedure set_from_command_line_real
        module procedure set_from_command_line_double
        module procedure set_from_command_line_integer
        module procedure set_from_command_line_logical
    end interface set_from_command_line

contains

    function relErr(x, x_ex) 
        implicit none
        real(tu_prec) :: x
        real(tu_prec) :: x_ex
        real(tu_prec) :: relErr
        relErr = abs((x-x_ex)/x_ex)
    end function relErr

    function relErr2(x, x_ex) 
        implicit none
        real(tu_prec) :: x(2)
        real(tu_prec) :: x_ex(2)
        real(tu_prec) :: relErr2
        relErr2 = sqrt(((x(1)-x_ex(1))*(x(1)-x_ex(1))+ &
                       (x(2)-x_ex(2))*(x(2)-x_ex(2))) &
                   /(x_ex(1)*x_ex(1)+x_ex(2)*x_ex(2)))
    end function relErr2

    subroutine relL2Err2(flag, x, x_ex, res)
        implicit none
        integer, intent(in) :: flag
        real(tu_prec), intent(in), optional :: x(2)
        real(tu_prec), intent(in), optional :: x_ex(2)
        real(tu_prec), intent(out), optional :: res

        real(tu_prec), save :: nom, den
        
        select case(flag)
            case(0)
                nom = 0.0d0
                den = 0.0d0
            case(1)
                nom = nom + (x(1)-x_ex(1))**2 + &
                            (x(2)-x_ex(2))**2 
                den = den + (x_ex(1))**2 + &
                            (x_ex(2))**2 
            case default
                res = sqrt(nom/den)
        end select 
    end subroutine relL2Err2

   
    function relErr3(x, x_ex) 
        implicit none
        real(tu_prec) :: x(3)
        real(tu_prec) :: x_ex(3)
        real(tu_prec) :: relErr3
        relErr3 = sqrt(((x(1)-x_ex(1))*(x(1)-x_ex(1))+ &
                       (x(2)-x_ex(2))*(x(2)-x_ex(2))+ &
                       (x(3)-x_ex(3))*(x(3)-x_ex(3))) &
                   /(x_ex(1)*x_ex(1)+x_ex(2)*x_ex(2)+x_ex(3)*x_ex(3)))
    end function relErr3

    subroutine relL2Err3(flag, x, x_ex, res)
        implicit none
        integer, intent(in) :: flag
        real(tu_prec), intent(in), optional :: x(3)
        real(tu_prec), intent(in), optional :: x_ex(3)
        real(tu_prec), intent(out), optional :: res

        real(tu_prec), save :: nom, den
        
        select case(flag)
            case(0)
                nom = 0.0d0
                den = 0.0d0
            case(1)
                nom = nom + (x(1)-x_ex(1))**2 + &
                            (x(2)-x_ex(2))**2 + &
                            (x(3)-x_ex(3))**2 
                den = den + (x_ex(1))**2 + &
                            (x_ex(2))**2 + &
                            (x_ex(3))**2
            case default
                res = sqrt(nom/den)
        end select 
    end subroutine relL2Err3
    
    subroutine relL2Err(flag, x, x_ex, res)
        implicit none
        integer :: flag
        real(tu_prec), intent(in), optional :: x
        real(tu_prec), intent(in), optional :: x_ex
        real(tu_prec), intent(out), optional :: res

        real(tu_prec), save :: nom, den
        
        select case(flag)
            case(0)  ! initialize
                nom = 0.0d0
                den = 0.0d0
            case(1) ! accumulate
                nom = nom + (x-x_ex)**2
                den = den + x_ex**2
            case default !finalize
                res = sqrt(nom/den)
        end select
    end subroutine relL2Err


subroutine set_from_command_line_double(varname, var)
    implicit none
    character(len=*), intent(in) :: varname
    real(kind=8), intent(inout) :: var

    integer :: k, pos
    integer, parameter :: buf_len = 256
    character(len=buf_len) :: buf
    real(kind=8)  :: var_old

    do k=1, iargc()
        call getarg(k, buf)
        pos= scan(buf, "=")
        if (buf(1:pos-1)==varname) then
            var_old = var
            read (buf(pos+1:buf_len),*) var
            write (*,*) "changed '", trim(varname), "' from",  var_old, "to", var
        end if
    end do    
end subroutine set_from_command_line_double

subroutine set_from_command_line_real(varname, var)
    implicit none
    character(len=*), intent(in) :: varname
    real(kind=4), intent(inout) :: var

    integer :: k, pos
    integer, parameter :: buf_len = 256
    character(len=buf_len) :: buf
    real(kind=4)  :: var_old

    do k=1, iargc()
        call getarg(k, buf)
        pos= scan(buf, "=")
        if (buf(1:pos-1)==varname) then
            var_old = var
            read (buf(pos+1:buf_len),*) var
            write (*,*) "changed '", trim(varname), "' from",  var_old, "to", var
        end if
    end do    
end subroutine set_from_command_line_real

    
subroutine set_from_command_line_integer(varname, var)
    implicit none
    character(len=*), intent(in) :: varname 
    integer, intent(inout) :: var 

    integer :: k, pos
    integer, parameter :: buf_len = 256
    character(len=buf_len) :: buf
    integer :: var_old

    do k=1, iargc()
        call getarg(k, buf)
        pos= scan(buf, "=")
        if (buf(1:pos-1)==varname) then
            var_old = var
            read (buf(pos+1:buf_len),*) var
            write (*,*) "changed '", trim(varname), "' from",  var_old, "to", var
        end if
    end do    
end subroutine set_from_command_line_integer

subroutine set_from_command_line_logical(varname, var)
    implicit none
    character(len=*), intent(in) :: varname
    logical, intent(inout) :: var

    integer :: k
    integer, parameter :: buf_len = 256
    character(len=buf_len) :: buf

    var = .false.
    do k=1, iargc()
        call getarg(k, buf)
        if(trim(buf)==trim(varname)) then
             var = .true.
             write (*,*) "option '", trim(varname), "'"
        end if
    end do    
end subroutine set_from_command_line_logical

#ifdef USE_SINGLE_PRECISION
end module test_utilities_f
#else
end module test_utilities
#endif
