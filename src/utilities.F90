        !> Module for auxiliar routines.
        !! @author Diego T. Volpatto
        module mUtilities

            implicit none

            contains

!*************************************************************************************

                !> Generate points between x1 and x2 equally spaced in
                !! x(i). Same idea of numpy subroutine.
                !! @param x1         interval lower bound
                !! @param x2         interval upper bound
                !! @param nintv      num of intervals
                !! @param x          vector to assemble the values
                subroutine linspace(x1, x2, nintv, x)

                    implicit none

                    real*8  :: x1, x2, h
                    integer :: nintv, i
                    real*8, allocatable  :: x(:)

                    allocate(x(nintv+1)); x = 0.0d0
                    
                    h = (x2 - x1)/nintv
                    do i=1,nintv+1
                        x(i) = x1 + h*(i-1.0)
                    enddo

                end subroutine

!*************************************************************************************

                !> A function to test purpose.
                !! @param x     input coordinate
                real*8 function f1(x)
                    
                    implicit none
                    real*8 :: x

                    !f1 = x + x**2.0d0 + x**3.d0
                    !f1 = dsin(x)
                    f1 = 1.d0/dsqrt(x)
                    return  

                end function f1

!*************************************************************************************

                !> Subroutine that computes gaussian quadrature of f1.
                !! @param n         quadrature order
                !! @param x1        integral lower bound
                !! @param x2        integral upper bound
                subroutine quad1(n, x1, x2)

                    use mshapeFunctions, only: xi, w

                    implicit none

                    real*8 :: intval, x1, x2, dx, f, x
                    integer :: n, i

                    dx = (x2-x1)/2.d0
                    intval = 0.0d0
                    do i=1,n
                        x = x1 + (1.d0/2.d0)*(x2-x1)*(1.d0+xi(i,n))
                        f = f1(x)
                        intval = intval + f*w(i,n)*dx
                    end do
                    
                    write(*,*) "intval =", intval

                end subroutine

!*************************************************************************************

                !> Check if shpf1d works properly.
                !! @param n         element node numbers
                !! @param nelem     num of discrete intervals
                !! @param x         master element's coordinates
                subroutine test_shpf1d(n, nelem, x)

                    use mshapeFunctions, only: shpf1d

                    implicit none

                    integer :: n, nelem, i
                    real*8  :: x(nelem+1), shl(n), dshl(n)
                    real*8  :: a(nelem+1), b(nelem+1)

                    ! Clear
                    shl = 0.d0; dshl = 0.d0
                    a = 0.d0; b = 0.d0
                    do i=1,nelem+1
                        call shpf1d(x(i),n,shl,dshl)
                        a(i) = sum(shl)
                        b(i) = sum(dshl)
                        write(*,*) "a =", a(i), "b =", b(i)
                    enddo

                end subroutine

!*************************************************************************************

                !> Prints in the screen a matrix A(n,m)
                !! @param A     A matrix
                !! @param n     Number of lines of A
                !! @param m     Number of colunms of A
                !! @author      Diego Volpatto
                subroutine print_matrix(A, n, m)

                    implicit none

                    integer :: n, m
                    real*8 :: A(n,m)

                    integer :: i, j

                    do i=1,n
                    do j=1,m
                    write(*,*) i, j, A(i,j)
                    enddo
                    enddo

                endsubroutine

!*************************************************************************************

                !> Routine to check if non-linear iteration converged.
                !! @param u         Current solution vector
                !! @param uprev     Previous iteration solution vector
                !! @param nnodes    Number of nodes
                !! @param flag      Flag for convergence checked
                !! @author      Diego Volpatto
                subroutine check_conv(u, uprev, nnodes, flag)

                    implicit none

                    integer :: nnodes
                    real*8 :: u(nnodes), uprev(nnodes)
                    logical :: flag

                    integer :: i
                    real*8 :: norm, norm_up, norm_down, tol

                    tol = 1.d-4

                    norm_up = 0.d0; norm_down = 0.d0

                    do i=1,nnodes
                    norm_up = norm_up + (u(i)-uprev(i))**2.d0
                    norm_down = norm_down + (0.5d0*(u(i)+uprev(i)))**2.d0
                    enddo

                    norm_up=dsqrt(norm_up); norm_down=dsqrt(norm_down)
                    norm = norm_up/(1.d0+norm_down)
                    write(*,*) "Stop criteria norm: ", norm

                    if (norm .le. tol) then
                    flag = .true.;!stop
                    endif

                endsubroutine

        end module
