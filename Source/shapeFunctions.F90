        module mshapeFunctions
            
            implicit none

            real*8  ::  xi(4,4), w(4,4)

            contains
            
            !> Gauss quadrature data set routine.
            subroutine setint

            !... This routine defines the values of the parameters
            !... required for the numerical integration of element
            !... matrices and vectors.

                implicit none

                ! Gaussian quadrature order 1

                xi(1,1) = 0.0
                w(1,1)  = 2.0

                ! Gaussian quadrature order 2

                xi(1,2) = -1./dsqrt(3.0d0)
                xi(2,2) = -xi(1,2)
                w(1,2)  = 1.0
                w(2,2)  = w(1,2)

                ! Gaussian quadrature order 3

                xi(1,3) = -dsqrt(3.0d0/5.0d0)
                xi(2,3) = 0.0
                xi(3,3) = -xi(1,3)
                w(1,3)  = 5.0/9.0
                w(2,3)  = 8.0/9.0
                w(3,3)  = w(1,3)

                ! Gaussian quadrature order 4

                xi(1,4) = -0.8611363116
                xi(2,4) = -0.3399810436
                xi(3,4) = -xi(2,4)
                xi(4,4) = -xi(1,4)
                w(1,4)  = 0.3478548451
                w(2,4)  = 0.6521451549
                w(3,4)  = w(2,4)
                w(4,4)  = w(1,4)

            end subroutine

!****************************************************************************************
            
            !> Calculates the values of the shape functions and their
            !> derivatives.
            !! @param xl        specified value of master element coord
            !! @param n         number of element nodes
            !! @param psi       shape function values
            !! @param dpsi      derivatives shape functions values
            subroutine shpf1d(xl,n,psi,dpsi)

            !... Calculates the values of the shape functions psi and
            !... their derivatives dpsi with respect to the master
            !... element coordinate at a specified value xl.
                implicit none

                integer  ::  n
                real*8   ::  xl, psi(n), dpsi(n)

                ! Testing shape function order
                if (n .lt. 2 .or. n .gt. 4) then
                    write(*,*) "Error in call to shape", n
                    stop
                endif

                ! Clear
                psi = 0.d0; dpsi = 0.d0

                ! Linear shape functions
                if (n .eq. 2) then
                    psi(1) = 0.5*(1.0-xl)
                    psi(2) = 0.5*(1.0+xl)
                    dpsi(1)= -0.5
                    dpsi(2)= 0.5
                    return
                endif

                ! Quadratic shape functions
                if (n .eq. 3) then
                    psi(1) = xl*(xl-1.0)*0.5
                    psi(2) = 1.0-xl**2.0
                    psi(3) = xl*(xl+1.0)*0.5
                    dpsi(1)= xl-0.5
                    dpsi(2)= -2.0*xl
                    dpsi(3)= xl+0.5
                    return
                endif

                ! Cubic shape functions
                if (n .eq. 4) then
                    psi(1) = 9.0/16.0*(1.0/9.0-xl**2.0)*(xl-1.0)
                    psi(2) = 27.0/16.0*(1.0-xl**2.0)*(1.0/3.0-xl)
                    psi(3) = 27.0/16.0*(1.0-xl**2.0)*(1.0/3.0+xl)
                    psi(4) = -9.0/16.0*(1.0/9.0-xl**2.0)*(1.0+xl)
                    dpsi(1)= -9.0/16.0*(3.0*xl**2.0-2.0*xl-1.0/9.0)
                    dpsi(2)= 27.0/16.0*(3.0*xl**2.0-2.0/3.0*xl-1.0)
                    dpsi(3)= 27.0/16.0*(-3.0*xl**2.0-2.0/3.0*xl+1.0)
                    dpsi(4)= -9.0/16.0*(-3.0*xl**2.0-2.0*xl+1.0/9.0)
                    return
                endif

            end subroutine

        end module
