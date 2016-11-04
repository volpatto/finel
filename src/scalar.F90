        !> Contains variables and subroutine related to a general scalar
        !! problem.
        !! @author Diego T. Volpatto
        module mscalar

            implicit none

            contains

                !> Computes a master element contribution.
                !! @param x1    Lower bound coordinate of element
                !! @param x2    Upper bound coordinate of element
                !! @param n     Nodal points number of element
                !! @param ni    Order of integration rule
                !! @param mat   Material number
                subroutine localElem(x1,x2,n,ni,mat)

                    use mshapeFunctions, only: xi, w, shpf1d

                    implicit none

                    real*8 :: x1, x2
                    integer :: n, ni, mat

                    real*8 :: psi(n), dpsi(n), x, dx
                    real*8 :: eleft(n,n), eright(n)
                    real*8 :: xk, xb, xc, xf
                    integer :: i, j, l

                    xk = 1.d0; xb = 1.d0; xc = 0.d0
                    xf = 0.0

                    dx = (x2-x1)/2.d0

                    ! Initialize element arrays
                    eleft = 0.0d0; eright = 0.d0
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,ni
                        x = x1+(1.d0+xi(l,ni))*dx
                        call shpf1d(xi(l,ni),n,psi,dpsi)
                        do i=1,n
                            eright(i) = eright(i)+psi(i)*xf*w(l,ni)*dx
                            do j=1,n
                               eleft(i,j) = eleft(i,j) + & 
                                   (xk*dpsi(i)*dpsi(j)/dx/dx + &
                                   xc*psi(i)*dpsi(j)/dx + &
                                   xb*psi(i)*psi(j))*w(l,ni)*dx
                                write(*,*) "i =", i,"j =",j, "eleft =",eleft(i,j)
                            enddo
                            write(*,*) "i =", i, "eright =", eright(i)
                        enddo
                    enddo

                end subroutine

        end module
