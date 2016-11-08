        !> Contains variables and subroutine related to a general scalar
        !! problem.
        !! @author Diego T. Volpatto
        module mscalar

            implicit none

            contains

                !> Computes a master element contribution -- 1D.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine localElem(mesh_, scalar_, nel)
                !subroutine localElem1D(x1,x2,n,eleft,eright,ni,mat)

                    use mshapeFunctions, only: xi, w, shpf1d
                    use meshStructure
                    use scalarStructure

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    integer :: nel

                    real*8 :: psi(mesh_%nen), dpsi(mesh_%nen), xl, dx
                    !real*8 :: eleft(n,n), eright(n)
                    real*8 :: x1, x2
                    real*8 :: xk, xb, xc, xf, pi
                    integer :: i, j, l, i1, i2

                    pi = 4.0d0*datan(1.0d0) 

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xf = 0.0

                    i1 = mesh_%gnode(nel,1)+1; !print*, i1
                    i2 = mesh_%gnode(nel,2)+1; !print*, i2
                    x1 = mesh_%x(i1); !print*, "x1=", x1
                    x2 = mesh_%x(i2); !print*,"x2=", x2

                    dx = (x2-x1)/2.d0; !print*, "dx =", dx

                    ! Initialize element arrays
                    scalar_%lhelem = 0.0d0; scalar_%rhelem = 0.d0
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        !print*, xl
                        xf = dsin(pi*xl)*dcos(pi*xl)
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        do i=1,mesh_%nen
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                psi(i)*xf*w(l,mesh_%nintp)*dx
                            do j=1,mesh_%nen
                               scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xk*dpsi(i)*dpsi(j)/dx/dx + &
                                   xc*psi(i)*dpsi(j)/dx + &
                                   xb*psi(i)*psi(j))*w(l,mesh_%nintp)*dx
                                !write(*,*) "i =", i,"j =",j, &
                                    !"eleft=",scalar_%lhelem(i,j)
                            enddo
                            !write(*,*) "i =", i, "eright =",scalar_%rhelem(i)
                        enddo
                    enddo

                end subroutine

        end module
