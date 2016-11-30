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
                        !xf = dsin(pi*xl)*dcos(pi*xl)
                        xf = 6.0d0*xl
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

                !> Computes a master element contribution -- 2D.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine localElem2D(mesh_, scalar_, nel)

                    use mshapeFunctions, only: xi, w, shpf1d
                    use meshStructure
                    use scalarStructure

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    integer :: nel

                    real*8 :: psi(mesh_%nen), dpsi(mesh_%nen), xl, dx
                    real*8 :: dxds(2,2), dsdx(2,2), detJ, fac
                    real*8 :: dpsix(mesh_%nen), dpsiy(mesh_%nen)
                    !real*8 :: eleft(n,n), eright(n)
                    real*8 :: x1, x2
                    real*8 :: xk, xb, xc, xf, pi
                    integer :: i, j, l, i1, i2, k

                    pi = 4.0d0*datan(1.0d0) 

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xf = 0.0

                    i1 = mesh_%gnode(nel,1)+1; !print*, i1
                    i2 = mesh_%gnode(nel,2)+1; !print*, i2
                    x1 = mesh_%x(i1); !print*, "x1=", x1
                    x2 = mesh_%x(i2); !print*,"x2=", x2

                    !dx = (x2-x1)/2.d0; !print*, "dx =", dx

                    ! Initialize element arrays
                    scalar_%lhelem = 0.0d0; scalar_%rhelem = 0.d0
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        !print*, xl
                        !xf = dsin(pi*xl)*dcos(pi*xl)
                        xf = 6.0d0*xl
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        ! Calculate DXDS
                        do i=1,2
                        do j=1,2
                        dxds(i,j) = 0.0d0
                        do k=1,mesh_%nen
                            dxds(i,j) = dxds(i,j)+dpsi(k,j)*x(i,k)
                        enddo
                        enddo
                        enddo
                        ! Calculate DSDX
                        detJ = dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
                        if (detJ .le. 0.0d0) then
                            print*, "Bad Jacobian =", detJ; stop
                        endif
                        dsdx(1,1)=dxds(2,2)/detJ
                        dsdx(2,2)=dxds(1,1)/detJ
                        dsdx(1,2)=-dxds(1,2)/detJ
                        dsdx(2,1)=-dxds(2,1)/detJ
                        ! Calculate d(psi)/dx
                        do i=1,mesh_%nen
                        dpsix(i)=dpsi(i,1)*dsdx(1,1)+dpsi(i,2)*dsdx(2,1)
                        dpsiy(i)=dpsi(i,1)*dsdx(1,2)+dpsi(i,2)*dsdx(2,2)
                        enddo
                        ! Accumulate integration point values of integrals
                        fac=detJ*w(l)
                        do i=1,mesh_%nen
                        scalar_%eright(i)=scalar_%eright(i)+ &
                                        xf*psi(i)*fac
                        do j=i,mesh_%nen
                        scalar_%eleft(i,j)=scalar_%eleft(i,j)+ &
                           fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           xb*psi(i)*psi(j))
                        enddo
                        enddo 
                    enddo
                    ! Calculate lower symmetric part of EK
                    do i=1,mesh_%nen
                    do j=1,i
                    scalar_%eright(i,j) = scalar_%eright(j,i)
                    enddo
                    enddo

                end subroutine

        end module
