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
                    real*8 :: x1, x2, xl2
                    real*8 :: xk, xb, xc, xf, pi, fac, uup
                    integer :: i, j, l, i1, i2
                    integer :: inodes(mesh_%nen)

                    pi = 4.0d0*datan(1.0d0) 

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xk = 1.d-8
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xb = 1.d0
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xf = 0.0

                    inodes=mesh_%gnode(nel,:)+1; !print*, inodes;
                    x1 = mesh_%x(1,inodes(1)); !print*, "x1=", x1
                    x2 = mesh_%x(1,inodes(2)); !print*,"x2=", x2

                    dx = (x2-x1)/2.d0; !print*, "dx =", dx;stop

                    ! Initialize element arrays
                    scalar_%lhelem = 0.0d0; scalar_%rhelem = 0.d0
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        !print*, xl
                        !xl = dx
                        !xf = dsin(pi*xl)*dcos(pi*xl)/xl
                        xf = xl
                        !xf = 6.d0*xl
                        !xf = (4.d0/3.d0)*xl
                        !xf = xl
                        !xf = 0.d0
                        !xf = xi(l,mesh_%nintp)
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        fac = w(l,mesh_%nintp)*dx
                        if (scalar_%transient .eq. 0) then
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                psi(i)*xf*fac
                            if (mesh_%geokind .eq. "radial") then
                            scalar_%rhelem(i)=scalar_%rhelem(i)*xl
                            endif
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xk*dpsi(i)*dpsi(j)/dx/dx + &
                                   xc*psi(i)*dpsi(j)/dx + &
                                   xb*psi(i)*psi(j))*fac
                                if (mesh_%geokind .eq. "radial") then
                                scalar_%lhelem(i,j)=scalar_%lhelem(i,j)*xl
                                endif
                            enddo
                        enddo
                        else
                        ! Compute previous solution contribution
                        uup = 0.d0
                        do i=1,mesh_%nen
                        uup = uup + psi(i)*scalar_%u_prev(inodes(i))
                        enddo
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                psi(i)*xf*fac*scalar_%dt+ &
                                uup*psi(i)*fac
                            if (mesh_%geokind .eq. "radial") then
                            scalar_%rhelem(i)=scalar_%rhelem(i)*xl
                            endif
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xk*dpsi(i)*dpsi(j)/dx/dx + &
                                   xc*psi(i)*dpsi(j)/dx + &
                                   xb*psi(i)*psi(j))*fac*scalar_%dt+ &
                                   psi(i)*psi(j)*fac
                                if (mesh_%geokind .eq. "radial") then
                                scalar_%lhelem(i,j)=scalar_%lhelem(i,j)*xl
                                endif
                            enddo
                        enddo
                        endif
                    enddo

                end subroutine

                !> Computes a master element contribution -- 2D.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine localElem2D(mesh_, scalar_, nel)

                    use mshapeFunctions, only: xit, xiq, wt, wq, shpf2d
                    use meshStructure
                    use scalarStructure

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    integer :: nel

                    real*8 :: psi(mesh_%nen), dpsi(mesh_%nsd,mesh_%nen)
                    real*8 :: dxds(2,2), dsdx(2,2), detJ, fac
                    real*8 :: dpsix(mesh_%nen), dpsiy(mesh_%nen)
                    real*8, allocatable :: xi(:,:), w(:)
                    !real*8 :: eleft(n,n), eright(n)
                    real*8 :: xl(mesh_%nsd)
                    real*8 :: xk, xb, xc, xf, pi
                    integer :: i, j, l, k
                    integer :: inodes(mesh_%nen)
                    real*8 :: xnodes(mesh_%nsd,mesh_%nen)
                    real*8 :: uup

                    pi = 4.0d0*datan(1.0d0)

                    ! Setting global nodes vector to local
                    inodes=mesh_%gnode(nel,:); !print*, inodes;i
                    if (mesh_%meshgen .eq. "easymesh") inodes=inodes+1 
                    do i=1,mesh_%nsd
                    do j=1,mesh_%nen
                    xnodes(i,j) = mesh_%x(i,inodes(j))
                    enddo
                    enddo

                    ! Setting quadrature points and weights due to
                    ! element geometry
                    if (mesh_%nen==3 .or. mesh_%nen==6) then 
                        allocate(xi(2,4)); xi=xit; 
                        allocate(w(4)); w=wt; 
                    endif
                    if (mesh_%nen==4 .or. mesh_%nen==8 .or.mesh_%nen==9) then 
                        allocate(xi(2,9)); xi=xiq; 
                        allocate(w(9)); w=wq; 
                    endif

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xk = 1.0d-8
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xb = 1.0d0
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xc = 0.0d0
                    xf = 0.0

                    ! Initialize element arrays
                    scalar_%lhelem = 0.0d0; scalar_%rhelem = 0.d0
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        ! Gauss point to integrate source/sink term
                        do i=1,mesh_%nsd
                        xl(i)=xi(i,l)
                        enddo
                        !xf=dsin(pi*mesh_%xV(nel))*dsin(pi*mesh_%yV(nel))
                        !xf = 1.0d0
                        call shpf2d(xi(mesh_%nsd,l),mesh_%nen,psi,dpsi)
                        ! Calculate DXDS
                        do i=1,2
                        do j=1,2
                        dxds(i,j) = 0.0d0
                        do k=1,mesh_%nen
                            dxds(i,j) = dxds(i,j)+dpsi(j,k)*xnodes(i,k)
                        enddo
                        enddo
                        enddo
                        ! Calculate DSDX
                        detJ = dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
                        if (detJ .le. 0.0d0) then
                            print*, "Bad Jacobian =", detJ, "Elem =",nel; stop
                        endif
                        dsdx(1,1)=dxds(2,2)/detJ
                        dsdx(2,2)=dxds(1,1)/detJ
                        dsdx(1,2)=-dxds(1,2)/detJ
                        dsdx(2,1)=-dxds(2,1)/detJ
                        ! Calculate d(psi)/dx
                        do i=1,mesh_%nen
                        dpsix(i)=dpsi(1,i)*dsdx(1,1)+dpsi(2,i)*dsdx(2,1)
                        dpsiy(i)=dpsi(1,i)*dsdx(1,2)+dpsi(2,i)*dsdx(2,2)
                        enddo
                        ! Accumulate integration point values of integrals
                        fac=detJ*w(l)
                        !xf=dsin(pi*mesh_%xV(nel))*dsin(pi*mesh_%yV(nel))
                        xf=1.d0
                        !xf=mesh_%xV(nel)
                        if (scalar_%transient .eq. 0) then
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i)=scalar_%rhelem(i)+ &
                                        xf*psi(i)*fac
                        ! Compute stiffness matrix - lhs
                        do j=i,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           xb*psi(i)*psi(j))
                        enddo
                        enddo
                        else 
                        ! Compute previous solution contribution
                        uup = 0.d0
                        do i=1,mesh_%nen
                        uup = uup + psi(i)*scalar_%u_prev(inodes(i))
                        enddo
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i)=scalar_%rhelem(i)+ &
                                        scalar_%dt*xf*psi(i)*fac+ &
                                        fac*psi(i)*uup
                        ! Compute stiffness matrix - lhs
                        do j=i,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           scalar_%dt*fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           xb*psi(i)*psi(j))
                        enddo
                        enddo
                        endif
                    enddo
                    ! Calculate lower symmetric part of EK
                    do i=1,mesh_%nen
                    do j=1,i
                    scalar_%lhelem(i,j) = scalar_%lhelem(j,i)
                    enddo
                    enddo
                    deallocate(xi); deallocate(w)
                end subroutine

                !> Computes a master element contribution in fracture
                !! case -- 1D.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine fracElem(mesh_, scalar_, nel, t)

                    use mshapeFunctions, only: xi, w, shpf1d
                    use meshStructure
                    use scalarStructure
                    use mtermo, only: Z_p

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    integer :: nel
                    real*8 :: t

                    real*8 :: psi(mesh_%nen), dpsi(mesh_%nen), xl, dx
                    !real*8 :: eleft(n,n), eright(n)
                    real*8 :: x1, x2, xl2, Ri, Ro
                    real*8 :: xk, xb, xc, xf, pi, fac, uup, uu, graduu
                    integer :: i, j, l, i1, i2, tmax
                    integer :: inodes(mesh_%nen)

                    ! Physical parameters declaration
                    real*8 :: dimL, p0, kappa, phi, alpha, beta, z, zp, mu

                    pi = 4.0d0*datan(1.0d0) 

                    ! Physical parameters
                    !dimL = 50.0
                    kappa = 1.d-14
                    mu = 1.d-5
                    phi = 0.25
                    p0 = 2.0d7
                    alpha = phi
                    beta = kappa/mu
                    !Ri = minval(mesh_%x(1,:)); 
                    Ro = maxval(mesh_%x(1,:))
                    Ri = mesh_%x(1,1)
                    tmax = maxval(scalar_%tprint)

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xk = beta
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xf = 0.0

                    inodes=mesh_%gnode(nel,:)+1; !print*, inodes;
                    x1 = mesh_%x(1,inodes(1)); !print*, "x1=", x1
                    x2 = mesh_%x(1,inodes(2)); !print*,"x2=", x2

                    dx = (x2-x1)/2.d0; !print*, "dx =", dx;stop

                    ! Initialize element arrays
                    scalar_%lhelem = 0.0d0; scalar_%rhelem = 0.d0
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        !print*, xl
                        xf = (2.0d9/(Ri-Ro))* &
                            (xl - Ro)* &
                            (dexp(-5.0d2*t/(tmax*scalar_%dt)))
                        !xf = 1.d0
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        fac = w(l,mesh_%nintp)*dx
                        if (scalar_%transient .eq. 0) then
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                psi(i)*xf*fac
                            if (mesh_%geokind .eq. "radial") then
                            scalar_%rhelem(i)=scalar_%rhelem(i)*xl
                            endif
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xk*dpsi(i)*dpsi(j)/dx/dx + &
                                   xc*psi(i)*dpsi(j)/dx + &
                                   xb*psi(i)*psi(j))*fac
                                if (mesh_%geokind .eq. "radial") then
                                scalar_%lhelem(i,j)=scalar_%lhelem(i,j)*xl
                                endif
                            enddo
                        enddo
                        else
                        ! Compute previous solution contribution
                        uup = 0.d0
                        uu = 0.d0
                        graduu = 0.d0
                        do i=1,mesh_%nen
                        uup = uup + psi(i)*scalar_%u_prev(inodes(i))
                        uu = uu + psi(i)*scalar_%u_prev_it(inodes(i))
                        graduu=graduu+dpsi(i)*scalar_%u_prev_it(inodes(i))
                        enddo
                        call Z_p(uu,z)
                        call Z_p(uup,zp)

                        !write(*,*) "uu =", uu
                        !write(*,*) "uup =", uup
                        !write(*,*) "z =", z
                        !write(*,*) "zp =", zp

                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                psi(i)*xf*fac*scalar_%dt+ &
                                alpha*(uup/zp)*psi(i)*fac
                            if (mesh_%geokind .eq. "radial") then
                            scalar_%rhelem(i)=scalar_%rhelem(i)*xl
                            endif
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xk*(uu/z)*dpsi(i)*dpsi(j)/dx/dx + &
                                   xc*psi(i)*dpsi(j)/dx + &
                                   xb*psi(i)*psi(j))*fac*scalar_%dt + &
                                   alpha*(psi(i)/z)*psi(j)*fac
                                if (mesh_%geokind .eq. "radial") then
                                scalar_%lhelem(i,j)=scalar_%lhelem(i,j)*xl
                                endif
                            enddo
                        enddo
                        endif
                    enddo

                end subroutine

                !> Computes stabilizing GGLS terms.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine stabGGLS(mesh_, scalar_, nel)

                    use mshapeFunctions, only: xi, w, shpf1d
                    use mshapeFunctions, only: xit, xiq, wt, wq, shpf2d
                    use meshStructure
                    use scalarStructure
                    use mUtilities, only: hk3_2d

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    integer :: nel

                    real*8 :: dpsi2(mesh_%nsd,mesh_%nen)
                    real*8 :: dxds(2,2), dsdx(2,2), detJ
                    real*8 :: dpsix(mesh_%nen), dpsiy(mesh_%nen)
                    real*8, allocatable :: xi2(:,:), w2(:)
                    real*8 :: psi(mesh_%nen), dpsi(mesh_%nen), xl, dx
                    real*8 :: x1, x2, xl2(mesh_%nsd)

                    real*8 :: xk, xb, xc, xf, pi, fac
                    integer :: i, j, l, i1, i2, k
                    integer :: inodes(mesh_%nen)
                    real*8 :: xnodes(mesh_%nsd,mesh_%nen)

                    ! Stabilizing GGLS variables
                    real*8 :: alpha, epst, tau, dxf, h

                    pi = 4.0d0*datan(1.0d0)

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xk = 1.d-8
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xb = 1.d0
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xf = 0.0

                    if (mesh_%nsd .eq. 1) then

                    inodes=mesh_%gnode(nel,:)+1; !print*, inodes;
                    x1 = mesh_%x(1,inodes(1)); !print*, "x1=", x1
                    x2 = mesh_%x(1,inodes(2)); !print*,"x2=", x2

                    dx = (x2-x1)/2.d0; !print*, "dx =", dx;stop
                    h = 2.d0*dx

                    else
                    ! Setting global nodes vector to local
                    inodes=mesh_%gnode(nel,:); !print*, inodes;
                    if (mesh_%meshgen .eq. "easymesh") inodes=inodes+1 
                    do i=1,mesh_%nsd
                    do j=1,mesh_%nen
                    xnodes(i,j) = mesh_%x(i,inodes(j))
                    enddo
                    enddo

                    call hk3_2d(xnodes,h)

                    ! Setting quadrature points and weights due to
                    ! element geometry
                    if (mesh_%nen==3 .or. mesh_%nen==6) then 
                        allocate(xi2(2,4)); xi2=xit; 
                        allocate(w2(4)); w2=wt; 
                    endif
                    if (mesh_%nen==4 .or. mesh_%nen==8 .or.mesh_%nen==9) then 
                        allocate(xi2(2,9)); xi2=xiq; 
                        allocate(w2(9)); w2=wq; 
                    endif
                    endif

                    ! Initialize element arrays
                    psi = 0.0d0; dpsi = 0.d0; dpsi2 = 0.d0

                    ! Stabilizing mesh parameters
                    alpha = (xb*(h**2.d0))/(6.d0*xk)
                    ! Exact eps form
                    !epst = (dcosh(dsqrt(6.d0*alpha))+2.d0)/ &
                        !(dcosh(dsqrt(6.d0*alpha))-1.d0) &
                        !- 1.d0/alpha
                    ! Asymptotic eps
                    if (alpha .ge. 8.d0) then 
                        epst = 1.d0
                    else if (alpha .ge. 1.d0) then
                        epst = 0.064d0*alpha+0.49
                    else
                        epst = 0.0d0
                    endif
                    tau = (epst*(h**2.0d0))/(6.d0*xb)
                    !write(456,1234) "nel =", nel, "alpha =", alpha,&
                            !"tau =", tau, "eps =", epst, "h =", h

                    if (mesh_%nsd .eq. 1) then
                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        xf = xl
                        dxf = 1.d0
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        fac = w(l,mesh_%nintp)*dx
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            ! Galerkin/Gradient-Least-Squares
                            ! Linear/Quad
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                dxf*tau*xb*dpsi(i)/dx*fac
                            if (mesh_%geokind .eq. "radial") then
                            scalar_%rhelem(i)=scalar_%rhelem(i)*xl
                            endif
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               ! Galerkin/Gradient-Least-Squares
                               ! Linear/Quad
                                scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xb**2.d0*tau*dpsi(i)*dpsi(j)/dx/dx)*fac
                                if (mesh_%geokind .eq. "radial") then
                                scalar_%lhelem(i,j)=scalar_%lhelem(i,j)*xl
                                endif
                            enddo
                        enddo
                    enddo
                
                else
                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        ! Gauss point to integrate source/sink term
                        do i=1,mesh_%nsd
                        xl2(i)=xi2(i,l)
                        enddo
                        call shpf2d(xi2(mesh_%nsd,l),mesh_%nen,psi,dpsi2)
                        ! Calculate DXDS
                        do i=1,2
                        do j=1,2
                        dxds(i,j) = 0.0d0
                        do k=1,mesh_%nen
                            dxds(i,j) = dxds(i,j)+dpsi2(j,k)*xnodes(i,k)
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
                        dpsix(i)=dpsi2(1,i)*dsdx(1,1)+dpsi2(2,i)*dsdx(2,1)
                        dpsiy(i)=dpsi2(1,i)*dsdx(1,2)+dpsi2(2,i)*dsdx(2,2)
                        enddo
                        ! Accumulate integration point values of integrals
                        fac=detJ*w2(l)
                        xf=1.d0
                        dxf=0.d0
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                dxf*tau*xb*(dpsix(i)+dpsiy(i))*fac
                        ! Compute stiffness matrix - lhs
                        do j=i,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           fac*(xb**2.d0*tau*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j)))
                        enddo
                        enddo
                    enddo
                    ! Calculate lower symmetric part of EK
                    do i=1,mesh_%nen
                    do j=1,i
                    scalar_%lhelem(i,j) = scalar_%lhelem(j,i)
                    enddo
                    enddo
                    deallocate(xi2); deallocate(w2)
                endif

1234            format(1(1x,(a),1x,(i0)),4(1x,(a),1x,(f0.5)))
                end subroutine

        end module
