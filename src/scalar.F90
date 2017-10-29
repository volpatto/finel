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
                    real*8 :: xk, xb, xc, xf, pi, v(2)
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

                    !print*, "mat = ", mesh_%mat(nel)
                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    !xk = 1.0d-8
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    !xb = 1.0d0
                    !xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    v(1) = scalar_%mat(mesh_%mat(nel),3); !print*, v(1)
                    v(2) = scalar_%mat(mesh_%mat(nel),4); !print*, v(2)
                    xf = 0.0d0

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
                        !xf=mesh_%yV(nel)
                        !print*, xf
                        if (scalar_%transient .eq. 0) then
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i)=scalar_%rhelem(i)+ &
                                        xf*psi(i)*fac
                        ! Compute stiffness matrix - lhs
                        ! I changed j=i,nen to j=1,nen (Diego, 23/10/2017)
                        do j=1,mesh_%nen
                        !scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           !fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           !xb*psi(i)*psi(j))
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           psi(i)*(v(1)*dpsix(j)+v(2)*dpsiy(j)) + &
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
                        do j=1,mesh_%nen
                        !scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           !scalar_%dt*fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           !xb*psi(i)*psi(j))
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           scalar_%dt*fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
                           psi(i)*(v(1)*dpsix(j)+v(2)*dpsiy(j)) + &
                           xb*psi(i)*psi(j)) + &
                           psi(i)*psi(j)*fac
                        enddo
                        enddo
                        endif
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

                !> Computes stabilizing SUPG terms.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine stabSUPG(mesh_, scalar_, nel)

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

                    real*8 :: xk, xb, xc, v(2), xf, pi, fac
                    integer :: i, j, l, i1, i2, k
                    integer :: inodes(mesh_%nen)
                    real*8 :: xnodes(mesh_%nsd,mesh_%nen)

                    ! Stabilizing SUPG variables
                    real*8 :: m_k, v_norm2, Pe_k, eps_k, tau_k, h

                    pi = 4.0d0*datan(1.0d0)

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    if (mesh_%nsd .eq. 1) xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    if (mesh_%nsd .eq. 2) then
                        v(1) = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                        v(2) = scalar_%mat(mesh_%mat(nel),4); !print*, xc
                    endif
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
                    v_norm2 = norm2(v) 
                    m_k = 1.d0/3.d0 ! For linear element
                    Pe_k = m_k*v_norm2*h/(2.d0*xk)
                    if (Pe_k .ge. 1) then
                        eps_k = 1.d0
                    else
                        eps_k = Pe_k
                    endif
                    tau_k = h*eps_k/(2.d0*v_norm2)

                    if (mesh_%nsd .eq. 1) then
                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        xf = xl
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        fac = w(l,mesh_%nintp)*dx
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            ! SUPG
                            ! Linear
                            scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                tau_k*xc*xf*dpsi(i)/dx*fac
                            if (mesh_%geokind .eq. "radial") then
                            scalar_%rhelem(i)=scalar_%rhelem(i)*xl
                            endif
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               ! SUPG
                               ! Linear
                                scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + & 
                                   (xc**2.d0*tau_k*dpsi(i)*dpsi(j)/dx/dx)*fac
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
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i) = scalar_%rhelem(i)+ &
                                xf*tau_k*(v(1)*dpsix(i)+v(2)*dpsiy(i))*fac
                        ! Compute stiffness matrix - lhs
                        do j=1,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           fac*(tau_k*(v(1)**2.d0*dpsix(i)*dpsix(j)+ v(2)**2.d0*dpsiy(i)*dpsiy(j) &
                           + v(1)*v(2)*dpsix(j)*dpsiy(i) + v(2)*v(1)*dpsiy(j)*dpsix(i)))
                        enddo
                        enddo
                    enddo
                    deallocate(xi2); deallocate(w2)
                endif

1234            format(1(1x,(a),1x,(i0)),4(1x,(a),1x,(f0.5)))
                end subroutine
                
                !> Computes stabilizing CAU terms to add in SUPG stabilizing.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine stabCAU(mesh_, scalar_, nel)

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

                    real*8 :: xk, xb, xc, v(2), xf, pi, fac
                    integer :: i, j, l, i1, i2, k
                    integer :: inodes(mesh_%nen)
                    real*8 :: xnodes(mesh_%nsd,mesh_%nen)
                    real*8 :: graduu, graduux, graduuy, uu, fuu

                    ! Stabilizing CAU variables
                    real*8 :: h, Res_phi, len_uv, len_grad, taup_s, eps_p, Pe_p, inner_ugrad
                    real*8 :: m_k, v_norm2, Pe_k, eps_k, tau_k, tau_c, p_order, alpha_c
                    real*8, parameter :: zero=1.d-30

                    pi = 4.0d0*datan(1.0d0)

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    if (mesh_%nsd .eq. 1) xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    if (mesh_%nsd .eq. 2) then
                        v(1) = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                        v(2) = scalar_%mat(mesh_%mat(nel),4); !print*, xc
                    endif
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
                    v_norm2 = norm2(v) 
                    m_k = 1.d0/3.d0 ! For linear element
                    Pe_k = m_k*v_norm2*h/(2.d0*xk)
                    if (Pe_k .ge. 1) then
                        eps_k = 1.d0
                    else
                        eps_k = Pe_k
                    endif
                    tau_k = h*eps_k/(2.d0*v_norm2)

                    if (mesh_%nsd .eq. 1) then
                    ! Begin integration loop
                    do l=1,mesh_%nintp
                        xl = x1+(1.d0+xi(l,mesh_%nintp))*dx
                        xf = xl
                        call shpf1d(xi(l,mesh_%nintp),mesh_%nen,psi,dpsi)
                        fac = w(l,mesh_%nintp)*dx
                        ! Compute previous iteration contribution
                        graduu = 0.d0
                        Res_phi = 0.d0
                        do i=1,mesh_%nen
                        graduu=graduu+dpsi(i)*scalar_%u_prev_it(inodes(i))
                        Res_phi = Res_phi + & 
                            xb*psi(i)*scalar_%u_prev_it(inodes(i))+ &
                            xc*dpsi(i)*scalar_%u_prev_it(inodes(i)) - & 
                            xf*psi(i)
                        enddo
                        if (graduu.lt.zero) goto 100
                        len_uv = dabs(Res_phi)/dabs(graduu)
                        if (len_uv.ge.dabs(xc)) then
                            tau_c = 0.d0
                        else
                            tau_c = tau_k*(dabs(xc)/len_uv-1.d0)
                        endif
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                            ! Compute stiffness matrix - lhs    
                            do j=1,mesh_%nen
                               ! CAU
                               ! Linear Element
                                scalar_%lhelem(i,j) = scalar_%lhelem(i,j) + &
                                    (tau_c*len_uv**2.d0*dpsi(j)*dpsi(i)/dx/dx)*fac
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
                        ! Compute previous iteration contribution
                        uu = 0.d0
                        graduux = 0.d0
                        graduuy = 0.d0
                        fuu = 0.d0
                        do i=1,mesh_%nen
                        uu=uu+psi(i)*scalar_%u_prev_it(inodes(i))
                        fuu=fuu+psi(i)*xf
                        graduux=graduux+dpsix(i)*scalar_%u_prev_it(inodes(i))
                        graduuy=graduuy+dpsiy(i)*scalar_%u_prev_it(inodes(i))
                        Res_phi = Res_phi + & 
                            xb*psi(i)*scalar_%u_prev_it(inodes(i))+ &
                            (v(1)*dpsix(i)*scalar_%u_prev_it(inodes(i))+&
                            v(2)*dpsiy(i)*scalar_%u_prev_it(inodes(i)))*& 
                            scalar_%u_prev_it(inodes(i)) - & 
                            xf*psi(i)
                        enddo
                        len_grad = dsqrt(graduux**2.d0+graduuy**2.d0)
                        if (len_grad.lt.zero) goto 100
                        Res_phi = xb*uu + &
                            (v(1)*graduux + v(2)*graduuy) - fuu
                        len_uv = dabs(Res_phi)/len_grad
                        ! ************** Stabilizing parameter *******************
                        p_order = 1.d0
                        Pe_p = v_norm2/(2.d0*xk*p_order)*h
                        if ((1.d0-1/Pe_p).lt.zero) then
                            eps_p = 0.d0
                        else
                            eps_p = 1.d0-1.d0/Pe_p
                        endif
                        taup_s = h*eps_p/(2.d0*v_norm2*p_order)
                        inner_ugrad = v(1)*graduux + v(2)*graduuy
                        if (inner_ugrad/Res_phi.lt.1.d0) then
                           alpha_c = 1.d0
                        else
                           alpha_c = inner_ugrad/Res_phi
                        endif
                        !if (len_uv.ge.v_norm2) then
                            !tau_c = 0.d0
                        !else
                            !tau_c = taup_s*(v_norm2/len_uv-1.d0)
                        !endif
                        if ((v_norm2/len_uv-alpha_c).lt.zero) then
                            tau_c = 0.d0
                        else
                            tau_c = taup_s*(v_norm2/len_uv-alpha_c)
                        endif
                        ! *******************************************************
                        !write(*,*) "Res_phi =", Res_phi
                        !write(*,*) "len_grad =", len_grad
                        !write(*,*) "len_uv =", len_uv
                        !write(*,*) "tau_c =", tau_c
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        ! Compute stiffness matrix - lhs
                        do j=1,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                            (fac*tau_c*len_uv**2.d0)*(dpsix(j)*dpsix(i)+dpsiy(j)+dpsiy(i))
                        enddo
                        enddo
                    enddo
                    deallocate(xi2); deallocate(w2)
100             endif

1234            format(1(1x,(a),1x,(i0)),4(1x,(a),1x,(f0.5)))
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
                    !xk = 1.d-8
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    !xb = 1.d0
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
                        do j=1,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           fac*(xb**2.d0*tau*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j)))
                        enddo
                        enddo
                    enddo
                    deallocate(xi2); deallocate(w2)
                endif

1234            format(1(1x,(a),1x,(i0)),4(1x,(a),1x,(f0.5)))
                end subroutine
                
                !> Computes a master element contribution in a Convection-Diffusion-Reaction 
                !! problem with a given velocity -- 2D.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @param nel      [in] Index of current element
                !! @author Diego Volpatto
                subroutine localElemCDR2D(mesh_, scalar_, scalar2_, nel)

                    use mshapeFunctions, only: xit, xiq, wt, wq, shpf2d
                    use meshStructure
                    use scalarStructure

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_, scalar2_
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
                    real*8 :: uup, vx, vy, gradpx, gradpy
                    real*8, parameter :: kappa = 1.0d0
                    real*8, parameter :: phi = 0.50

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

                    !xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    xk = 1.0d-2
                    !xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xb = 0.0d0
                    !xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xc = 0.0d0
                    xf = 0.0
                    ! Advection velocity components
                    vx = 0.d0; vy = 0.d0

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
                        xf=0.d0
                        !xf=mesh_%xV(nel)
                        if (scalar_%transient .eq. 0) then
                        ! Compute velocity by post-processing
                        gradpx = 0.d0
                        gradpy = 0.d0
                        do i=1,mesh_%nen
                        gradpx = gradpx + dpsix(i)*scalar2_%u(inodes(i))
                        gradpy = gradpy + dpsiy(i)*scalar2_%u(inodes(i))
                        enddo
                        vx = -kappa*gradpx; vy = -kappa*gradpy
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i)=scalar_%rhelem(i)+ &
                                        xf*psi(i)*fac
                        ! Compute stiffness matrix - lhs
                        do j=i,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
			   psi(i)*(vx*dpsix(j)+vy*dpsiy(j)) + &
                           xb*psi(i)*psi(j))
                        enddo
                        enddo
                        else
                        ! Compute previous solution contribution
                        uup = 0.d0
                        do i=1,mesh_%nen
                        uup = uup + psi(i)*scalar_%u_prev(inodes(i))
                        enddo
                        ! Compute velocity by post-processing
                        gradpx = 0.d0
                        gradpy = 0.d0
                        do i=1,mesh_%nen
                        gradpx = gradpx + dpsix(i)*scalar2_%u(inodes(i))
                        gradpy = gradpy + dpsiy(i)*scalar2_%u(inodes(i))
                        enddo
                        vx = -kappa*gradpx; vy = -kappa*gradpy
                        !write(*,*) vx, vy
                        ! Compute load vector - rhs
                        do i=1,mesh_%nen
                        scalar_%rhelem(i)=scalar_%rhelem(i)+ &
                                        scalar_%dt*xf*psi(i)*fac+ &
                                        fac*phi*psi(i)*uup
                        ! Compute stiffness matrix - lhs
                        do j=1,mesh_%nen
                        scalar_%lhelem(i,j)=scalar_%lhelem(i,j)+ &
                           scalar_%dt*fac*(xk*(dpsix(i)*dpsix(j)+dpsiy(i)*dpsiy(j))+&
			   psi(i)*(vx*dpsix(j)+vy*dpsiy(j)) + &
                           xb*psi(i)*psi(j)) + &
                           phi*psi(i)*psi(j)*fac
                        enddo
                        enddo
                        endif
                    enddo
                    !do i=1,mesh_%nen
                        !write(*,*) (scalar_%lhelem(i,j),j=1,mesh_%nen)
                    !enddo
                    deallocate(xi); deallocate(w)
                end subroutine

        end module
