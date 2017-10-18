        !> Contains post-processing routines such as gradient
        !! computation and a-posteriori error routines.
        !! @author Diego T. Volpatto
        module mpostProc

            implicit none

            contains

                !> Post-processing gradient by node mean.
                !! case -- 1D.
                !! @param mesh_    [in/out] A mesh structure
                !! @param scalar_  [in/out] A scalar structure
                !! @author Diego Volpatto
                subroutine gradeval1D(mesh_, scalar_)

                    use mshapeFunctions, only: xi, shpf1d
                    use meshStructure
                    use scalarStructure
                    !use mtermo, only: Z_p

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    !integer :: nel
                    !real*8 :: t
                    integer :: nnodes, nelem

                    real*8 :: psi(mesh_%nen), dpsi(mesh_%nen), xl, dx
                    !real*8 :: eleft(n,n), eright(n)
                    real*8 :: x1, x2, xl2, Ri, Ro
                    real*8 :: xk, xb, xc, xf, pi, fac, uup, uu, graduu
                    integer :: i, j, l, i1, i2, tmax
                    integer :: inodes(mesh_%nen), nel
                    real*8 :: xnodes(mesh_%nsd, mesh_%nen)
                    real*8 :: incidence(mesh_%nnodes)

                    ! Physical parameters declaration
                    !real*8 :: dimL, p0, kappa, phi, alpha, beta, z, zp, mu

                    !pi = 4.0d0*datan(1.0d0) 

                    ! Physical parameters
                    !dimL = 50.0
                    !kappa = 1.d-14
                    !mu = 1.d-5
                    !phi = 0.25
                    !p0 = 2.0d7
                    !alpha = phi
                    !beta = kappa/mu
                    !Ri = minval(mesh_%x); Ro = maxval(mesh_%x)
                    !tmax = maxval(scalar_%tprint)

                    scalar_%incidence = 0.d0
                    scalar_%grad = 0.d0

                    do nel=1,mesh_%nelems

                    xk = scalar_%mat(mesh_%mat(nel),1); !print*, xk
                    !xk = beta
                    xb = scalar_%mat(mesh_%mat(nel),2); !print*, xb
                    xc = scalar_%mat(mesh_%mat(nel),3); !print*, xc
                    xf = 0.0

                    inodes=mesh_%gnode(nel,:)+1; !print*, inodes;
                    do i=1,mesh_%nsd
                    do j=1,mesh_%nen
                    xnodes(i,j) = mesh_%x(i,inodes(j))
                    enddo
                    enddo
                    psi = 0.0d0; dpsi = 0.d0

                    ! Begin integration loop
                    do l=1,mesh_%nen
                        call shpf1d(xi(l,mesh_%nen),mesh_%nen,psi,dpsi)
                        !call shpf1d(xnodes(1,l),mesh_%nen,psi,dpsi)
                        ! Compute previous solution contribution
                        uup = 0.d0
                        uu = 0.d0
                        graduu = 0.d0
                        do i=1,mesh_%nen
                        !uup = uup + psi(i)*scalar_%u_prev(inodes(i))
                        uu = uu + psi(i)*scalar_%u(inodes(i))
                        graduu=graduu+dpsi(i)*scalar_%u(inodes(i))
                        enddo
                        scalar_%grad(inodes(l)) = scalar_%grad(inodes(l)) - & 
                                    1.d0*graduu
                        scalar_%incidence(inodes(l)) = scalar_%incidence(inodes(l))+1.d0
                    enddo
                    enddo

                    do i=1,mesh_%nnodes
                    scalar_%grad(i)=scalar_%grad(i)/ & 
                            scalar_%incidence(i)
                    enddo

                end subroutine

                !> Computes mass flux over mesh nodes and interface
                !! fracture/well.
                !! @param mesh_     A mesh structure
                !! @param scalar_   A scalar structure
                !! @param CumulVol  Cumulative volume produced
                !! @param flPrev    Previous mass flux in boundary
                !! @author Diego Volpatto
                subroutine mflux1D(mesh_, scalar_, CumulVol, flPrev)

                    use meshStructure
                    use scalarStructure
                    use mtermo, only: Z_p

                    implicit none

                    !integer :: n, ni, mat
                    type(mesh) :: mesh_
                    type(scalarStructureSystem) :: scalar_
                    real*8 :: CumulVol, flPrev

                    real*8 :: p, z, gradp, MM, T, kappa, mu, Rgas, Ri
                    real*8 :: kB, nAv
                    integer :: i
                    real*8 :: pi, volInst, flInterface, InterfaceArea

                    pi = 4.0d0*datan(1.0d0)

                    ! Physical parameters
                    Ri = mesh_%x(1,1)
                    Ri = 1.d0
                    kappa = 1.d-14
                    mu = 1.d-5
                    kB = 1.381d-23
                    nAv = 6.0221415d23
                    Rgas = kB*nAv 
                    T = 352.55
                    MM = 16.01d-3

                    do i=1,mesh_%nnodes
                    p = scalar_%u(i)
                    call Z_p(p,z)
                    gradp = scalar_%grad(i)
                    if (mesh_%geokind .eq. "radial") then
                    scalar_%massflux(i) = -Ri*(MM/(Rgas*T))*(kappa/mu)*&
                        (p/z)*gradp
                    else
                    scalar_%massflux(i) = -(MM/(Rgas*T))*(kappa/mu)*&
                        (p/z)*gradp
                    endif
                    enddo

                    InterfaceArea = 2.d0*pi*Ri
                    flInterface = scalar_%massflux(1)*InterfaceArea
                    volInst = (flInterface + flPrev)/2.d0*scalar_%dt
                    CumulVol = CumulVol + volInst
                    flPrev = flInterface

                endsubroutine

        end module
