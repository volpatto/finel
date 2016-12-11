        !> Processor module to compute, assemble and solve the system.
        !! @author Diego T. Volpatto
        module mprocessor
            
            implicit none
            
            contains

            !> Form and assemble Ku = F system
            !! @param mesh_     [in/out] A mesh structure
            !! @param scalar_   [in/out] A scalar structure
            !! @param t         [in] Current simulation time
            !! @author Diego Volpatto
            subroutine formKF(mesh_, scalar_, t)

                use meshStructure
                use scalarStructure
                use mscalar,    only: localElem, localElem2D
                use mscalar,    only: fracElem

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                integer :: i1, i2, i3, nel
                real*8 :: t

                do nel = 1,mesh_%nelems

                    if (mesh_%nsd==1) then
                    !call localElem(mesh_, scalar_, nel)
                    call fracElem(mesh_, scalar_, nel, t)
                    else
                    call localElem2D(mesh_, scalar_, nel)
                    endif
                    call assmb(mesh_, scalar_, nel)

                enddo

            endsubroutine

            !> Assemble element stiffness matrix and load vector to
            !! global stiffness matrix and load vector, respectively.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @param nel       Index of current element
            !! @author Diego Volpatto
            subroutine assmb(mesh_, scalar_, nel)

                use meshStructure
                use scalarStructure

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                integer :: nel
                integer :: i, j, ig, jg

                do i=1,mesh_%nen

                   ig = mesh_%gnode(nel,i)+1 
                   
                   ! Assemble global load vector
                   scalar_%rhsys(ig) = scalar_%rhsys(ig) + &
                       scalar_%rhelem(i)

                   do j=1,mesh_%nen

                    jg = mesh_%gnode(nel,j)+1
                    
                    ! Assemble global stiffness matrix
                    scalar_%lhsys(ig,jg) = scalar_%lhsys(ig,jg) + &
                        scalar_%lhelem(i,j)

                    enddo

                 enddo

            endsubroutine

            !> Apply Dirichlet Boundary Condition.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @param n         Node index of BC
            !! @author Diego Volpatto
            subroutine drchlt(mesh_, scalar_, n)

                use meshStructure
                use scalarStructure

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                real*8 :: val
                integer :: n 

                integer :: i

                val = scalar_%vbc(mesh_%flagnode(n))

                do i=1,mesh_%nnodes

                ! Modify right hand side vector
                scalar_%rhsys(i) = scalar_%rhsys(i)-scalar_%lhsys(i,n)*val

                ! Modify row n and column n of stiffness matrix
                scalar_%lhsys(i,n) = 0.0d0
                scalar_%lhsys(n,i) = 0.0d0
                enddo

                scalar_%lhsys(n,n) = 1.0d0
                scalar_%rhsys(n) = val

            endsubroutine

            !> Apply Neumann Boundary Condition -- 1D.
            !! Prescribe -k(x)u = vbc.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @param n         Node index of BC
            !! @author Diego Volpatto
            subroutine neumann(mesh_, scalar_, n)

                use meshStructure
                use scalarStructure

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                integer :: n
                real*8 :: val, valmax, valmin

                val = scalar_%vbc(mesh_%flagnode(n))
                if (mesh_%geokind .eq. "rectangular") then
                    if (n.eq.1) scalar_%rhsys(n) = scalar_%rhsys(n) + val;
                    if (n.eq.mesh_%nnodes) scalar_%rhsys(n) = scalar_%rhsys(n) - val;
                else
                    valmax = maxval(mesh_%x)*val; 
                    valmin = minval(mesh_%x(1,1:))*val
                    if (n.eq.1) scalar_%rhsys(n) = scalar_%rhsys(n) +valmin;
                    if (n.eq.mesh_%nnodes) scalar_%rhsys(n) =scalar_%rhsys(n)-valmax;
                endif

            endsubroutine

            !> Modify Ku=F system to incorporate BC data.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @author Diego Volpatto
            subroutine applybc(mesh_, scalar_)

                use meshStructure
                use scalarStructure

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_

                integer :: i
                real*8 :: val

                ! Searching nodes with BC
                do i=1,mesh_%nnodes
                if (mesh_%flagnode(i) .ne. 0) then
                    !call drchlt(mesh_, scalar_,i, val(mesh_%flagnode(i)))
                    if (scalar_%kbc(mesh_%flagnode(i)).eq.1) then
                        call drchlt(mesh_, scalar_,i)
                    endif
                    if (scalar_%kbc(mesh_%flagnode(i)).eq.2) then
                        call neumann(mesh_, scalar_,i)
                    endif
                endif
                enddo

            endsubroutine

            !> Executes Gauss reduction and forward substitution solving
            !! Ku=F.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            subroutine solver(mesh_, scalar_)

                use meshStructure
                use scalarStructure
                use msolver, only: tri, rhsub

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_

                call tri(scalar_%lhsys,mesh_%nnodes);
                call rhsub(scalar_%lhsys,scalar_%u,scalar_%rhsys,mesh_%nnodes)

            endsubroutine

            !> Processor routine phase.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @author Diego T. Volpatto
            subroutine processor(mesh_, scalar_)

                use meshStructure
                use scalarStructure
                use mshapeFunctions, only: setint, setint2
                use mIO
                use mUtilities, only: check_conv, factor_picard
                use mUtilities, only: scaling_picard
                use mtermo, only: init_zcoef

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                
                real*8 :: t1, t2, t
                integer :: tstep, i, j
                integer, parameter :: niter = 5000
                logical :: flagit

                ! Picard parameters
                real*8 :: delta, omega, alpha, eps, omega_min, rho
                real*8 :: deltap

                if (mesh_%nsd==1) call setint
                if (mesh_%nsd==2) call setint2
                
                call init_zcoef

                ! Set the initial condition
                scalar_%u_prev = 2.0d7

                write(iout,*) "Processor procedures start:"
                write(*,*) "Processor procedures start:"
                write(iout,*)
                write(*,*)

                ! Begin loop in time
                i = 1
                do tstep = 1,scalar_%nsteps

                write(iout,*) "----------------------------------------"
                write(*,*) "----------------------------------------"
                write(iout,3333) "Time step = ", tstep, "t = ", scalar_%dt*tstep
                write(*,3333) "Time step =", tstep, "t =", scalar_%dt*tstep
                write(iout,*)
                write(*,*)

                t = scalar_%dt*tstep

                if (scalar_%linflag .eq. 1) then
                ! Linear problem mount and solver procedures

                call cpu_time(t1)
                call formKF(mesh_, scalar_,t)
                call applybc(mesh_, scalar_)
                call cpu_time(t2)
                write(iout,*) "Time elapsed (s) to mount system Ku=F:", (t2-t1)
                write(*,*) "Time elapsed (s) to mount system Ku=F:", (t2-t1)

                call cpu_time(t1)
                call solver(mesh_, scalar_)
                call cpu_time(t2)
                write(iout,*) "Time elapsed (s) to solve system Ku=F:", (t2-t1)
                write(*,*) "Time elapsed (s) to solve system Ku=F:", (t2-t1)
                
                ! Clear stiffness matrix and load vector
                scalar_%lhsys = 0.d0; scalar_%rhsys = 0.d0

                else
                ! Non-linear problem mount and solver by Picard
                ! iteration.

                ! Adaptative Picard method parameters initial values
                omega = 1.0d0
                omega_min = 0.5d0
                alpha = 0.9d0
                rho = 0.7d0
                eps = 1.d-5

                flagit = .false.
                scalar_%u_prev_it = scalar_%u_prev

                j = 0
                do while ((j.le.niter) .and. (flagit .eqv. .false.))
                
                j = j + 1
                write(iout,4444) "***** Picard iteration",j,"*****"
                write(*,4444) "***** Picard iteration",j,"*****"
                
                call cpu_time(t1)
                call formKF(mesh_, scalar_,t)
                call applybc(mesh_, scalar_)
                call cpu_time(t2)
                write(iout,*) "Time elapsed (s) to mount system Ku=F:", (t2-t1)
                write(*,*) "Time elapsed (s) to mount system Ku=F:", (t2-t1)
                
                call cpu_time(t1)
                call solver(mesh_, scalar_)
                call cpu_time(t2)
                write(iout,*) "Time elapsed (s) to solve system Ku=F:", (t2-t1)
                write(*,*) "Time elapsed (s) to solve system Ku=F:", (t2-t1)

                ! Check if iteration converged
                call check_conv(scalar_%u,scalar_%u_prev_it,mesh_%nnodes, &
                            eps, delta, flagit)
                !write(*,*) scalar_%u(1:5)
                !write(*,*) scalar_%u_prev_it(1:5)

                ! Adaptative Picard procedures
                ! ********************************************************

                ! Compute underrelaxation factor
                call factor_picard(alpha, delta, eps, omega_min, omega)

                ! Update previous iteration with underrelaxation factor
                write(*,5555) "Picard underrelaxation factor: ", omega
                write(iout,5555) "Picard underrelaxation factor: ", omega
                write(*,5555) "Relative error: ", delta
                write(iout,5555) "Relative error: ", delta
                scalar_%u_prev_it = omega*scalar_%u+(1.d0-omega)*scalar_%u_prev_it

                ! Rescaling shape factor for underrelaxation
                call scaling_picard(i, delta, deltap, eps, rho,&
                    omega, omega_min, alpha)

                ! Update maximum change error
                deltap = delta

                ! ********************************************************

                ! Clear K and F to reassemble in next iteration
                scalar_%lhsys = 0.d0; scalar_%rhsys = 0.d0;

                write(iout,*)
                write(*,*)

                enddo

                endif

                ! Record solution to a file
                if (tstep .eq. scalar_%tprint(i)) then
                    call print_files(mesh_, scalar_, i)
                    i = i + 1
                endif

                ! Update previous solution
                scalar_%u_prev = scalar_%u; scalar_%u = 0.d0

                enddo

                write(iout,*) "----------------------------------------"
                write(*,*) "----------------------------------------"
                write(iout,*)
                write(*,*)

3333            format(1x,(a,1x),i0,(1x,a),(f0.10))
4444            format(1x,(a,1x),i0,(1x,a))
5555            format(1x,(a,1x),f0.10)

            endsubroutine

        endmodule
