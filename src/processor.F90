        !> Processor module to compute, assemble and solve the system.
        !! @author Diego T. Volpatto
        module mprocessor
            
            implicit none
            
            contains

            !> Form and assemble Ku = F system
            !! @param mesh_     [in/out] A mesh structure
            !! @param scalar_   [in/out] A scalar structure
            !! @author Diego Volpatto
            subroutine formKF(mesh_, scalar_)

                use meshStructure
                use scalarStructure
                use mscalar,    only: localElem 

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                integer :: i1, i2, i3, nel

                do nel = 1,mesh_%nelems

                    call localElem(mesh_, scalar_, nel)
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
            subroutine neumann1D(mesh_, scalar_, n)

                use meshStructure
                use scalarStructure

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                integer :: n
                real*8 :: val

                val = scalar_%vbc(mesh_%flagnode(n))
                if (n.eq.1) scalar_%rhsys(n) = scalar_%rhsys(n) + val;
                if (n.eq.mesh_%nnodes) scalar_%rhsys(n) = scalar_%rhsys(n) - val;

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
                        call neumann1D(mesh_, scalar_,i)
                    endif
                endif
                enddo

            endsubroutine

            subroutine solver(mesh_, scalar_)

                use meshStructure
                use scalarStructure
                use msolver, only: tri, rhsub

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_

                call tri(scalar_%lhsys,mesh_%nnodes)
                call rhsub(scalar_%lhsys,scalar_%u,scalar_%rhsys,mesh_%nnodes)

            endsubroutine

            !> Processor routine phase.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @author Diego T. Volpatto
            subroutine processor(mesh_, scalar_)

                use meshStructure
                use scalarStructure
                use mshapeFunctions, only: setint

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_

                call setint
                call formKF(mesh_, scalar_)
                call applybc(mesh_, scalar_)
                call solver(mesh_, scalar_)

            endsubroutine

        endmodule
