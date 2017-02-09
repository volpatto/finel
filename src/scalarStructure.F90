!> Module that contains the data structure of a general scalar problem
!! @author Diego T. Volpatto
        module scalarStructure

            !> Variables and characteristic data for a scalar problem.
            type scalarStructureSystem

                real*8, allocatable :: u(:), & !< Solution vector
                                    u_prev(:), & !< Previous time step solution vector
                                    u_prev_it(:), & !< Previous non-linear solution vector
                                    lhelem(:,:), & !< Element left-hand system
                                    rhelem(:), & !< Element right-hand system
                                    lhsys(:,:), &!< Global left-hand system
                                    rhsys(:), & !< Global right-hand system
                                    vbc(:), & !< BC values vector
                                    mat(:,:) !< Material properties values 
                integer*4, allocatable :: kbc(:), & !< BC kind
                                       tprint(:) ! Time step index to print

                integer :: transient !< Transient's flag
                integer :: nsteps !< Number of time steps
                real*8 :: dt !< Time step
                integer :: linflag !< Flag to indicate if problem is linear or not
                integer :: stabm !< Flag to indicate which numerical method is employed

                real*8, allocatable :: grad(:), & !< Post-processed gradient
                                    incidence(:), & !< Node incidence to grad computation
                                    massflux(:) !< Mass flux over mesh nodes

            end type

            contains

            !> Routine to allocate and clear the Ku = F system.
            !! @param scalar_   [in/out] A general scalar structure
            !! @param n         [in] Number of global nodes
            subroutine mallocGlobalKF(scalar_,n)

                implicit none

                type(scalarStructureSystem) :: scalar_
                integer :: n

                allocate(scalar_%u(n)); scalar_%u = 0.0d0
                allocate(scalar_%grad(n)); scalar_%grad = 0.0d0
                allocate(scalar_%massflux(n)); scalar_%massflux = 0.0d0
                allocate(scalar_%incidence(n)); scalar_%incidence = 0.0d0
                if (scalar_%transient .eq. 1) then 
                allocate(scalar_%u_prev(n)); scalar_%u_prev = 0.0d0
                endif
                if (scalar_%linflag .eq. 0) then
                allocate(scalar_%u_prev_it(n)); scalar_%u_prev_it = 0.0d0
                endif
                allocate(scalar_%lhsys(n,n)); scalar_%lhsys = 0.0d0
                allocate(scalar_%rhsys(n)); scalar_%rhsys = 0.0d0

                print*, "Global KF memory allocated"

            endsubroutine

            !> Routine to allocate and clear the element KF.
            !! @param scalar_   [in/out] A general scalar structure
            !! @param n         [in] Number of element nodes
            subroutine mallocElemKF(scalar_,n)

                implicit none

                type(scalarStructureSystem) :: scalar_
                integer :: n

                allocate(scalar_%lhelem(n,n)); scalar_%lhelem = 0.0d0
                allocate(scalar_%rhelem(n)); scalar_%rhelem = 0.0d0

                print*, "Element KF memory allocated"
                write(*,*)

            endsubroutine

        end module
