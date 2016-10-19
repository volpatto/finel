!> Module that contains the data structure of a general scalar problem
!! @author Diego T. Volpatto
        module scalarStruture

            !> Variables and characteristic data for a scalar problem.
            type scalarStructureSystem

                real*8, allocatable :: u(:) !< Solution vector
                !> Global left-hand system
                real*8, pointer :: lhsys(:,:)=>null()
                !> Global right-hand system
                real*8, pointer :: rhsys(:)=>null()

                integer :: numat !< Number of materials

            end type

        end module
