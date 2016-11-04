!> Module that contains the data structure of a mesh associate to a
!! problem.
!! @author Diego T. Volpatto
        module meshStruture

            !> Data type for a mesh.
            type mesh

                integer :: numat !< Number of materials
                integer :: nsd !< Number os spatial
                integer :: nnodes !< Number of nodes
                integer :: nelems !< Number of elements
                real*8, allocatable :: x(:), & !< x coordinates nodes
                                    y(:) !< y coordinates nodes
                integer*4, allocatable :: flagnode(:) !< boundary flag

                integer :: nelem !< Number of elements
                real*8, allocatable :: xV(:),& !< Circumcenter Elem xcoor 
                                    yV(:) !< Circumcenter Elem ycoor
                integer*4, allocatable :: gnode(:,:),& !< Global node
                                    mat(:) !< Element material kind
                integer*4, allocatable :: ei(:),& !< i-opposite element
                                    ej(:),& !< j-opposite element
                                    ek(:) !< k-opposite element
                integer*4, allocatable :: si(:),& !< Opposite i-side
                                    sj(:),& !< Opposite j-side
                                    sk(:) !< Opposite k-side
            end type

            contains

!****************************************************************************************

                !> Routine that allocate memory to node data.
                !! @param meshStrct     [in/out] mesh structure to allocate
                !! @param n     [in] number of nodes
                !! @author Diego T. Volpatto
                subroutine mallocNodes(meshStrct, n)

                    implicit none

                    type(mesh) :: meshStrct
                    integer :: n 

                    allocate(meshStrct%x(n)); meshStrct%x = 0.0d0
                    allocate(meshStrct%y(n)); meshStrct%y = 0.0d0
                    allocate(meshStrct%flagnode(n)); meshStrct%flagnode = 0

                endsubroutine

!****************************************************************************************

                !> Routine that allocate memory to element data.
                !! @param meshStrct     [in/out] mesh structure to allocate
                !! @param n    [in] number of elements
                !! @author Diego T. Volpatto
                subroutine mallocElem(meshStrct, n)

                    implicit none

                    type(mesh) :: meshStrct
                    integer :: n 

                    allocate(meshStrct%gnode(n,3)); meshStrct%gnode = 0
                    allocate(meshStrct%ei(n)); meshStrct%ei = 0
                    allocate(meshStrct%ej(n)); meshStrct%ej = 0
                    allocate(meshStrct%ek(n)); meshStrct%ek = 0
                    allocate(meshStrct%si(n)); meshStrct%si = 0
                    allocate(meshStrct%sj(n)); meshStrct%sj = 0
                    allocate(meshStrct%sk(n)); meshStrct%sk = 0
                    allocate(meshStrct%xV(n)); meshStrct%xV = 0.0d0
                    allocate(meshStrct%yV(n)); meshStrct%yV = 0.0d0
                    allocate(meshStrct%mat(n)); meshStrct%mat = 0

                endsubroutine
        end module
