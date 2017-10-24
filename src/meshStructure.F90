!> Module that contains the data structure of a mesh associate to a
!! problem.
!! @author Diego T. Volpatto
        module meshStructure

            !> Data type for a mesh.
            type mesh

                integer :: numat !< Number of materials
                integer :: nsd !< Number of spatial
                integer :: nintp !< Number of integration points
                integer :: nnodes !< Number of nodes
                integer :: nelems !< Number of elements
                integer :: nen !< Number of element's nodes
                real*8, allocatable :: x(:,:) !< nodes coordinates
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
                character(50) :: geokind, & !< Coordinate frame kind
                                filename, & !< Mesh file
                                meshgen !< Mesh generator program
            end type

            contains

!****************************************************************************************

                !> Routine that allocate memory to node data.
                !! @param meshStrct     [in/out] mesh structure to allocate
                !! @param n     [in] number of nodes
                !! @author Diego T. Volpatto
                subroutine mallocNodes(meshStrct)

                    implicit none

                    type(mesh) :: meshStrct
                    integer :: n

                    n = meshStrct%nnodes

                    allocate(meshStrct%x(2,n+1)); meshStrct%x = 0.0d0
                    allocate(meshStrct%flagnode(n+1)); meshStrct%flagnode = 0

                endsubroutine

!****************************************************************************************

                !> Routine that allocate memory to element data.
                !! @param meshStrct     [in/out] mesh structure to allocate
                !! @param n    [in] number of elements
                !! @author Diego T. Volpatto
                subroutine mallocElem(meshStrct)

                    implicit none

                    type(mesh) :: meshStrct
                    integer :: n 

                    n = meshStrct%nelems

                    allocate(meshStrct%gnode(n+1,meshStrct%nen)); meshStrct%gnode = 0
                    allocate(meshStrct%ei(n+1)); meshStrct%ei = 0
                    allocate(meshStrct%ej(n+1)); meshStrct%ej = 0
                    allocate(meshStrct%ek(n+1)); meshStrct%ek = 0
                    allocate(meshStrct%si(n+1)); meshStrct%si = 0
                    allocate(meshStrct%sj(n+1)); meshStrct%sj = 0
                    allocate(meshStrct%sk(n+1)); meshStrct%sk = 0
                    allocate(meshStrct%xV(n+1)); meshStrct%xV = 0.0d0
                    allocate(meshStrct%yV(n+1)); meshStrct%yV = 0.0d0
                    allocate(meshStrct%mat(n+1)); meshStrct%mat = 0

                endsubroutine
        end module
