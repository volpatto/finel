        !> Module for setup phase by IO procedures.
        !! @author Diego T. Volpatto
        module msetup
            
            implicit none
            
            contains
            !> Reads parameters from input file.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            subroutine setupPhase(mesh_, scalar_)

                use mIO
                use mInputReader
                use meshStructure
                use scalarStructure

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                character(len=50) :: keyword_name, default_title_value
                integer :: i, j

                call readInputFileDS

                ! Reads case title
                keyword_name = "title"
                default_title_value = "unknown title"
                call readStringKeywordValue(keyword_name, title, default_title_value)
                write(iout,*) "Title problem: ", title
                write(iout,*)

                ! Reads number os spatial dimensions
                keyword_name = "nsd"
                call readIntegerKeywordValue(keyword_name,mesh_%nsd, 0_4)
                write(iout,*) "Spatial dimensions: ", mesh_%nsd
                write(iout,*)

                ! Reads integration order
                keyword_name = "nintp"
                call readIntegerKeywordValue(keyword_name,mesh_%nintp, 0_4)
                write(iout,*) "Number of gauss quadrature points: ",mesh_%nintp
                write(iout,*)

                ! Reads number of element's nodes
                keyword_name = "nen"
                call readIntegerKeywordValue(keyword_name,mesh_%nen, 0_4)
                write(iout,*) "Number of element's nodes:",mesh_%nen
                write(iout,*)

                ! Reads number of materials
                keyword_name = "numat"
                call readIntegerKeywordValue(keyword_name,mesh_%numat, 0_4)
                write(iout,*) "Number of materials:",mesh_%numat
                !write(iout,*)

                ! Reads material components properties
                keyword_name = "mat"
                call readRealMatrixValues(keyword_name,scalar_%mat, 0.0d0)
                do i=1,mesh_%numat
                write(iout,*) "Material:", i, (scalar_%mat(i,j), j=1,3)
                enddo
                write(iout,*)

                ! Reads boundary conditions data
                keyword_name = "vbc"
                call readBoundaryConditions(keyword_name,scalar_%kbc, &
                    scalar_%vbc, 0.0d0)

            endsubroutine

            !> Realizes preprocessor routines.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @author Diego T. Volpatto
            subroutine preprocessor(mesh_, scalar_)

                use meshStructure
                use scalarStructure
                use mIO

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_

                call setupPhase(mesh_,scalar_)
                call read_nodes(mesh_)
                call read_elems(mesh_)
                call mallocGlobalKF(scalar_, mesh_%nnodes)
                call mallocElemKF(scalar_, mesh_%nen)

            endsubroutine

        endmodule
