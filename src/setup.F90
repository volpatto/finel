        !> Module for setup phase by IO procedures.
        !! @author Diego T. Volpatto
        module msetup
            
            implicit none
            
            contains
            !> Reads parameters from input file.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @author Diego T. Volpatto
            subroutine setupPhase(mesh_, scalar_)

                use mIO
                use mInputReader
                use meshStructure
                use scalarStructure
                use mshapeFunctions, only: quadext

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_
                character(len=50) :: keyword_name, default_title_value
                character(len=50) :: default_file_value,default_geo_value
                character(len=50) :: default_filename, default_meshgen
                integer :: i, j

                call readInputFileDS

                ! Reads case title
                keyword_name = "title"
                default_title_value = "unknown title"
                call readStringKeywordValue(keyword_name, title, default_title_value)
                write(iout,*) "Title problem: ", title
                write(iout,*)
                write(*,*) "Title problem: ", title
                write(*,*)

                ! Reads output solution extension
                keyword_name = "file_format"
                default_file_value = ".dat"
                call readStringKeywordValue(keyword_name,file_format,default_file_value)
                file_format = trim(file_format)
                write(iout,1111) "Writing solution in format ", file_format
                write(iout,*)
                write(*,1111) "Writing solution in format ", file_format
                write(*,*)

                ! Reads what program generate the mesh
                keyword_name = "mesh_gen"
                default_file_value = "easymesh"
                call readStringKeywordValue(keyword_name,mesh_%meshgen,&
                            default_file_value)
                mesh_%meshgen = trim(mesh_%meshgen)
                if ((mesh_%meshgen.ne."easymesh").and.(mesh_%meshgen.ne."triangle")) then
                    write(iout,*) "Incompatible mesh generator"
                    write(*,*) "Incompatible mesh generator"
                    write(*,*)
                    stop
                endif

                ! Reads number os spatial dimensions
                keyword_name = "nsd"
                call readIntegerKeywordValue(keyword_name,mesh_%nsd, 0_4)
                write(iout,*) "Spatial dimensions: ", mesh_%nsd
                write(iout,*)
                write(*,*) "Spatial dimensions: ", mesh_%nsd
                write(*,*)
                
                if (mesh_%nsd .eq. 1) then
                ! Reads geometry frame kind if one-dimensional
                keyword_name = "geometry"
                default_geo_value = "rectangular"
                call readStringKeywordValue(keyword_name,mesh_%geokind,default_geo_value)
                mesh_%geokind = trim(mesh_%geokind)
                write(iout,1111) "Coordinate frame system ", mesh_%geokind
                write(iout,*)
                write(*,1111) "Coordinate frame system ", mesh_%geokind
                write(*,*)
                endif

                ! Reads if stabilized methods are employed
                keyword_name = "stabilizing"
                call readIntegerKeywordValue(keyword_name,scalar_%stabm,0)
                write(iout,*) "Numerical method:"
                write(iout,*) "(0) Std Galerkin (1) GGLS (2) SUPG (3) CAU:",scalar_%stabm
                write(iout,*)
                write(*,*) "Numerical method:"
                write(*,*) "(0) Std Galerkin (1) GGLS (2) SUPG (3) CAU:",scalar_%stabm
                write(*,*)

                ! Reads mesh input file case name
                keyword_name = "filename"
                default_filename = "case1"
                call readStringKeywordValue(keyword_name,mesh_%filename,default_filename)
                write(iout,1111) "Mesh generator: ", mesh_%meshgen
                if (mesh_%meshgen .eq. "easymesh") then
                write(iout,2222) 'Nodes file: '//trim(mesh_%filename)//'.n'
                write(iout,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.e'
                endif
                if (mesh_%meshgen .eq. "triangle") then
                write(iout,2222) 'Nodes file: '//trim(mesh_%filename)//'.node'
                write(iout,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.ele'
                endif
                write(iout,*)
                write(*,1111) "Mesh generator: ", mesh_%meshgen
                if (mesh_%meshgen .eq. "easymesh") then
                write(*,2222) 'Nodes file: '//trim(mesh_%filename)//'.n'
                write(*,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.e'
                endif
                if (mesh_%meshgen .eq. "triangle") then
                write(*,2222) 'Nodes file: '//trim(mesh_%filename)//'.node'
                write(*,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.ele'
                endif
                write(*,*); !stop

                ! Reads quadrature kind
                keyword_name = "quadext"
                call readIntegerKeywordValue(keyword_name,quadext, 0)
                write(iout,*) "Gaussian quadrature (0=normal,1=extended): ",quadext
                write(iout,*)
                write(*,*) "Gaussian quadrature (0=normal,1=extended): ",quadext
                write(*,*)

                ! Reads integration order
                keyword_name = "nintp"
                call readIntegerKeywordValue(keyword_name,mesh_%nintp, 0_4)
                write(iout,*) "Number of gauss quadrature points: ",mesh_%nintp
                write(iout,*)
                write(*,*) "Number of gauss quadrature points: ",mesh_%nintp
                write(*,*)

                ! Reads number of element's nodes
                keyword_name = "nen"
                call readIntegerKeywordValue(keyword_name,mesh_%nen, 0_4)
                write(iout,*) "Number of element's nodes:",mesh_%nen
                write(iout,*)
                write(*,*) "Number of element's nodes:",mesh_%nen
                write(*,*)

                ! Reads number of element's nodes
                keyword_name = "linflag"
                call readIntegerKeywordValue(keyword_name,scalar_%linflag, 1)
                write(iout,*) "(1) Linear analysis (0) Nonlinear:",scalar_%linflag
                write(iout,*)
                write(*,*) "(1) Linear analysis (0) Nonlinear:",scalar_%linflag
                write(*,*)

                ! Reads number of materials
                keyword_name = "numat"
                call readIntegerKeywordValue(keyword_name,mesh_%numat, 0_4)
                write(iout,*) "Number of materials:",mesh_%numat
                write(*,*) "Number of materials:",mesh_%numat
                !write(iout,*)

                ! Reads material components properties
                keyword_name = "mat"
                call readRealMatrixValues(keyword_name,scalar_%mat, 0.0d0)
                do i=1,mesh_%numat
                write(iout,*) "Material:", i, (scalar_%mat(i,j), j=1,3)
                write(*,*) "Material:", i, (scalar_%mat(i,j), j=1,3)
                enddo
                write(iout,*)
                write(*,*)

                ! Reads boundary conditions constant data
                keyword_name = "vbc"
                call readBoundaryConditions(keyword_name,scalar_%kbc, &
                    scalar_%vbc, 0.0d0)

                ! Reads transient/stationary mode
                keyword_name = "transient"
                call readIntegerKeywordValue(keyword_name,scalar_%transient, 0)
                write(iout,*) "Transient problem (0=no; 1=yes):",scalar_%transient
                write(*,*) "Transient problem (0=no; 1=yes):",scalar_%transient
                !write(iout,*) "Transient problem (0=no; 1=yes):",scalar_%transient
                !write(iout,*)

                if (scalar_%transient .eq. 1) then
                ! Reads time step value
                keyword_name = "dt"
                call readRealKeywordValue(keyword_name,scalar_%dt, 0.0d0)
                write(iout,*) "Time step:",scalar_%dt
                write(*,*) "Time step:",scalar_%dt
                !write(iout,*) "Time step:",scalar_%dt
                !write(iout,*)

                ! Reads number of time steps
                keyword_name = "nsteps"
                call readIntegerKeywordValue(keyword_name,scalar_%nsteps, 1)
                write(iout,*) "Number of time's steps:",scalar_%nsteps
                write(iout,*) "Total time:",scalar_%nsteps*scalar_%dt
                write(*,*) "Number of time's steps:",scalar_%nsteps
                write(*,*) "Total time:",scalar_%nsteps*scalar_%dt

                ! Reads time steps index that solution will be printed
                keyword_name = "tprint"
                call readIntegerArrayValues(keyword_name,scalar_%tprint, 0)
                !write(iout,*) "Time steps which solution will be record"
                write(iout,*) "Time steps of solution that will be record"
                write(*,*) "Time steps of solution that will be record"
                do i=1,size(scalar_%tprint)
                write(iout,*) scalar_%tprint(i)
                write(*,*) scalar_%tprint(i)
                enddo
                else
                    scalar_%nsteps = 1
                    write(iout,*) "Number of time's steps:",scalar_%nsteps
                    write(*,*) "Number of time's steps:",scalar_%nsteps
                    allocate(scalar_%tprint(1)); scalar_%tprint = 1
                    write(*,*) "Time steps of solution that will be record"
                    do i=1,size(scalar_%tprint)
                    write(iout,*) scalar_%tprint(i)
                    write(*,*) scalar_%tprint(i)
                    enddo
                endif
                write(iout,*);
                write(*,*);

1111            format(1x,2(a))
2222            format(1x,(a))

            endsubroutine

            !> Reads parameters from input file for Tracer problem.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @param scalar2_   Another scalar structure
            !! @author Diego T. Volpatto
            subroutine setupPhaseTracer(mesh_, scalar_, scalar2_)

                use mIO
                use mInputReader
                use meshStructure
                use scalarStructure
                use mshapeFunctions, only: quadext

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_, scalar2_
                character(len=50) :: keyword_name, default_title_value
                character(len=50) :: default_file_value,default_geo_value
                character(len=50) :: default_filename, default_meshgen
                integer :: i, j

                call readInputFileDS

                ! Reads case title
                keyword_name = "title"
                default_title_value = "unknown title"
                call readStringKeywordValue(keyword_name, title, default_title_value)
                write(iout,*) "Title problem: ", title
                write(iout,*)
                write(*,*) "Title problem: ", title
                write(*,*)

                ! Reads output solution extension
                keyword_name = "file_format"
                default_file_value = ".dat"
                call readStringKeywordValue(keyword_name,file_format,default_file_value)
                file_format = trim(file_format)
                write(iout,1111) "Writing solution in format ", file_format
                write(iout,*)
                write(*,1111) "Writing solution in format ", file_format
                write(*,*)

                ! Reads what program generate the mesh
                keyword_name = "mesh_gen"
                default_file_value = "easymesh"
                call readStringKeywordValue(keyword_name,mesh_%meshgen,&
                            default_file_value)
                mesh_%meshgen = trim(mesh_%meshgen)
                if ((mesh_%meshgen.ne."easymesh").and.(mesh_%meshgen.ne."triangle")) then
                    write(iout,*) "Incompatible mesh generator"
                    write(*,*) "Incompatible mesh generator"
                    write(*,*)
                    stop
                endif

                ! Reads number os spatial dimensions
                keyword_name = "nsd"
                call readIntegerKeywordValue(keyword_name,mesh_%nsd, 0_4)
                write(iout,*) "Spatial dimensions: ", mesh_%nsd
                write(iout,*)
                write(*,*) "Spatial dimensions: ", mesh_%nsd
                write(*,*)
                
                if (mesh_%nsd .eq. 1) then
                ! Reads geometry frame kind if one-dimensional
                keyword_name = "geometry"
                default_geo_value = "rectangular"
                call readStringKeywordValue(keyword_name,mesh_%geokind,default_geo_value)
                mesh_%geokind = trim(mesh_%geokind)
                write(iout,1111) "Coordinate frame system ", mesh_%geokind
                write(iout,*)
                write(*,1111) "Coordinate frame system ", mesh_%geokind
                write(*,*)
                endif

                ! Reads if stabilized methods are employed
                keyword_name = "stabilizing"
                scalar_%stabm = 0
                call readIntegerKeywordValue(keyword_name,scalar2_%stabm,0)
                write(iout,*) "Numerical method for tracer:"
                write(iout,*) "(0) Standard Galerkin (1) GGLS:",scalar2_%stabm
                write(iout,*)
                write(*,*) "Numerical method for tracer:"
                write(*,*) "(0) Standard Galerkin (1) GGLS:",scalar2_%stabm
                write(*,*)

                ! Reads mesh input file case name
                keyword_name = "filename"
                default_filename = "case1"
                call readStringKeywordValue(keyword_name,mesh_%filename,default_filename)
                write(iout,1111) "Mesh generator: ", mesh_%meshgen
                if (mesh_%meshgen .eq. "easymesh") then
                write(iout,2222) 'Nodes file: '//trim(mesh_%filename)//'.n'
                write(iout,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.e'
                endif
                if (mesh_%meshgen .eq. "triangle") then
                write(iout,2222) 'Nodes file: '//trim(mesh_%filename)//'.node'
                write(iout,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.ele'
                endif
                write(iout,*)
                write(*,1111) "Mesh generator: ", mesh_%meshgen
                if (mesh_%meshgen .eq. "easymesh") then
                write(*,2222) 'Nodes file: '//trim(mesh_%filename)//'.n'
                write(*,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.e'
                endif
                if (mesh_%meshgen .eq. "triangle") then
                write(*,2222) 'Nodes file: '//trim(mesh_%filename)//'.node'
                write(*,2222) 'Element and conectivity file: '//trim(mesh_%filename)//'.ele'
                endif
                write(*,*); !stop

                ! Reads quadrature kind
                keyword_name = "quadext"
                call readIntegerKeywordValue(keyword_name,quadext, 0)
                write(iout,*) "Gaussian quadrature (0=normal,1=extended): ",quadext
                write(iout,*)
                write(*,*) "Gaussian quadrature (0=normal,1=extended): ",quadext
                write(*,*)

                ! Reads integration order
                keyword_name = "nintp"
                call readIntegerKeywordValue(keyword_name,mesh_%nintp, 0_4)
                write(iout,*) "Number of gauss quadrature points: ",mesh_%nintp
                write(iout,*)
                write(*,*) "Number of gauss quadrature points: ",mesh_%nintp
                write(*,*)

                ! Reads number of element's nodes
                keyword_name = "nen"
                call readIntegerKeywordValue(keyword_name,mesh_%nen, 0_4)
                write(iout,*) "Number of element's nodes:",mesh_%nen
                write(iout,*)
                write(*,*) "Number of element's nodes:",mesh_%nen
                write(*,*)

                ! Reads number of element's nodes
                keyword_name = "linflag"
                call readIntegerKeywordValue(keyword_name,scalar_%linflag, 1)
                write(iout,*) "(1) Linear analysis (0) Nonlinear:",scalar_%linflag
                write(iout,*)
                write(*,*) "(1) Linear analysis (0) Nonlinear:",scalar_%linflag
                write(*,*)

                ! Reads number of materials
                keyword_name = "numat"
                call readIntegerKeywordValue(keyword_name,mesh_%numat, 0_4)
                write(iout,*) "Number of materials:",mesh_%numat
                write(*,*) "Number of materials:",mesh_%numat
                !write(iout,*)

                ! Reads material components properties
                keyword_name = "mat"
                call readRealMatrixValues(keyword_name,scalar_%mat, 0.0d0)
                do i=1,mesh_%numat
                write(iout,*) "Material:", i, (scalar_%mat(i,j), j=1,3)
                write(*,*) "Material:", i, (scalar_%mat(i,j), j=1,3)
                enddo
                write(iout,*)
                write(*,*)

                ! Reads boundary conditions constant data
                keyword_name = "vbc_flow"
                call readBoundaryConditions(keyword_name,scalar_%kbc, &
                    scalar_%vbc, 0.0d0)
                keyword_name = "vbc_tracer"
                call readBoundaryConditions(keyword_name,scalar2_%kbc, &
                    scalar2_%vbc, 0.0d0)

                ! Reads transient/stationary mode
                keyword_name = "transient"
                call readIntegerKeywordValue(keyword_name,scalar2_%transient, 0)
                write(iout,*) "Transient problem (0=no; 1=yes):",scalar2_%transient
                write(*,*) "Transient problem (0=no; 1=yes):",scalar2_%transient
                scalar_%transient = 0
                !write(iout,*) "Transient problem (0=no; 1=yes):",scalar_%transient
                !write(iout,*)

                if (scalar2_%transient .eq. 1) then
                ! Reads time step value
                keyword_name = "dt"
                call readRealKeywordValue(keyword_name,scalar2_%dt, 0.0d0)
                write(iout,*) "Time step:",scalar2_%dt
                write(*,*) "Time step:",scalar2_%dt
                !write(iout,*) "Time step:",scalar_%dt
                !write(iout,*)

                ! Reads number of time steps
                keyword_name = "nsteps"
                call readIntegerKeywordValue(keyword_name,scalar2_%nsteps, 1)
                write(iout,*) "Number of time's steps:",scalar2_%nsteps
                write(iout,*) "Total time:",scalar2_%nsteps*scalar2_%dt
                write(*,*) "Number of time's steps:",scalar2_%nsteps
                write(*,*) "Total time:",scalar2_%nsteps*scalar2_%dt
                scalar_%nsteps = 1

                ! Reads time steps index that solution will be printed
                keyword_name = "tprint"
                allocate(scalar_%tprint(1)); scalar_%tprint = 1
                call readIntegerArrayValues(keyword_name,scalar2_%tprint, 0)
                !write(iout,*) "Time steps which solution will be record"
                write(iout,*) "Time steps of solution that will be record"
                write(*,*) "Time steps of solution that will be record"
                do i=1,size(scalar2_%tprint)
                write(iout,*) scalar2_%tprint(i)
                write(*,*) scalar2_%tprint(i)
                enddo
                else
                    scalar2_%nsteps = 1
                    write(iout,*) "Number of time's steps:",scalar2_%nsteps
                    write(*,*) "Number of time's steps:",scalar2_%nsteps
                    allocate(scalar_%tprint(1)); scalar2_%tprint = 1
                    write(*,*) "Time steps of solution that will be record"
                    do i=1,size(scalar2_%tprint)
                    write(iout,*) scalar2_%tprint(i)
                    write(*,*) scalar2_%tprint(i)
                    enddo
                endif
                write(iout,*);
                write(*,*);

1111            format(1x,2(a))
2222            format(1x,(a))

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

            !> Realizes preprocessor routines for Tracer problem.
            !! @param mesh_     A mesh structure
            !! @param scalar_   A scalar structure
            !! @param scalar_   Another scalar structure
            !! @author Diego T. Volpatto
            subroutine preprocessorTracer(mesh_, scalar_, scalar2_)

                use meshStructure
                use scalarStructure
                use mIO

                implicit none

                type(mesh) :: mesh_
                type(scalarStructureSystem) :: scalar_, scalar2_

                call setupPhaseTracer(mesh_,scalar_,scalar2_)
                call read_nodes(mesh_)
                call read_elems(mesh_)
                call mallocGlobalKF(scalar_, mesh_%nnodes)
                call mallocGlobalKF(scalar2_, mesh_%nnodes)
                call mallocElemKF(scalar_, mesh_%nen)
                call mallocElemKF(scalar2_, mesh_%nen)

            endsubroutine
        endmodule
