        !> A FINite ELement program for general purpose problems. 
        !! The present code is based in the book "Finite Elements: An
        !! Introduction" wrote by Eric Becker, Graham Carey and Tinsley
        !! Oden.
        !!
        !! Due to the evolution of Fortran programming language, the
        !! code developed here incorporates several changes comparing to
        !! the original given in the book cited before. Modular paradigm
        !! was employed, as well a little of derived data structure.
        !!
        !! Implementations by Diego T. Volpatto.
        !! email: volpatto@lncc.br or dtvolpatto@gmail.com
        !! @author Diego Tavares Volpatto
        program finel

          use mUtilities
          use msetup
          use meshStructure
          use scalarStructure
          use mprocessor
          use mIO
          use mscalar 

          implicit none

          type(mesh) :: malha
          type(scalarStructureSystem) :: potencial
          real*8 :: t1, t2

          call cpu_time(t1)
          call openFiles
          call preprocessor(malha, potencial)
          call processor(malha, potencial)
          call cpu_time(t2)
          write(iout,*) "Execution elapsed time in minutes:", (t2-t1)/60.d0
          write(*,*) "Execution elapsed time in minutes:", (t2-t1)/60.d0
          call closeFiles

        end program
