        !> A FINite ELement program for general purpose problems. 
        !! The present is based in the book "Finite Elements: An
        !! Introduction" wrote by Eric Becker, Graham Carey and Tinsley
        !! Oden.
        !!
        !! Due to the evolution of Fortran programming language, the
        !! code developed here incorporate several changes comparing to
        !! the original given in the book cited before. Modular paradigm
        !! was employed, as well a little of derived data structure.
        !!
        !! Implementations by Diego T. Volpatto.
        !! email: volpatto@lncc.br or dtvolpatto@gmail.com
        !! @author Diego Tavares Volpatto
        program finel

          use mUtilities
          use msetup
          use meshStruture
          use mshapeFunctions,  only: setint, xi, w
          use mscalar, only: localElem 

          implicit none

          real*8, allocatable :: xeps(:)
          real*8 :: pi
          integer :: nelem
          type(mesh) :: malha

          !nelem = 20

          !pi = 4.d0*datan(1.d0)
          !call linspace(-1.0d0, 1.0d0, nelem, xeps)
          !call test_shpf1d(4,nelem,xeps)
          !call setint
          !call localElem(0.d0,1.d0,2,4,1)
          !call quad1(4, 0.d0, 5.d0)
          call read_nodes(malha)

        end program
