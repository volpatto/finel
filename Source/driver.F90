        program finel

          use mUtilities
          use mshapeFunctions,  only: setint, xi, w

          implicit none

          real*8, allocatable :: xeps(:)
          real*8 :: pi
          integer :: nelem

          nelem = 20

          pi = 4.d0*datan(1.d0)
          call linspace(-1.0d0, 1.0d0, nelem, xeps)
          call test_shpf1d(4,nelem,xeps)
          !call setint
          !call quad1(4, 0.d0, 5.d0)

        end program
