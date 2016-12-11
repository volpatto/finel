        !> Module for shape functions computations and relate
        !! operations.
        !! @author Diego T. Volpatto
        module mshapeFunctions
            
            implicit none
            
            
            real*8  ::  xi(20,20), &  !< Gauss point integration 1D
                        w(20,20), &   !< Gauss weights 1D
                        xiq(2,9), & !< Gauss points to rectangular element
                        wq(9),& !< Gauss weights to rectangular element
                        xit(2,4), & !< Gauss points to triangle element
                        wt(4) !< Gauss weights to triangle element

            integer :: quadext

            contains
            
            !> Gauss quadrature data set routine - 1D.
            subroutine setint

            !... This routine defines the values of the parameters
            !... required for the numerical integration of element
            !... matrices and vectors.

                implicit none

                integer :: i

                ! Clear
                w = 0.d0; xi = 0.d0

                if (quadext .eq. 0) then
                ! Gaussian quadrature order 1

                xi(1,1) = 0.0
                w(1,1)  = 2.0

                ! Gaussian quadrature order 2

                xi(1,2) = -1./dsqrt(3.0d0)
                xi(2,2) = -xi(1,2)
                w(1,2)  = 1.0
                w(2,2)  = w(1,2)

                ! Gaussian quadrature order 3

                xi(1,3) = -dsqrt(3.0d0/5.0d0)
                xi(2,3) = 0.0
                xi(3,3) = -xi(1,3)
                w(1,3)  = 5.0/9.0
                w(2,3)  = 8.0/9.0
                w(3,3)  = w(1,3)

                ! Gaussian quadrature order 4

                xi(1,4) = -0.8611363116
                xi(2,4) = -0.3399810436
                xi(3,4) = -xi(2,4)
                xi(4,4) = -xi(1,4)
                w(1,4)  = 0.3478548451
                w(2,4)  = 0.6521451549
                w(3,4)  = w(2,4)
                w(4,4)  = w(1,4)

                ! Gaussian quadrature order 5

                w(1,5) = .236926885056189
                w(2,5) = .478628670499366
                w(3,5) = .568888888888888
                w(4,5) = .478628670499366
                w(5,5) = .236926885056189
                xi(1,5)=-.906179845938664
                xi(2,5)=-.538469310105683
                xi(3,5)= 0.d0
                xi(4,5)= .538469310105683
                xi(5,5)= .906179845938664

                ! Gaussian quadrature order 6

                 w(1,6) = .171324492397170
                 w(2,6) = .360761573048139
                 w(3,6) = .467913934572691
                 w(4,6) = .467913934572691
                 w(5,6) = .360761573048139
                 w(6,6) = .171324492397170
                 xi(1,6)=-.932469514203152
                 xi(2,6)=-.661209386466265
                 xi(3,6)=-.238619186083197
                 xi(4,6)= .238619186083197
                 xi(5,6)= .661209386466365
                 xi(6,6)= .932469514203152

                 ! Gaussian quadrature order 7

                 w(1,7) = .129484966168870
                 w(2,7) = .279705391489277 
                 w(3,7) = .381830050505119
                 w(4,7) = .417959183673469
                 w(5,7) = .381830050505119
                 w(6,7) = .279705391489277
                 w(7,7) = .129484966168870
                 xi(1,7)=-.949107912342759
                 xi(2,7)=-.741531185599394
                 xi(3,7)=-.405845151377397
                 xi(4,7)= .0d0
                 xi(5,7)= .405845151377397
                 xi(6,7)= .741531185599394
                 xi(7,7)= .949107912342759

                 ! Gaussian quadrature order 8

                 w(1,8) = .101228536290376
                 w(2,8) = .222381034453374
                 w(3,8) = .313706645877887
                 w(4,8) = .362683783378362
                 w(5,8) = .362683783378362
                 w(6,8) = .313706645877887
                 w(7,8) = .222381034453374
                 w(8,8) = .101228536290376
                 xi(1,8)=-.960289856497536
                 xi(2,8)=-.796666477413627
                 xi(3,8)=-.525532409916329
                 xi(4,8)=-.183434642495650
                 xi(5,8)= .183434642495650
                 xi(6,8)= .525532409916329
                 xi(7,8)= .796666477413627
                 xi(8,8)= .960289856497536

                else if (quadext .eq. 1) then

                ! Extend Gaussian quadrature log singular (Diego)
                xi(1,20) = 0.058684713389643455348
                xi(2,20) = 0.20085590974338672323
                xi(3,20) = 0.34102346665699661535
                xi(4,20) = 0.47499755685218773898
                xi(5,20) = 0.59933105712959326204
                xi(6,20) = 0.71098630175917040439
                xi(7,20) = 0.80729488074114865580
                xi(8,20) = 0.88597995933835106670
                xi(9,20) = 0.94519165897723331249
                xi(10,20)= 0.98353845345472765660
                w(1,20) = 0.40077930960551466949
                w(2,20) = 0.22788836401795107375
                w(3,20) = 0.14812360023344707567
                w(4,20) = 0.096547816551015655829
                w(5,20) = 0.060649749122039659929
                w(6,20) = 0.035609402612321070364
                w(7,20) = 0.018803993942658554996
                w(8,20) = 0.0083797607340355934749
                w(9,20) = 0.0027604277995848179492
                w(10,20)= 0.00045757538143182854107
                do i=11,20
                xi(i,20) = -xi(i-10,20)
                w(i,20) = w(i-10,20)
                enddo

                ! Universal Gauss quadrature left-singular order 3
                xi(1,6) = 0.074330911998429636136
                xi(2,6) = 0.49067720343094899076
                xi(3,6) = 0.93020808634535177637
                w(1,6) = 0.23536348707782220590
                w(2,6) = 0.59474998941244362493
                w(3,6) = 0.74068285030463078840
                xi(4,6) = -xi(1,6)
                xi(5,6) = -xi(2,6)
                xi(6,6) = -xi(3,6)
                w(4,6) = w(1,6)
                w(5,6) = w(2,6)
                w(6,6) = w(3,6)
                
                ! Universal Gauss quadrature left-singular order 4
                xi(1,8) = 0.04043445929275625351
                xi(2,8) = 0.30329580016011919663
                xi(3,8) = 0.68411187359652082596
                xi(4,8) = 0.96067483342440115801
                w(1,8) = 0.13725763004144154931
                w(2,8) = 0.37180716249909688774
                w(3,8) = 0.50325894862261033029
                w(4,8) = 0.55847258563174785190
                do i=5,8
                xi(i,8) = -xi(i-4,8)
                w(i,8) = w(i-4,8)
                enddo

                else
                    print*, "Unknown quadrature kind"; stop
                endif

            end subroutine

      !> Gauss quadrature data set routine - 2D.
      SUBROUTINE SETINT2
!**************************************************************
!
!     SET INTEGRATING CONSTANTS
!
!**************************************************************
      implicit none

!.....QUADRILATERAL ELEMENTS
!.....GAUSS INTEGRATION ORDER 3*3
      XIQ(1,1)=-DSQRT(3.d0/5.d0)
      XIQ(2,1)=XIQ(1,1)
      XIQ(1,2)=0.d0
      XIQ(2,2)=XIQ(1,1)
      XIQ(1,3)=-XIQ(1,1)
      XIQ(2,3)=XIQ(1,1)
      XIQ(1,4)=XIQ(1,1)
      XIQ(2,4)=0.d0
      XIQ(1,5)=0.d0
      XIQ(2,5)=0.d0
      XIQ(1,6)=-XIQ(1,1)
      XIQ(2,6)=0.d0
      XIQ(1,7)=XIQ(1,1)
      XIQ(2,7)=-XIQ(1,1)
      XIQ(1,8)=0.d0
      XIQ(2,8)=-XIQ(1,1)
      XIQ(1,9)=-XIQ(1,1)
      XIQ(2,9)=-XIQ(1,1)
      WQ(1)=25.d0/81.d0
      WQ(2)=40.d0/81.d0
      WQ(3)=WQ(1)
      WQ(4)=WQ(2)
      WQ(5)=64.d0/81.d0
      WQ(6)=WQ(2)
      WQ(7)=WQ(1)
      WQ(8)=WQ(2)
      WQ(9)=WQ(1)
!.....TRIANGULAR ELEMENTS
!.....GAUSS INTEGRATION ORDER 3
      XIT(1,1)=1.d0/3.d0
      XIT(2,1)=XIT(1,1)
      XIT(1,2)=2.d0/15.d0
      XIT(2,2)=11.d0/15.d0
      XIT(1,3)=XIT(1,2)
      XIT(2,3)=XIT(1,2)
      XIT(1,4)=XIT(2,2)
      XIT(2,4)=XIT(1,2)
      ! My weights. (Diego)
      WT(1)=-27.d0/96.d0
      WT(2)=25.d0/96.d0
      WT(3)=WT(2)
      WT(4)=WT(2)

      ENDSUBROUTINE

!****************************************************************************************
            
            !> Calculates the values of the shape functions and their
            !! derivatives - 1D.
            !! @param xl        [in] specified value of master element coord
            !! @param n         [in] number of element nodes
            !! @param psi       [out] shape function values
            !! @param dpsi      [out] derivatives shape functions values
            !! @author Diego Volpatto
            subroutine shpf1d(xl,n,psi,dpsi)

            !... Calculates the values of the shape functions psi and
            !... their derivatives dpsi with respect to the master
            !... element coordinate at a specified value xl.
                implicit none

                integer  ::  n
                real*8   ::  xl, psi(n), dpsi(n)

                ! Testing shape function order
                if (n .lt. 2 .or. n .gt. 4) then
                    write(*,*) "Error in call to shape", n
                    stop
                endif

                ! Clear
                psi = 0.d0; dpsi = 0.d0

                ! Linear shape functions
                if (n .eq. 2) then
                    psi(1) = 0.5*(1.0-xl)
                    psi(2) = 0.5*(1.0+xl)
                    dpsi(1)= -0.5
                    dpsi(2)= 0.5
                    return
                endif

                ! Quadratic shape functions
                if (n .eq. 3) then
                    psi(1) = xl*(xl-1.0)*0.5
                    psi(2) = 1.0-xl**2.0
                    psi(3) = xl*(xl+1.0)*0.5
                    dpsi(1)= xl-0.5
                    dpsi(2)= -2.0*xl
                    dpsi(3)= xl+0.5
                    return
                endif

                ! Cubic shape functions
                if (n .eq. 4) then
                    psi(1) = 9.0/16.0*(1.0/9.0-xl**2.0)*(xl-1.0)
                    psi(2) = 27.0/16.0*(1.0-xl**2.0)*(1.0/3.0-xl)
                    psi(3) = 27.0/16.0*(1.0-xl**2.0)*(1.0/3.0+xl)
                    psi(4) = -9.0/16.0*(1.0/9.0-xl**2.0)*(1.0+xl)
                    dpsi(1)= -9.0/16.0*(3.0*xl**2.0-2.0*xl-1.0/9.0)
                    dpsi(2)= 27.0/16.0*(3.0*xl**2.0-2.0/3.0*xl-1.0)
                    dpsi(3)= 27.0/16.0*(-3.0*xl**2.0-2.0/3.0*xl+1.0)
                    dpsi(4)= -9.0/16.0*(-3.0*xl**2.0-2.0*xl+1.0/9.0)
                    return
                endif

            end subroutine

            !> Calculates the values of the shape functions and their
            !! derivatives - 2D.
            !! @param xl        [in] specified value of master element coord
            !! @param n         [in] number of element nodes
            !! @param psi       [out] shape function values
            !! @param dpsi      [out] derivatives shape functions values
            !! @author Diego Volpatto
      SUBROUTINE shpf2d(XL,N,PSI,DPSI)
!**************************************************************
!
!     SHAPE FUNCTIONS;
!
!*************************************************************
      implicit none
      real*8 :: PSI(9),DPSI(2,9)
      real*8 :: XL(2)
      integer :: n
      real*8 :: x,y,xx,yy,xy,y2,x2,xy2,xxy,xyy

      x  = xl(1); y = xl(2)
      xx = x*x; yy = y*y; xy = x*y
      xxy= xx*y; xyy = x*yy
      y2 = 2.0d0*y; x2 = 2.0d0*x; xy2 = xy*2.0d0

      if (n .eq. 9) then
!.....SHAPE FUNCTIONS FOR
!.....QUADRILATERAL 9-NODE ELEMENTS
          PSI(1)=0.25d0*(xx-x)*(yy-y)
          PSI(2)=0.25d0*(xx+x)*(yy-y)
          PSI(3)=0.25d0*(xx+x)*(yy+y)
          PSI(4)=0.25d0*(xx-x)*(yy+y)
          PSI(5)=0.5d0*(1.d0-xx)*(yy-y)
          PSI(6)=0.5d0*(xx+x)*(1.d0-yy)
          PSI(7)=0.5d0*(1.-xx)*(yy+y)
          PSI(8)=0.5d0*(xx-x)*(1.d0-yy)
          PSI(9)=(1.d0-xx)*(1.d0-yy)
          DPSI(1,1)=(0.5d0*x-0.25d0)*(yy-y)
          DPSI(2,1)=(0.5d0*y-0.25d0)*(xx-x)
          DPSI(1,2)=(0.5d0*x+0.25d0)*(yy-y)
          DPSI(2,2)=(0.5d0*y-0.25)*(xx+x)
          DPSI(1,3)=(0.5d0*x+0.25d0)*(yy+y)
          DPSI(2,3)=(0.5d0*y+0.25d0)*(xx+x)
          DPSI(1,4)=(0.5d0*x-0.25d0)*(yy+y)
          DPSI(2,4)=(0.5d0*y+0.25d0)*(xx-x)
          DPSI(1,5)=-x*(yy-y)
          DPSI(2,5)=(y-0.5d0)*(1.d0-xx)
          DPSI(1,6)=(x+0.5d0)*(1.d0-yy)
          DPSI(2,6)=-y*(xx+x)
          DPSI(1,7)=-x*(yy+y)
          DPSI(2,7)=(y+0.5d0)*(1.d0-xx)
          DPSI(1,8)=(x-0.5d0)*(1.d0-yy)
          DPSI(2,8)=-y*(xx-x)
          DPSI(1,9)=-x2*(1.d0-yy)
          DPSI(2,9)=-y2*(1.d0-xx)
      else if (n .eq. 8) then
!.....SHAPE FUNCTIONS FOR
!.....QUADRILATERAL 8-NODE ELEMENTS
          PSI(1)=0.25d0*(-1.d0+xy+xx+yy-xxy-xyy)
          PSI(2)=0.5d0*(1.d0-y-xx+xyy)
          PSI(3)=0.25d0*(-1.d0-xy+xx+yy-xxy+xyy)
          PSI(4)=0.5d0*(1.d0+x-yy-xyy)
          PSI(5)=0.25d0*(-1.d0+xy+xx+yy+xxy+xyy)
          PSI(6)=0.5d0*(1.d0+y-xx-xxy)
          PSI(7)=0.25d0*(-1.d0-xy+xx+yy+xxy-xyy)
          PSI(8)=0.5d0*(1.d0-x-yy+xyy)
          DPSI(1,1)=0.25d0*(y+x2-xy2-yy)
          DPSI(2,1)=0.25d0*(x+y2-xy2-xx)
          DPSI(1,2)=-x+xy
          DPSI(2,2)=0.5d0*(-1.d0+xx)
          DPSI(1,3)=0.25d0*(-y+x2-xy2+yy)
          DPSI(2,3)=0.25d0*(-x+y2+xy2-xx)
          DPSI(1,4)=0.5d0*(1.d0-yy)
          DPSI(2,4)=-y-xy
          DPSI(1,5)=0.25d0*(y+x2+xy2+yy)
          DPSI(2,5)=0.25d0*(x+y2+xy2+xx)
          DPSI(1,6)=-x-xy
          DPSI(2,6)=0.5d0*(1.d0-xx)
          DPSI(1,7)=0.25d0*(-y+x2+xy2-yy)
          DPSI(2,7)=0.25d0*(-x+y2-xy2+xx)
          DPSI(1,8)=0.5d0*(-1.d0+yy)
          DPSI(2,8)=-y+xy
      else if (n .eq. 4) then
!.....SHAPE FUNCTIONS FOR
!.....QUADRILATERAL 4-NODE ELEMENTS
          PSI(1)=0.25d0*(1.d0-x-y+xy)
          PSI(2)=0.25d0*(1.d0+x-y-xy)
          PSI(3)=0.25d0*(1.d0+x+y+xy)
          PSI(4)=0.25d0*(1.d0-x+y-xy)
          DPSI(1,1)=0.25d0*(-1.d0+y)
          DPSI(2,1)=0.25d0*(-1.d0+x)
          DPSI(1,2)=0.25d0*(1.d0-y)
          DPSI(2,2)=0.25d0*(-1.d0-x)
          DPSI(1,3)=0.25d0*(1.d0+y)
          DPSI(2,3)=0.25d0*(1.d0+x)
          DPSI(1,4)=0.25d0*(-1.d0-y)
          DPSI(2,4)=0.25d0*(1.d0-x)
      else if (n .eq. 3) then
!.....SHAPE FUNCTIONS FOR
!.....TRIANGULAR 3-NODES ELEMENTS
          PSI(1)=(1.d0-x-y)
          PSI(2)=x
          PSI(3)=y
          DPSI(1,1)=-1.d0
          DPSI(2,1)=-1.d0
          DPSI(1,2)=1.d0
          DPSI(2,2)=0.d0
          DPSI(1,3)=0.d0
          DPSI(2,3)=1.d0
      else if (n .eq. 6) then
!.....SHAPE FUNCTIONS FOR
!.....TRIANGULAR 6-NODES ELEMENTS
          PSI(1)=2.d0*(1.d0-x-y)*(0.5d0-x-y)
          PSI(2)=2.d0*x*(x-0.5d0)
          PSI(3)=2.d0*y*(y-0.5d0)
          PSI(4)=4.d0*(1.d0-x-y)*x
          PSI(5)=4.d0*xy
          PSI(6)=4.d0*y*(1.d0-x-y)
          DPSI(1,1)=-3.d0+4.d0*(x+y)
          DPSI(2,1)=DPSI(1,1)
          DPSI(1,2)=4.d0*x-1.d0
          DPSI(2,2)=0.0d0
          DPSI(1,3)=0.d0
          DPSI(2,3)=4.d0*y-1.
          DPSI(1,4)=4.d0-8.d0*x-4.d0*y
          DPSI(2,4)=-4.d0*x
          DPSI(1,5)=4.d0*y
          DPSI(2,5)=4.d0*x
          DPSI(1,6)=-4.d0*y
          DPSI(2,6)=4.d0-4.d0*x-8.d0*y
      ELSE
          WRITE(*,100) N; stop
      endif

100   FORMAT(' ERROR IN CALL TO SHAPE,N='I3)
      ENDSUBROUTINE
!
!
!

        end module
