        !> Contains subroutine to compute numerical solution of linear
        !! systems Ax=b. 
        !!
        !! The present module has a general purpose such
        !! that the routines here intend to be independent of others
        !! module.
        !!
        !! @author Diego T. Volpatto
        module msolver

            implicit none

            contains
            
            !> Applies Gauss reduction in A(n,n) to obtain a superior
            !! triangular equivalent form.
            !! @param A     [in/out]A matrix A(n,n)
            !! @param n     [in]Number of rows/columns of matrix A
            !! @author Diego T. Volpatto
            subroutine tri(A, n)

                implicit none

                integer :: n
                real*8 :: A(n,n)

                integer :: n1, j1, i, j, k
                real*8 :: fac
                real*8, parameter :: tiny_val = 1.d-30

                n1 = n-1
                ! Eliminate degree of freedom i
                do i=1,n1

                ! Check for excessively small pivot
                if (dabs(A(i,i)).lt.tiny_val) then
                    print*, i, A(i,i)
                    print*, "Reduction failed due to small pivot"
                    stop
                endif

                j1 = i+1
                ! Modify row j
                do j=j1,n
                if (A(j,i).eq.0.0d0) cycle
                fac = A(j,i)/A(i,i)
                do k=j1,n
                A(j,k)=A(j,k)-A(i,k)*fac
                enddo
                enddo
                enddo

            endsubroutine

            !> Does the forward substitution on the right-side-hand.
            !! @param A     [in]A matrix A(n,n)
            !! @param x     [out]Solution vector
            !! @param b     [in/out]RHS-vector
            !! @param n     [in]Number of solution points
            !! @author Diego T. Volpatto
            subroutine rhsub(A, x, b, n)

                implicit none

                real*8 :: A(n,n), x(n), b(n)
                integer :: n
                integer :: n1, j1, i, j, ib

                n1 = n-1
                ! Begin forward reduction of right hand side
                do i=1,n1
                j1=i+1
                do j=j1,n
                b(j)=b(j)-b(i)*A(j,i)/A(i,i)
                enddo
                enddo
                ! Begin back substitution
                x(n)=b(n)/A(n,n)
                do i=1,n1
                ib=n-i
                j1=ib+1
                do j=j1,n
                b(ib)=b(ib)-A(ib,j)*x(j)
                enddo
                x(ib)=b(ib)/A(ib,ib)
                enddo

            endsubroutine

            subroutine solverGaussSeidel(A, x, b, x0, n)

                implicit none
                
                real*8 :: A(n,n), x(n), b(n), x0(n)
                integer :: n
		integer :: k
		integer, parameter :: kmax = 1.d2
                real*8 :: sum_ax

		x = x0
		do k=1,kmax
		!do i=1,n
			!sum_ax = 0.d0
			!do j=1,n
			     !if (j .ne. i) sum_ax = A(i,j)*x(j)
			!enddo
			!x(i) = (b(i)-sum_ax)/A(i,i)
		!enddo
                if (k .eq. kmax) stop "kmax";
		enddo
                
            endsubroutine
        endmodule
