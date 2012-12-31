module MyUtility 
! pack of utilities
! last modified : Feb 10, 2008
! Written by Yoki Okawa


integer, parameter :: db = selected_real_kind(p=13,r=200)


real(kind=db), parameter  :: two = 2.0_db
real(kind=db), parameter  :: one = 1.0_db
real(kind=db), parameter  :: p5  = 0.50_db
real(kind=db), parameter  :: zero = 0.0_db
real(kind=db), parameter  :: p1  = 0.10_db
real(kind=db), parameter  :: p01  = 0.01_db
real(kind=db), parameter  :: p0001 = 0.0001_db
real(kind=db), parameter  :: p00001 = 0.00001_db

contains
	subroutine MyError(str)
	! display error message and stop 
	implicit none
	CHARACTER(LEN=*), INTENT(IN) :: str
	write(*,*) 'error ', str
	write(*,*) 'press any key to continue'
	read(*,*) 
	
	return
	end subroutine MyError
	
	subroutine MatrixWrite(X,str)
	! write matrix as a matrix
	! Input:
	! X: input of matrix
	! 
	! Optional Input:
	! str: (1f10.3) type string to specify the format of output.
	
	implicit none
	real(kind=db), intent(IN),dimension(:,:) :: X
	character(len=*), intent(IN), optional :: str
	integer :: i,j,n,m
	character(len=100) :: fmt
	
	! set default value of str
	if (.not. present(str)) then
		fmt  = '(1G12.6)'
	else
		fmt = str
	end if
	
	n = size(x,1) 
	m = size(x,2)
	
	do i = 1, n
		do j = 1, m
			write(*,fmt,advance ='no') X(i,j)
		end do
		write(*,*) ' '
	end do
	
	return
	end subroutine MatrixWrite
	
	subroutine VectorWrite(vec,str)
	! write a vector as a horizontal vector
	! Input:
	! X: input of matrix
	! 
	! Optional Input:
	! str: (1f10.3) type string to specify the format of output.
	
	implicit none
	real(kind=db), intent(IN),dimension(:) :: vec
	character(len=*), intent(IN), optional :: str
	integer :: i,j,n,m
	character(len=100) :: fmt
	
	! set default value of str
	if (.not. present(str)) then
		fmt  = "(' ' , 1F20.10)"
	else
		fmt = str
	end if
	
	n = size(vec) 
	
	do i = 1, n
		write(*,fmt) vec(i)
	end do

	
	return
	end subroutine VectorWrite
	
	
	pure function outer(x,y)
	! caclulate the outer product 
	! outer = x* y^t if x is column vector
	implicit none
	real(kind=db), intent(IN),dimension(:) :: x, y
	real(kind=db), dimension(size(x),size(y)) :: outer
	outer = spread(x,dim=2,ncopies=size(y)) * &
		spread(y,dim=1,ncopies=size(x))
	end function outer
	
	pure function norm(x)
	! calculate L^2 norm of X
	implicit none
	real(kind=db),intent(IN), dimension(:) :: x
	real(kind=db) ::  norm
	
	norm = ( dot_product(x,x) ) ** (0.50_db)
	end function norm
	
	pure function eye(M)
	! returns identity matrix of dimension M
	implicit none
	integer,intent(IN) :: m
	real(kind=db), dimension(M,M) ::  eye
	integer :: i, j
	
	forall (i=1:M, j =1:M)
		eye(i,j) = 0
	end forall
	
	forall (i = 1:M)
		eye(i,i) = 1
	end forall
	
	end function eye
	
	pure function diag(M)
	! returns a vector which contains diagonal element of M
	implicit none
	real(kind=db),intent(IN), dimension(:,:) :: M
	real(kind=db), dimension(min(size(M,1),size(M,2))) ::  diag
	integer :: n, i
	
	n= min(size(M,1),size(M,2))
	
	forall (i=1:n)
		diag(i) = M(i,i)
	end forall
	end function diag
	


end module MyUtility
