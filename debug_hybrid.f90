program mydebug
use hybrid
! call test1()
! call qrtest()

! call finitedifftest
! call doglegtest
call FsolveHybridTest
write(*,*) 'press any key to continue'
read(*,*) 


end program mydebug





subroutine doglegtest
use myutility; use hybrid
! subroutine to test dogleg

implicit none
real(kind=db), dimension(3,3) :: jacob, Q, R
real(kind=db), dimension(3) :: x0, p, Qtf
external  :: funs2
real(kind=db), dimension(3) ::  fval
real(kind=db) :: delta
integer :: flag
! real(kind=db), dimension(2,2) :: GetJacobian

x0(1) = 0.5_db
x0(2) = 1.0_db
x0(3) = 1.5_db
delta = 0.10_db

call funs2(x0,fval)
write(*,*) 'funs2([0.50, 1.00. 1.50 ]):'
write(*,*) fval

x0 = x0+1
call funs2(x0,fval)
call GetJacobian(jacob,  funs2, x0, 0.0001_db,fval)
call QRfactorization(jacob,Q,R)

Qtf = matmul(transpose(Q), fval)
call dogleg(p,Q,R,delta,Qtf,flag)

write(*,*) 'p:'
write(*,*) p
write(*,*) 'flag:'
write(*,*) flag

end subroutine doglegtest



subroutine finitedifftest()
use myutility; use hybrid, only : GetJacobian
! test finite difference
implicit none

real(kind=db), dimension(2,2) :: jacob
real(kind=db), dimension(2) :: x0
external  :: funs1
real(kind=db), dimension(2) ::  fval
! real(kind=db), dimension(2,2) :: GetJacobian

x0(1) = 2.0_db
x0(2) = 3.0_db

call funs1(x0, fval)
call GetJacobian(jacob,  funs1, x0, 0.001_db,fval)

call MatrixWrite(jacob)

end subroutine finitedifftest

subroutine funs1(x, fval0)
use myutility; 
implicit none
real(kind=db), intent(IN), dimension(:) :: x
real(kind=db), intent(OUT), dimension(:) :: fval0

fval0(1) = x(2) - x(1)**2
fval0(2) = 2 - x(1) - x(2)

end subroutine funs1



subroutine funs3(x, fval0)
use myutility; 
implicit none
real(kind=db), intent(IN), dimension(:) :: x
real(kind=db), intent(OUT), dimension(:) :: fval0

fval0(1) = 10_db*(x(2) - x(1)**2)
fval0(2) = 2 - x(1) - x(2)

end subroutine funs3


subroutine funs2(x, fval0)
use myutility; 
! a little difficult function
! solution x = [0.50, 1.00. 1.50 ]  (+ 2pi*n)

implicit none
real(kind=db), intent(IN), dimension(:) :: x
real(kind=db), intent(OUT), dimension(:) :: fval0

fval0(1) = 1.20_db * sin(x(1)) -1.40_db*cos(x(2))+ 0.70_db*sin(x(3)) &  
			- 0.517133908732486_db

fval0(2) = 0.80_db * cos(x(1)) -0.50_db*sin(x(2))+ 1.00_db*cos(x(3)) &
			- 0.352067758776053_db

fval0(3) = 3.50_db * sin(x(1)) -4.25_db*cos(x(2))+ 2.80_db*cos(x(3))  &
			+ 0.4202312501553165_db

end subroutine funs2


subroutine qrtest
use myutility; use hybrid, only : QRfactorization, QRupdate
! test QR factorization and update

implicit none
real(kind=db), dimension(5,5) ::  Q,  A3
real(kind=db), dimension(5,5) :: A, R ,A2, A4, A5, A6, B1, B2
real(kind=db), dimension(5) :: u,v
INTEGER              :: isize
INTEGER,ALLOCATABLE  :: iseed(:)

! set a seed. 
CALL RANDOM_SEED(SIZE=isize)
ALLOCATE( iseed(isize) )
CALL RANDOM_SEED(GET=iseed)
iseed = 1
CALL RANDOM_SEED(PUT = iseed)  

CALL RANDOM_NUMBER(A)           ! generate  random number
A = A - 0.5
write(*,*) 'A:'
call MatrixWrite( A )

Q = 1
R = 1
call QRfactorization(A,Q,R)

A2 = matmul(Q , R)
A3 = matmul(Q , transpose(Q))

write(*,*) 'Q:'
call MatrixWrite( Q )
write(*,*) 'R:'
call MatrixWrite( R )
write(*,*) 'Q*R:'
call MatrixWrite( A2)
write(*,*) "Q*Q':"
call MatrixWrite( A3)
A3 = A2 - A
write(*,*) "Q*R - A:"
call MatrixWrite( A3)

! update test
call RANDOM_NUMBER(u)
call RANDOM_NUMBER(v)
u = (u-0.5) * 1
v = (v-0.5) * 1

write(*,*) ' '
write(*,*) 'update A '
call MatrixWrite( A+outer(u,v))
call QRupdate(Q,R,u,v)


write(*,*) 'update Q:'
call MatrixWrite( Q )
write(*,*) 'update R:'
call MatrixWrite( R )
write(*,*) "Q*Q':"
call MatrixWrite( matmul(Q , transpose(Q)))
write(*,*) "Q*R - A:"
call MatrixWrite( matmul(Q , R) - (A + outer(u,v)))


! write(*,*) "Q*Q':"
! call MatrixWrite( A3)
! 

end subroutine qrtest


subroutine temp
! temp place to put code


end subroutine temp

