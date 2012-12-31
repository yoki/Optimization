module hybrid
use MyUtility
! This module is for solving system of nonlinear equations using modified Powell's 
! Hybrid method. Original code is from in MINPACK.

! subroutine list
! FsolveHybrid    : main subroutine to perform the modified Powell's Hybrid
!                   method. It calls UpdateDelta, dogleg, QRfactorization, QRupdate
!                   and GetJacobian.
! UpdateDelta     : It updates Delta, the size of the trust region
! dogleg          : It finds the next step p within the trust region
! QRfactorization : It performs QR factorization
! QRupdate        : It updates QR factorization after the Broyden's update
! applyGivens     : It apply Givens transformation. It is called by QRupdate.
! 
! FsolveHybridTest: Sample subroutine to illustrate the usage of FsolveHybrid
!                   It contains funstest1 and funstest2. 
!                   Just call FsolveHybrid from the main program. 

! Last Updated : Feb 17, 2009


contains
	subroutine FsolveHybrid(fun, x0, xout, xtol, info, fvalout, JacobianOut, &
		JacobianStep, display, MaxFunCall, factor, NoUpdate,deltaSpeed)
	! Solve system of nonlinear equations using modified Powell's Hybrid method. 
	implicit none
		INTERFACE
			SUBroutine  fun (x, fval0)
			use myutility; 
			implicit none
			real(kind=db), intent(IN), dimension(:)  :: x    
			real(kind=db), intent(OUT), dimension(:) :: fval0
			end subroutine fun
		END INTERFACE

	!. Declearation of variables
	real(kind=db), intent(IN), dimension(:)  :: x0         ! Initial value
	real(kind=db), intent(OUT),dimension(size(x0)) :: xout ! solution 
	real(kind=db), intent(IN),  optional     :: xtol       ! Relative tol of X
	integer      , intent(OUT), optional     :: info       ! Information on output


	! output function value and Jacobian 
	real(kind=db), intent(OUT), dimension(size(x0)),          optional :: fvalOut
	real(kind=db), intent(OUT), dimension(size(x0),size(x0)), optional :: JacobianOut


	real(kind=db), intent(IN), optional :: JacobianStep ! Relative step for Derivative
	real(kind=db), intent(IN), optional :: factor       ! inital value of delta       
	integer      , intent(IN), optional :: display      ! Controls display on iteration
	integer      , intent(IN), optional :: MaxFunCall   ! Max number of Function call
	integer      , intent(IN), optional :: NoUpdate     ! Jacobian Recalculation info
	real(kind=db), intent(IN), optional :: deltaSpeed   ! How fast Delta should decrease
	
	integer  :: n                  ! number of variables
	integer  :: IterationCount     ! number of iterations
	integer  :: FunctionCount      ! number of function calls
	integer  :: GoodJacobian       ! number of concective sucessfull iterations
	integer  :: BadJacobian        ! number of concective failing iterations
	integer  :: SlowJacobian       ! Degree of slowness after repeated Jacobian Update
	integer  :: SlowFunction       ! number of concective failure of improving 
	integer  :: info2              ! information for the termination
	integer  :: i                  ! index for the loop
	integer  :: display2           ! varible to control displaying information
	integer  :: MaxFunCall2        ! Maximum number of Function call
	integer  :: NoUpdate2          ! variable to control Jacobian Update 
	integer  :: DirectionFlag      ! Flag for the output of dogleg
	integer  :: UpdateJacobian     ! Mark for the updating Jacobian 

	real(kind=db) :: temp              
	real(kind=db) :: Delta             ! size of trust region
	real(kind=db) :: pnorm             ! norm of step
	real(kind=db) :: ActualReduction   ! 1 - norm(fvalold)/norm(fvalnew)
	real(kind=db) :: ReductionRatio    ! ActuanReduction / PredictedReduction
	real(kind=db) :: xtol2             ! torelance of x
	real(kind=db) :: JacobianStep2     ! Finite difference step size
	real(kind=db) :: factor2           ! initial value of delta
	real(kind=db) :: DeltaOld          ! Used for display purpose
	real(kind=db) :: deltaSpeed2       ! How fast Delta should decrease
	

	real(kind=db), dimension(size(x0))  :: xbest         ! best x so far
	real(kind=db), dimension(size(x0))  :: xold          ! x befor update
	real(kind=db), dimension(size(x0))  :: xnew          ! xold + p
	real(kind=db), dimension(size(x0))  :: fvalbest      ! fun(xbest)
	real(kind=db), dimension(size(x0))  :: fvalpredicted ! fun(xold)+ J*p
	real(kind=db), dimension(size(x0))  :: fvalold       ! fun(xold) 
	real(kind=db), dimension(size(x0))  :: fvalnew       ! fun(xold+p)
	real(kind=db), dimension(size(x0))  :: p             ! predicted direction
	real(kind=db), dimension(size(x0))  :: Psidiag       ! Normaization coefs
	real(kind=db), dimension(size(x0))  :: Qtfval        ! Q^T * fvalbest

	real(kind=db), dimension(size(x0),size(x0)) :: Q,R    ! results of QR factorization
	real(kind=db), dimension(size(x0),size(x0)) :: J      ! Jacobian
	real(kind=db), dimension(size(x0),size(x0)) :: PsiInv ! Inverse of nomarization mat
	
	CHARACTER(LEN=79) :: st1, st2 ! used for output
	CHARACTER(LEN=6)  :: st6      ! used for output
	CHARACTER(LEN=12) :: st12     ! used for output

	
	!. Initialize values
	! counters
	IterationCount = 0
	FunctionCount  = 0
	GoodJacobian   = 0
	BadJacobian    = 0
	SlowFunction   = 0
	SlowJacobian   = 0
	info2          = 0
	n              = size(x0)

	! Set default values for optional inputs
	xtol2         = p0001
	display2      = 0
	MaxFunCall2   = n*100
	NoUpdate2     = 0
	factor2       = 100
	deltaSpeed2   = 0.25_db      

	if (present(xtol))         xtol2        = xtol
	if (present(display))      display2     = display
	if (present(MaxFunCall))   MaxFunCall2  = MaxFunCall
	if (present(NoUpdate))     NoUpdate2    = NoUpdate
	if (present(factor))       factor2      = factor
	if (present(deltaSpeed))   deltaSpeed2  = deltaSpeed
	

	JacobianStep2 = xtol * p1
	if (present(JacobianStep)) JacobianStep2  = JacobianStep

	! Jacobian and function values
	xbest = x0
	call fun(xbest, fvalbest) ! output: fvalbest
	FunctionCount = FunctionCount + 1
	call GetJacobian(J, fun, xbest, JacobianStep2,fvalbest) ! output : J 
	FunctionCount = FunctionCount + n
	call QRfactorization(J,Q,R) ! output: Q, R
	Qtfval = matmul(transpose(Q),fvalbest)

	! calculate normalzation matrix Psi
	PsiInv  = 0
	do i = 1, n
		! normalization factor is R(i,i) unless R(i,i) = 0
		temp = 1
		if(R(i,i) /= zero) temp = abs(R(i,i))
		PsiInv(i,i) = 1/temp
		PsiDiag(i) = temp
	end do

	! calculate initial value of Delta
	Delta = factor2 * norm(PsiDiag * xbest)
	if(Delta == 0 ) Delta = 1

	! check initial guess is good or not
	if (norm(fvalbest) == zero) info2 = 1
	
	! display first line
	if(display2 ==1) then
		write(*,*) 'FsolveHybrid:'
	    st1= "       Norm      Actual  Trust-Region     Step  Jacobian   Direction"
     	st2= " iter  f(x)    Reduction    Size          Size  Recalculate Type"
		write(*,*) st1
		write(*,*) st2
		st1 = "('   0 ',1G11.5)" 
		write(*,st1) norm(fvalbest)
	end if 


	! ****************************
	!         main loop
	! ****************************
	do
		IterationCount = IterationCount + 1
		
		! old values are values at the start of the iteration
		fvalold = fvalbest
		xold    = xbest
		
		!. *** calculate the best direction *** 
		call dogleg(p,Q,matmul(R,PsiInv),Delta,Qtfval,DirectionFlag)
		! output: p, DirectionFlag
		p = matmul(PsiInv, p)
		
		!. update the trust region
		call fun(xold + p, fvalnew)
		FunctionCount = FunctionCount +1
		fvalpredicted = fvalbest + matmul(Q,matmul(R,p))
		DeltaOld = Delta
		call UpdateDelta(Delta,GoodJacobian,BadJacobian,&
			ActualReduction, ReductionRatio,fvalold,fvalnew,&
			fvalpredicted, PsiDiag*p,deltaSpeed2 )
		! output: Delta, GoodJacobian, BadJacobian, 
		!         ActualReduction, ReductionRatio
		
		! get the best value so far
		if(norm(fvalnew) < norm(fvalold) .and. ReductionRatio > p0001) then
			xbest = xold +p
			fvalbest = fvalnew
		end if

		
		!. *** Check convergence ***
		! Sucessful Convergence
		if(Delta < xtol*norm(PsiDiag*xbest) .or. norm(fvalbest) == 0) info2 = 1
		
		! Too much function call
		if(FunctionCount > MaxFunCall2) info = 2
		
		! tol is too small 
		if(Delta < 100 * epsilon(Delta) * norm(PsiDiag * xbest)) info2 = 3
		
		! Not successful based on Jacobian
		if(ActualReduction > p1) SlowJacobian = 0
		if(SlowJacobian == 5) info2 = 4
		! if Jacobian is recalculated every time, we do not performe this test
		if(noupdate2  == 1) SlowJacobian = 0
		
		! Not sucessful based on Function Value
		SlowFunction = SlowFunction + 1
		if(ActualReduction > p01) SlowFunction = 0
		if( SlowFunction == 10) info2 = 5
		
		
		!.***  Update Jacobian *** 
		pnorm = norm(p)
		UpdateJacobian = 0
		if(BadJacobian == 2 .or. pnorm == 0 .or. noupdate == 1) then
			! calculate Jacobian using finite difference
			call GetJacobian(J, fun, xbest, JacobianStep2,fvalbest) ! output : J 
			FunctionCount = FunctionCount + n
			call QRfactorization(J,Q,R) ! output: Q, R
			Qtfval = matmul(transpose(Q),fvalbest)

			! recalculate normalzation matrix Psi
			do i = 1, n
				! normalization factor is R(i,i) unless R(i,i) = 0
				temp = 1
				if(R(i,i) /= zero) temp = abs(R(i,i))
				PsiInv(i,i) = min(PsiInv(i,i), 1/temp)
				PsiDiag(i) = 1/PsiInv(i,i) 
			end do
			
			! take care of counts
			BadJacobian = 0
			SlowJacobian = SlowJacobian +1
			UpdateJacobian = 1
		else if (ReductionRatio > p0001) then
			! Broyden's Rank 1 update
			call QRupdate(Q,R,fvalnew - fvalpredicted,p/((pnorm)**2))
			Qtfval = matmul(transpose(Q),fvalbest)
		end if
		
		! display iteration
		if(display2 ==1) then
			st1 = "(1I4,' ',1G11.5,' ',1G11.5,' ', 1G11.5,' ', 1G11.5,' ',  2A)" 
			st6 = "      "
			if(UpdateJacobian == 1) st6 = ' Yes  '
			select case (DirectionFlag)
				case (1)
					st12 = 'Newton'
				case (2)
					st12 = 'Cauchy'
				case (3)
					st12 = 'Combination'
			end select
			write(*,st1) IterationCount, norm(fvalbest), 100.0_db * ActualReduction, &
				 DeltaOld, norm(PsiDiag *p),st6, st12
		end if 

		
		! exit check 
		if( info2 /= 0) exit 
	end do

	!. prepare output
	xout = xbest
	fvalOut = fvalbest
	JacobianOut = J
	
	! display result
	if(display2 ==1) then
		select case (info2)
			case (0)
				write(*,*) 'Bad input'
			case (1)
				if(norm(fvalbest)>p1) then
					write(*,*) ' Trust Region shrinks enough so no progress is possible.'
					write(*,*) ' Make sure function value is close to zero enough'
				else
					write(*,*) ' Sucessful convergence'
				end if
			case (2)
				write(*,*) ' Too much Function call'
			case (3)
				write(*,*) ' Tol too small'
			case (4)
				write(*,*) ' Too much Jacobian '
			case (5)
				write(*,*) ' Slow Objective function Improvement'
		end select
		if (norm(fvalbest)<xtol .and. info2 /= 1) then
 		   write(*,*) ' Note: Function value is fairly small.'
 		   write(*,*) '       although it is not required torelance.'
 		   write(*,*) '       It might be converged.'
		end if       
	end if 
	
	end subroutine FsolveHybrid

	pure subroutine UpdateDelta(Delta,GoodJacobian,BadJacobian, &
	 	ActualReduction, ReductionRatio, oldfval,newfval,predictedfval,p, deltaSpeed2)
	! update Trust region Delta
	
	implicit none
	real(kind=db), intent(INOUT) 	:: Delta			! Trust region size
	integer, intent(INOUT) 			:: GoodJacobian		! number of sucessful iteration
	integer, intent(INOUT) 			:: BadJacobian		! number of bad iteration
	real(kind=db), intent(IN), dimension(:) :: oldfval	! f(xold)
	real(kind=db), intent(IN), dimension(:) :: newfval	! f(xold+p)
	real(kind=db), intent(IN), dimension(:) :: predictedfval	! f(xold) + J^T*p
	real(kind=db), intent(IN), dimension(:) :: p		! suggested direction
	real(kind=db), intent(IN)  :: deltaSpeed2    		! How fast delta decreases
	real(kind=db), intent(OUT) :: ActualReduction		! reduction for actual fval
	real(kind=db), intent(OUT) :: ReductionRatio		! ActRed/PredictedRed
	
	real(kind=db) :: pnorm					! norm(p)
	
	real(kind=db) :: PredictedReduction		! scaled reduction for predicted fval

	! prepare comparison of values
	PredictedReduction = 1 - (norm(predictedfval) / norm(oldfval))**2
	ActualReduction    = 1 - (norm(newfval)       / norm(oldfval))**2
	pnorm = norm(p)
	
	! calculate the ratio of actual to predicted reduction
	ReductionRatio = 0 ! special case if PredictedReduction = 0
	if (PredictedReduction /= 0) &
		ReductionRatio = ActualReduction/PredictedReduction	  
	
	if (ReductionRatio < p1 ) then 
		! prediction was not good.  shrink trust region
		Delta = Delta * deltaSpeed2
		BadJacobian = BadJacobian +1
		GoodJacobian = 0
	else
		BadJacobian = 0
		GoodJacobian = 1+GoodJacobian
		if (GoodJacobian>1 .or. ReductionRatio > p5 ) then 
			! prediction was fair. expand trust egion
			Delta = max(Delta, two * pnorm) 
			
		else if (abs(1 -ReductionRatio) < p1 ) then
			! prediction was very good. (the ratio is close to one).
			! Expand trust region 
			Delta = two * pnorm
		end if
	end if	
	
	end subroutine updateDelta
	


	pure subroutine dogleg(p,Q,R,delta,Qtf,flag)
	! find linear combination of newton direction and steepest descent 
	! direction. 
	
	! flag : it indicates the type of the p (optional)
	! flag = 1  : newton direction
	! flag = 2  : Steepest descent direction
	! flag = 3  : Linear combination of bot
	
	implicit none
	real(kind=db), intent(out), dimension(:) :: p	! output direction
	real(kind=db), intent(in),  &
		dimension(size(p),size(p)) 	::  Q, R		! QR decomposition of Jacobian
	real(kind=db), intent(in)		:: delta        ! trust region parameter
	real(kind=db), intent(in),  &
				dimension(size(p))  :: Qtf			! Q^T * f(x)
	integer, intent(out), optional 	:: flag 		! flag about the type of p

	integer 							:: i, tempflag
	integer 							:: n		! number of variables = size(p)
	real(kind=db) 						:: gnorm, mu, nunorm, theta, mugnorm, temp
	real(kind=db) 						:: Jgnorm
	real(kind=db), dimension(size(p)) 	:: nu		! Newton direction
	real(kind=db), dimension(size(p)) 	:: g		! Steepest descent direction
	real(kind=db), dimension(size(p)) 	:: mug
	
	n = size(p)
	
	! ******************************
	! calculate newton direction
	! ******************************
	! prepare a small value in case diagonal element of R is zero
	temp = epsilon(mu) * maxval(abs(diag(R)))
	
	nu(n) = -1* Qtf(n) /temp       ! this is special value in case
	if (R(n,n) /= zero) nu(n) =  -1* Qtf(n) / R(n,n)  ! normal case
	
	! solve backwards 
	do i = n-1, 1, -1
		if (R(i,i)==0) then 
			! special value
			nu(i) = (-1*Qtf(i) - dot_product(R(i,i+1:n),nu(i+1:n))) / temp
		else
			! normal value
			nu(i) = (-1*Qtf(i) - dot_product(R(i,i+1:n),nu(i+1:n))) / R(i,i) 
		end if
	end do
	nunorm = norm(nu)

	
	if (nunorm < delta) then 
		! newton direction 
		p = nu
		tempflag = 1
	else
		! newton direction was not accepted. 
		g = - one * matmul(transpose(R),Qtf)   ! Steepest descent
		gnorm = norm(g)
		Jgnorm = norm(matmul(Q,matmul(R,g)))
		
		if (Jgnorm == 0) then 
			! special attention if steepest direction is zero
			p = delta * nu/nunorm
			flag = 3
			
		else if ((gnorm**2) *gnorm / (Jgnorm**2) > delta) then
			! accept steepest descent direction
			p = delta *g /gnorm
			tempflag = 2
		else
			! linear combination of both
			! calculate the weight of each direction
			mu = gnorm**2  / Jgnorm**2
			mug = mu *g
			mugnorm = norm(mug)
			theta = (delta**2 - mugnorm**2) / (dot_product(mug, nu-mug) +  &
				((dot_product(nu,mug)-delta**2)**2 + (nunorm**2-delta**2)  &
					* (delta**two - mugnorm**2))**p5)
			
			p = (1-theta) * mu * g + theta*nu
			tempflag = 3
			
		end if
	end if

	if (present(flag)) flag = tempflag
	end subroutine dogleg
	
	subroutine QRfactorization(A,Q,R)
	! Calculate QR factorizaton using Householder transformation. 
	! You can obtain better speed and stability by using LAPACK routine.	
	
	! It finds ortogonal matrix Q and upper triangular R such that
	!
	!               A  = Q * [R; ZeroMatrix]
	!
	! Arguments for this subroutine
	! A: m by n (m>=n) input matrix for the QR factorization to be computed 
	! Q: m by m output orthogonal matrix 
	! R: m by n output upper triangular matrix. 
	! 
	! Written by Yoki Okawa
	! Date: Feb 10, 08


	implicit none
	
	real(kind=db), INTENT(IN), dimension(:,:)    :: A
	real(kind=db), INTENT(INOUT), dimension(:,:) :: Q,R
	real(kind=db), allocatable,  dimension(:,:)  :: X2

	real(kind=db), dimension(size(A,1))     :: u
	real(kind=db), dimension(size(A,1),size(A,1))     :: P

	integer :: m, n, mQ1,mQ2, nR, mR,i,j
	
	m = size(A,1)		! number of rows in A
	n = size(A,2)		! number of columns in A
	
	! check size of the outputs
	mQ1 = size(Q,1)     ! number of rows in Q
	mQ2 = size(Q,2)		! number of columns in Q
	mR = size(R,1)      ! number of rows in R
	nR = size(R,2)		! number of columns in R
	if (n /= nR .or. m /= mQ1 .or. m /= mQ2  .or. m/=mR ) then 
	    call myerror &
	    ('QRfactorization : output matrix dimensions do not match with inputs')
	end if
	
	if (m<n) then 
		call myerror &
	('QRfactorization : number of rows must be equal or greater than number of columns')
	end if
	
	
	! main loop
	R = A
	Q = eye(M)
	do i = 1, n
		! checke if all elements are already zero
		u(i:m) = R(i:M,i)
		u(i) = u(i) - norm(R(i:M,i))
		if(norm(u) == 0 ) then
			continue
		end if
		
		! no need for nth column if m == n
		if( (m==n).and. (i==n) )then 
			exit 
		end if		
		
		P = eye(M) ! P is identity at the left top part
		
		! right bottom (m-i+1) by (m-i+1) part of P contains numbers
		P(i:m,i:m) = eye(M-i+1) - &
			two*outer(u(i:m),u(i:m))/( norm(u(i:m))**two)
		
		R  = matmul(P , R) ! eliminate column i
		Q = matmul(Q , P)  ! update Q 
		

	end do

	end subroutine QRfactorization
	
	pure subroutine QRupdate(Q,R,u,v)
	! update QR factorization when QR -> QR +u*v
	! Q: (inout) n by n Orthogonal matrix
	! R: (inout) n by n upper triangular matrix
	! u,v: n dimensional vector
	
	implicit none
	real(kind=db), intent(INOUT), dimension(:,:) :: Q,R
	real(kind=db), intent(IN), dimension(:)      :: u, v
 
 	integer :: N                        ! dimension of Q or R
 	real(kind=db), allocatable, dimension(:)   :: w  ! Qt * w
 	real(kind=db), allocatable, dimension(:,:) ::Qt  ! transpose(Q)
 	integer :: i, j
 	real(kind=db) :: s,c, t, wnorm       ! sin, cos, tan
 	
 	N = size(u)
 	allocate(w(N))
 	allocate(Qt(N,N)) 

 	Qt = transpose(Q)
	w = matmul(Qt,u)
 	wnorm = norm(w)   ! norm of w
 	
 	! make w to unit vector
 	do i = N, 2,-1
 		! calculate cos and sin
 		if (w(i-1) == zero) then 
 			c = zero
 			s = one
 		else
 			t = w(i)/w(i-1)
 			c = one /( (t**2+one)**p5)
 			s = t * c
 		end if
 		
 		call applyGivens(R,c,s,i-1,i)
 		call applyGivens(Qt,c,s,i-1,i)
 		w(i-1) = c *w(i-1) + s* w(i)
 		w(i) = c*w(i) - s*w(i-1)
	end do
 	
 	! update R
 	R(1,:) = R(1,:) + w(1) *v
 	
 	! Transform upper Hessenberg matrix R to upper triangular matrix
 	! H in the documentation is currentry R
 	do i = 1, N-1
 		if (R(i,i) == zero) then 
 			c = zero
 			s = one
 		else
 			t = R(i+1,i)/R(i,i)
 			c = one /( (t**2+one)**p5)
 			s = t * c
 		end if
 		call applyGivens(R,c,s,i,i+1)
 		call applyGivens(Qt,c,s,i,i+1)
	end do
 	
 	Q = transpose(Qt)
	
	end subroutine QRupdate
	
	pure subroutine applyGivens(A,c2,s2,i2,j2)
	! apply givens transformation with cos, sin, index i and j to matrix A.
	! A <- P * A, P: givens transformation. 
	! P(i2,j2) = -s2; P(j2, i2) = s2; P(i2,i2) = P(j2,j2) = c2	
	
	implicit none
	real(kind=db), intent(INOUT), dimension(:,:) :: A
	real(kind=db), intent(IN) :: c2, s2
	real(kind=db), allocatable, dimension(:) :: ai, aj
	integer, intent(IN) :: i2, j2
	integer :: N

	! store original input
	N = size(A,2)
	allocate(ai(N))
	allocate(aj(N))
	
	ai = A(i2,:)
	aj = A(j2,:)
	
	! only row i and row j changes
	A(i2,:) = A(i2,:) + (c2-1) * ai
	A(i2,:) = A(i2,:) + s2 * aj
	
	! change in row j
	A(j2,:) = A(j2,:) - s2 * ai
	A(j2,:) = A(j2,:) + (c2-1) * aj

	end subroutine applyGivens
	
	subroutine GetJacobian(Jacobian, fun, x0, xrealEPS,fval)
	! Calculate Jacobian using forward difference
	! Jacobian(i,j) : derivative of fun(i) with respect to x(j)
	! 
	! n : number of dimensions of fun
	! m : number of dimensions of x0
	! fun:      function to evaluate (actually, subroutine)
	! x0 :      position to evaluate
	! xreleps : relative amount of x to change
	! fval :    value of fun(x) ( used to save time)
	
	
	implicit none
	INTERFACE
		subroutine fun(x, fval0)
		use myutility; 
		implicit none
		real(kind=db), intent(IN), dimension(:)  :: x
		real(kind=db), intent(OUT), dimension(:) :: fval0
		end subroutine fun
	END INTERFACE
	
	real(kind=db), intent(in), dimension(:)		:: x0, fval
	real(kind=db), intent(out), &
		dimension(size(fval),size(x0) ) 		:: Jacobian
	real(kind=db), intent(in)             		:: xrealEPS
	
	real(kind=db), dimension(size(fval))		:: fval0, fval1 
	real(kind=db) 								:: xdx
	real(kind=db), dimension(size(x0)) 			:: xtemp
	integer :: j , m
	
	m = size(x0)  ! number of variables
	
	! main loop (make it forall for speed)
	do  j = 1,m
		! special treatment if x0 = 0
		if(x0(j) == 0) then
			xdx = 0.001_db
		else
			xdx = x0(j) * (one + xrealEPS)
		end if

		xtemp = x0
		xtemp(j) = xdx
		
		call fun(xtemp, fval1) ! evaluate function at xtemp
		Jacobian(:,j) = (fval1 - fval) / (xdx - x0(j))
	end do
	
	end subroutine GetJacobian
	
	
	subroutine FsolveHybridTest
	! test  subroutine FsolveHybrid

	implicit none
	real(kind=db), dimension(2,2) :: jacob, Q, R
	real(kind=db), dimension(2) :: xout2, fval
	integer ::fsolveinfo

	call FsolveHybrid( &
		fun          = funstest1, &         ! Function to be solved
		x0           =(/-1.2_db,1.0_db/), & ! Initial value
		xout         =xout2, &              ! output 
		xtol         =0.00001_db, &         ! error torelance
		info         =fsolveinfo, &         ! info for the solution
		fvalout      =fval,  &              ! f(xout)
		JacobianOut  =jacob, &              ! Jacobian at x = xout
		JacobianStep =0.000001_db,  &       ! Stepsize for the Jacobian
		display      =1,  &                 ! Control for the display
		MaxFunCall   = 1000, &              ! Max number of function call
		factor       =1.0_db, &             ! Initial value of delta
		NoUpdate     = 0)                   ! control for update of Jacobian

	write(*,*) ' '
	write(*,*) 'Solution:'
	call VectorWrite(xout2)
	write(*,*) ' '
	write(*,*) 'Function Value at the solution:'
	call VectorWrite(fval)
		
	contains
		subroutine funstest1(x, fval0)
		! badly scaled function if we start from (-1.2, 1.0)
		! solution x1 = x2 = 1
		implicit none
		real(kind=db), intent(IN), dimension(:) :: x
		real(kind=db), intent(OUT), dimension(:) :: fval0

		fval0(1) = 10_db*(x(2) - x(1)**2)
		fval0(2) = 1 - x(1) 

		end subroutine funstest1


		subroutine funstest2(x, fval0)
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

		end subroutine funstest2
	end subroutine FsolveHybridTest

end module hybrid



