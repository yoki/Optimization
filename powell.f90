 !. Update of Delta
dmult = 9.0d-1*fmin + 1.0d-1*fnp - fsq

! revise Delta
if(dmult.ge.0) then 
	! increase Delta
	sp = 0.0d0
	ss = 0.0d0
	do i = 1,n
		sp = sp + dabs(f(i)* (f(i)-w(nw+i)))
		ss = ss + (f(i)-w(nw+i))**2
	end do
	pj = 1. + dmult/ (sp+dsqrt(sp*sp+dmult*ss))
	spns = 4.0d0
	sp = dmin1(spns,tinc,pj)
	tinc = pj/sp
	dd = dmin1(dm,sp*dd)
else
	! decrease Delta
	dd = dmax1(dss,2.5d-1*dd)
	tinc = 1.0d0 
end if

! Interchange x and x+delta if F(x+delta) is less than F(x)
if( (dmult.ge. 0) .or. (fsq-fmin.lt.0)) then
	fmin = fsq
	do  i = 1,n
		sp = x(i)
		x(i) = w(nx+i)
		w(nx+i) = sp
		sp = f(i)
		f(i) = w(nf+i)
		w(nf+i) = sp
		w(nw+i) = -w(nw+i)
	end do
!!  	if(is-1.gt.0) go to 50 ** consider this line!
end if

! Calculate the change in F and X
do i = 1,n
    x(i) = x(i) - w(nx+i)
    f(i) = f(i) - w(nf+i)
continue


!. Keep linear independence

! calculate the change in x and its angle with the first direction (150)
dn = 0.0d0
sp = 0.0d0
do  i = 1,n
	f(i) = dmult*x(i) + anmult*f(i)
	dn = dn + f(i)*f(i)
	sp = sp + f(i)*w(nd+i)
continue
ds = 2.5d-1*dn
! test if extra step is needed for independence
if(w(ndc+1)-dtest.ge.0 .and. sp*sp-ds.lt.0) then
! take the extra step and update the direction matrix. 
	do  i = 1,n
		x(i) = w(nx+i) + dstep*w(nd+i)
		w(ndc+i) = w(ndc+i+1) + 1.0d0
	continue
	w(nd) = 1.0d0
	do  i = 1,n
		k = nd + i
		sp = w(k)
		if (n.lt.2) cycle
		do  j = 2,n
			w(k) = w(k+n)
			k = k + n
		continue
		w(k) = sp
	continue
	go to 1
else
	! express the new direction in terms of those of the 
	! direction matrix, and update the counts in w(ndc+1) etc
	sp = 0.0d0
	k = nd
	do 64 i = 1,n
		x(i) = dw
		dw = 0.0d0
		do  j = 1,n
			k = k + 1
			dw = dw + f(j)*w(k)
		continue
		if (is \=2) then
		    w(ndc+i) = w(ndc+i) + 1.0d0
		    sp = sp + dw*dw
		    if(sp-ds.le.0) cycle
		    is = 1
		    kk = i
		    x(1) = dw
		else
		    x(i) = dw
		end if
		w(ndc+i) = w(ndc+i+1) + 1.0d0
	continue
end if
w(nd) = 1.0d0

! reorder the direction so that kk is first
if(kk-1.ge.0) then
	ks = ndc + kk*n
	do  i = 1,n
	k = ks + i
	sp = w(k)
	do  j = 2,kk
	  w(k) = w(k-n)
	  k = k - n
	continue
	w(k) = sp
	continue
end if

   ! calculation of new direction
   70 do  i = 1,n
        w(nw+i) = 0.0d0
      continue
      sp = x(1)*x(1)
      k = nd
      if (n.lt.2) go to 99
      do 75 i = 2,n
        ds = dsqrt(sp* (sp+x(i)*x(i)))
        dw = sp/ds
        ds = x(i)/ds
        sp = sp + x(i)*x(i)
        do 76 j = 1,n
          k = k + 1
          w(nw+j) = w(nw+j) + x(i-1)*w(k)
          w(k) = dw*w(k+n) - ds*w(nw+j)
        continue
      continue
   99 sp = 1.0d0/dsqrt(dn)
      do 77 i = 1,n
        k = k + 1
        w(k) = sp*f(i)
      continue
! Calculate the next vector x and predict right hand side
   80 fnp = 0.0d0
      k = 0
      do 78 i = 1,n
        x(i) = w(nx+i) + f(i)
        w(nw+i) = w(nf+i)
        do 79 j = 1,n
          k = k + 1
          w(nw+i) = w(nw+i) + w(k)*f(j)
   79   continue
        fnp = fnp + w(nw+i)**2
   78 continue
      go to 1