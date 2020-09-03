 module var_rad
 double precision R
 end module var_rad


Program Circle 
	use var_rad
	implicit none
	double precision nintegral,vintegral, a , b
	double precision No, actual
	integer n
	
	
	write(*,*) 'Insert radius of Circle: '
	read(*,*) R

	a=0
	b=R
	n=1000000
	
	actual = 0.5*(2.4048)**(2)/(R**(2))
	call nsimpson(j,a,b,nintegral,n)
	no=1/(nintegral)
	
	
	call nsimpson2(H,a,b,vintegral,n)
	write (*,*) 'Ground State Approximation = ',no*vintegral
	write (*,*) 'The actual value is: ',actual
	write (*,*) 'The percent error is: ', ( ( (no*vintegral)-(actual) )/actual)*100 , ' %'
	
	





	CONTAINS

	subroutine nsimpson(j,a,b,nintegral,n)
		implicit none
		double precision j, a , b , nintegral ,s
		double precision h, x
		integer n, i

		if ((n/2)*2.ne.n) then
			n=n+1
		End IF

		s=0.0
		h=(b-a)/dble(n)
		do i=2, n-2, 2
			x=a+(i*h)
			s=s+(2.0*(j(x)) + 4.0*j(x+h) )
		End DO
		
		nintegral=(s+ j(a) + j(b) + 4.0*j(a+h))*h/3.0
	return
	End subroutine nsimpson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine nsimpson2(H,a,b,vintegral,n)
		implicit none
		double precision H, a , b , vintegral , s
		double precision hh, x
		integer n, i

		if ((n/2)*2.ne.n) then
			n=n+1
		End IF


		s=0.0
		hh=(b-a)/dble(n)
		do i=2, n-2, 2
			x=a+(i*hh)
			s=s+(2.0*H(x) + 4.0*H(x+hh))
		End DO
		a=10.0d0**(-10)
		vintegral=(s+ H(a)+ H(b)+ 4.0*H(a+hh) )*hh/3.0
	return
	End subroutine nsimpson2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	function f(x)
		use var_rad
		implicit none
		double precision f, x
		Real, Parameter :: Pi = 3.14159265359

		f = (1.0/6.0) + (1.0/3.0)*cos(1.2024*x/R) + (1.0/3.0)*cos(2.08262*x/R) + (1.0/6.0)*cos(2.4048*x/R) &
		- ( (1.0/6.0) + (1.0/3.0)*cos(1.2024) + (1.0/3.0)*cos(2.08262) + (1.0/6.0)*cos(2.4048) ) 
	return
	End 

	function j(x)
		implicit none
		double precision j, x
		
		j = x*f(x)*f(x)
	return
	End

	function g(x)
		implicit none
		double precision g, x
		g = ( f(x+ (10.0d0**(-3))) - f(x- (10.0d0**(-3))) )/(2*10.0d0**(-3)) 
	return
	End

	function s(x)
		implicit none
		double precision s,x
		s = ( f(x+ (10.0d0**(-3))) -2*f(x) + f(x- (10.0d0**(-3))) )/(10.0d0**(-6))
	return
	End


	function H(x)
		implicit none
		double precision H, x
		H = -(0.5)*f(x)*x*( s(x) + ( (1.0/x)*g(x)))
	return
	End
		

End program circle
 
		
		

	
