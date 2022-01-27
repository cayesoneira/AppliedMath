subroutine iter(f,eps,nitmax,x0)
implicit none
real,intent(inout)::eps,nitmax,x0
real::x1
real::f
integer::i

i=0
do
	x1=f(x0)

	if(abs(x1-x0)<eps*(1+abs(x0))then
		print*,'El punto fijo es',x1,'y su residuo es',f(x1)
		stop
	end if

	if(i>nitmax)then
		print*,'Alcanzado el número máximo de iteraciones.'
		print*,'El  punto fijo tal y como estamos sería',x1
		stop
	end if

	x0=x1
	i=i+1
	print*,'Iteración',i,'x=',x1
end do

end subroutine iter
