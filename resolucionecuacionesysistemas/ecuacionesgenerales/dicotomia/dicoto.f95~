subroutine dicoto(f,a,b,del,eps,z)
implicit none
real,intent(in)::a,b,del,eps
real::x,y,t,q,alp
real,intent(inout)::z
real::f
integer::i

x=a
y=b
i=0

do
	z=(x+y)/2
	alp=(y-x)/(2**i)

	if(f(z)==0)then
		print*,'La raíz es',z
		stop
	end if

	if(f(z)*f(y)<0)then
		t=z
		q=y
	end if

	if(f(x)*f(z)<0)then
		t=x
		q=z
	end if

	x=t
	y=q

	!Test de parada.
	if(abs(f(z))<del.or.abs(alp)<eps)then
		print*,'La aproximación de la raíz es',z

		if(abs(f(z))<del)then
			print*,'Se detuvo por acercamiento de F(x0) al 0.'
		end if

		if(abs(alp)<eps)then
			print*,'Se detuvo por acercamiento de x0 a la raíz.'
		end if

		print*,'El residuo es',f(z)
		stop
	end if

	i=i+1
end do





end subroutine dicoto
