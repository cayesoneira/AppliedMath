program cholesky
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa,at
real(kind=clreal),allocatable,dimension(:)::u,b,bb,r,auxxx,v
integer,allocatable,dimension(:)::ip
real(kind=clreal)::deter,aux,auxx,suma
integer::n,i,j,k

!TOMA DE DATOS
print*,'CHOLESKY'
print*,'Introducir dimensión de la matriz:'
read*,n
allocate(a(n,n),aa(n,n),u(n),b(n),bb(n),ip(n),r(n),auxxx(n),at(n,n),v(n))
print*,'Introducir la matriz a por filas:'
do i=1,n
	read*,a(i,1:n)
end do
do i=1,n
	print*,a(i,1:n)
end do
print*,'Introducir el vector b:'
read*,b
print*,b
ip=(/(i,i=1,n)/)
print*,'Matriz A:'
do i=1,n
	print*,a(i,:)
end do
print*,'Vector b:'
print*,b
aa=a
bb=b

!MÉTODO
deter=1.
do i=1,n
	do j=1,i-1
		a(i,j)=a(i,j)-sum(a(i,1:j-1)*a(j,1:j-1))
		a(i,j)=a(i,j)/a(j,j)
	end do
	a(i,i)=a(i,i)-sum(a(i,1:i-1)*a(i,1:i-1))
	if(a(i,i)<1.e-12)then
		print*,'Non ten descomposición.'
		stop
	end if
	a(i,i)=sqrt(a(i,i))
	deter=deter*a(i,i)
end do
deter=deter**2

print*,'Determinante de A:'
print*,deter

do i=1,n
	do j=1,n
		if(i<j)then
			a(i,j)=a(i,j)-a(i,j)
		end if
	end do
end do

print*,'Matriz B:'
do i=1,n
	print*,a(i,:)
end do

!RESOLUCIÓN DEL SISTEMA
b(1)=b(1)/a(1,1)
do i=2,n
	suma=0.
	do j=1,i-1
		suma=suma+a(i,j)*b(j)
	end do
	b(i)=(b(i)-suma)/a(i,i)
end do
print*,'Primera solución:'
print*,b
a=transpose(a)
print*,'Traspuesta de B:'
do i=1,n
	print*,a(i,:)
end do
b(n)=b(n)/a(n,n)
do i=n-1,1,-1
	suma=0.
	do j=i+1,n
		suma=suma+a(i,j)*b(j)
	end do
	b(i)=(b(i)-suma)/a(i,i)
end do
print*,'Solución final:'
print*,b
u=b
!RESIDUO
a=aa
b=bb
r=0.
auxxx=0.
do j=1,n
	auxxx=auxxx+a(:,j)*u(j)
end do
r=auxxx-b
print*,'Residuo:'
print*,r
print*,'Norma del residuo:'
print*,sqrt(dot_product(r,r))

!FINAL
deallocate(a,aa,u,b,bb,ip,r,v,auxxx)
end program cholesky
