program choleskytrid
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa
real(kind=clreal),allocatable,dimension(:)::u,b,bb,r,ad,au,v,auau,adad,x,y
integer,allocatable,dimension(:)::ip
real(kind=clreal)::deter,aux,auxx
integer::n,i,j,k

!TOMA DE DATOS
print*,'CHOLESKY TRIDIAGONAL'
print*,'Introducir dimensión de la matriz:'
read*,n
print*,n
allocate(ad(n),adad(n),au(n),auau(n),v(n),u(n),b(n),bb(n),r(n),x(n),y(n))
print*,'Introducir la diagonal superior de A rematándola con un 0:'
read*,au
print*,'Introducir la diagonal principal de A:'
read*,ad
print*,'Introducir el vector b:'
read*,b

print*,'Diagonal superior de A:'
print*,au
print*,'Diagonal principal de A:'
print*,ad
print*,'Vector b:'
print*,b

adad=ad
auau=au
bb=b

!MÉTODO
do i=1,n
	if(ad(i)<0)then
		print*,'No es definida positiva.'
		stop
	end if
end do
aux=0.
auxx=0.
x(1)=sqrt(ad(1))
do j=1,n-1
	y(j)=au(j)/x(j)
	x(j+1)=sqrt(ad(j+1)-y(j)**2)
end do
au=y
ad=x

print*,'Diagonal principal de B:'
print*,ad
print*,'Diagonal inferior de B:'
print*,au

deter=1.
do i=1,n
	deter=deter*ad(i)
end do

print*,'Determinante de A:'
print*,deter**2

!RESOLUCIÓN DEL SISTEMA
v(1)=b(1)/ad(1)
do i=2,n
	v(i)=(b(i)-au(i-1)*v(i-1))/ad(i)
end do
print*,'Primera solución:'
print*,v

u(n)=v(n)/ad(n)
do i=n-1,1,-1
	u(i)=(v(i)-au(i)*u(i+1))/ad(i)
end do
print*,'Solución final:'
print*,u

!RESIDUO
b=bb
ad=adad
au=auau

r=0.
r(1)=ad(1)*u(1)+au(1)*u(2)-b(1)
do i=2,n-1
	r(i)=au(i-1)*u(i-1)+ad(i)*u(i)+au(i)*u(i+1)-b(i)
end do
r(n)=au(n-1)*u(n-1)+ad(n)*u(n)-b(n)
print*,'Residuo:'
print*,r
print*,'Norma del residuo:'
print*,sqrt(dot_product(r,r))

!FINAL
deallocate(u,v,au,ad,auau,adad,b,bb,r)
end program choleskytrid
