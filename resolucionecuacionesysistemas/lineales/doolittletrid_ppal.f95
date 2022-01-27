program doolittletrid
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa
real(kind=clreal),allocatable,dimension(:)::u,b,bb,r,ad,al,au,v,auau,alal,adad
integer,allocatable,dimension(:)::ip,aux
real(kind=clreal)::deter
integer::n

!DECLARACIONES PARTICULARES
integer::i

!TOMA DE DATOS
print*,'DOOLITTLE TRIDIAGONAL'
print*,'Introducir dimensión de la matriz:'
read*,n
print*,n
allocate(ad(n),adad(n),au(n),auau(n),al(n),alal(n),v(n),u(n),b(n),bb(n),r(n),aux(n))
print*,'Introducir la diagonal superior de A rematándola con un 0:'
read*,au
print*,'Introducir la diagonal principal de A:'
read*,ad
print*,'Introducir la diagonal inferior de A rematándola con un 0:'
read*,al
print*,'Introducir el vector b'
read*,b

adad=ad
alal=al
auau=au
bb=b

!MÉTODO
deter=ad(1)
do i=1,n-1
	if(abs(ad(i))<1.e-12)then
		print*,'No existe la factorización LU.'
		stop
	end if
	al(i)=al(i)/ad(i)
	ad(i+1)=ad(i+1)-al(i)*au(i)
	deter=deter*ad(i+1)
end do
if(abs(ad(n))<1.e-12)then
	print*,'No existe la factorización LU.'
	stop
end if

print*,'Superdiagonal final:'
print*,au
print*,'Diagonal principal final:'
print*,ad
print*,'Subdiagonal final:'
print*,al

print*,'Determinante:'
print*,deter
!RESOLUCIÓN DEL SISTEMA
v(1)=b(1)
do i=2,n
	v(i)=b(i)-al(i-1)*v(i-1)
end do
print*,'Primera solución, vector v:'
print*,v

u(n)=v(n)/ad(n)
do i=n-1,1,-1
	u(i)=(v(i)-au(i)*u(i+1))/ad(i)
end do
print*,'Solución final, vector u:'
print*,u

!RESIDUO
b=bb
ad=adad
au=auau
al=alal
r=0.
aux=0.
r(1)=ad(1)*u(1)+au(1)*u(2)-b(1)
do i=2,n
	r(i)=al(i-1)*u(i-1)+ad(i)*u(i)+au(i)*u(i+1)-b(i)
end do
r(n)=al(n-1)*u(n-1)+ad(n)*u(n)-b(n)
print*,'Residuo:'
print*,r
print*,'Norma del residuo:'
print*,sqrt(dot_product(r,r))

!FINAL
deallocate(u,v,au,al,ad,auau,alal,adad,b,bb,r,aux)
end program doolittletrid
