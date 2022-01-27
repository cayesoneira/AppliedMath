program gauss
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa
real(kind=clreal),allocatable,dimension(:)::u,b,bb,r,aux
real(kind=clreal)::deter,piv,factor,suma
integer::k,n,i,j

!TOMA DE DATOS
print*,'GAUSS USUAL'
print*,'Introducir dimensión de la matriz:'
read*,n
print*,n
allocate(a(n,n),aa(n,n),u(n),b(n),bb(n),r(n),aux(n))
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

aa=a
bb=b

!MÉTODO
deter=1.
do k=1,n-1
	piv=a(k,k)
	!comprobacion de que el
	!k-esimo pivote no es nulo
	if(abs(piv)<1.e-12) then
		print*,'Pivote nulo en la etapa: ',k
		stop
	end if
	!actualizacion del determinante
	deter=deter*piv
	!eliminacion
	do i=k+1,n
		factor=a(i,k)/piv
		do j=k+1,n
			a(i,j)=a(i,j)-factor*a(k,j)
		end do
		b(i)=b(i)-factor*b(k)
	end do
end do
do i=1,n
	do j=1,n
		if(i>j)then
			a(i,j)=a(i,j)-a(i,j)
		end if
	end do
end do
!comprobacion de que el
!ultimo pivote no es nulo
if(abs(a(n,n))<1.e-12) then
	print*,'Pivote nulo en la etapa: ',n
	stop
end if
!fin del calculo del determinante
deter=deter*a(n,n)
print*,'Determinante de A:',deter
print*,'Matriz U triangular superior:'
do i=1,n
	print*,a(i,1:n)
end do

!RESOLUCIÓN DEL SISTEMA
do i=n,1,-1
	suma=0.
	do j=i+1,n
		suma=suma+a(i,j)*u(j)
	end do
	u(i)=(b(i)-suma)/a(i,i)
end do
print*,'Solución:'
print*,u

!RESIDUO
a=aa
b=bb
r=0.
aux=0.
do j=1,n
	b(1:n)=b(1:n)-a(1:n,j)*u(j)
end do
r=-b
print*,'Residuo:'
print*,r
print*,'Norma del residuo:'
print*,sqrt(dot_product(r,r))

!FINAL
deallocate(a,aa,u,b,bb,r,aux)
end program gauss
