program jacobi
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a
real(kind=clreal),allocatable,dimension(:)::u,b,r,uold,aux

!DECLARACIONES PARTICULARES
real(kind=clreal)::error,eps
integer::nitmax,i,n,iter,j

!TOMA DE DATOS
print*,'Introducir orden de la matriz:'
read*,n
print*,n
allocate(a(n,n),u(n),b(n),r(n),uold(n),aux(n))
print*,'Introducir la matriz A:'
do i=1,n
	read*,a(i,1:n)
	print*,a(i,1:n)
end do

print*,'Introducir el vector b:'
read*,b
print*,b

!MÉTODO
print*,'Fijar iteraciones máximas:'
read*,nitmax
print*,nitmax
print*,'Fijar épsilon de parada:'
read*,eps
print*,eps

do i=1,n
	if(abs(a(i,i))<1.e-12)then
		print*,'Elemento diagonal',i,'nulo.'
		stop
	end if
end do
do iter=1,nitmax
	uold=u
	do i=1,n
		u(i)=(b(i)-sum(a(i,1:i-1)*uold(1:i-1))-sum(a(i,i+1:n)*uold(i+1:n)))/a(i,i)
	end do
	error=sum(abs(u-uold))
	print*,'Iterante:',iter,u
	if(error<eps)then
		print*,'Se satisface el test de parada en la iteración',iter
		print*,'El error entre iterantes en norma 1 en la parada es:',error
		exit
	end if
end do

if(iter>nitmax)then
	print*,'Se efectuaron',nitmax,'iteraciones sin satisfacer el test de parada.'
end if

print*,'La solución aproximada es'
print*,u

!RESIDUO
r=0.
aux=0.
do j=1,n
	aux=aux+a(:,j)*u(j)
end do
r=aux-b
print*,'Residuo:'
print*,r
print*,'Norma del residuo:'
print*,sqrt(dot_product(r,r))

!FINAL
deallocate(a,u,b,r,aux)
end program jacobi
