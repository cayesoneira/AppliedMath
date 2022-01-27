program seidel
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa
real(kind=clreal),allocatable,dimension(:)::u,b,bb,uold,auxxx,r
integer,allocatable,dimension(:)::ip

!DECLARACIONES PARTICULARES
real(kind=clreal)::error,eps,ui
integer::nitmax,i,n,iter,j

!TOMA DE DATOS
print*,'Introducir dimensión de la matriz:'
read*,n
print*,n
allocate(a(n,n),aa(n,n),u(n),b(n),r(n),uold(n),auxxx(n))
print*,'Introducir la matriz A:'
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
	error=0.
	do i=1,n
		ui=(b(i)-sum(a(i,1:i-1)*u(1:i-1))-sum(a(i,i+1:n)*u(i+1:n)))/a(i,i)
		error=error+abs(ui-u(i))
		u(i)=ui
	end do
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
deallocate(a,aa,u,b,bb,r,uold,auxxx)
end program seidel
