program householder
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa
real(kind=clreal),allocatable,dimension(:)::u,b,bb,r,c,auxx,w
real(kind=clreal)::deter,s2,alfa,beta,p,q,aux,suma
integer::k,n,j,i

!TOMA DE DATOS
print*,'HOUSEHOLDER'
print*,'Introducir orden de la matriz:'
read*,n
print*,n
allocate(a(n,n),aa(n,n),u(n),b(n),bb(n),r(n),auxx(n),c(n),w(n))
print*,'Introducir la matriz A por filas:'
do i=1,n
	read*,a(i,1:n)
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
	s2=sum(a(k:n,k)*a(k:n,k))
	
	if(s2<1.e-12)then
		print*,'Vector nulo, etapa',k
		print*,'Matriz singular'
		stop
	else if(sum(a(k+1:n,k)*a(k+1:n,k))<1.e-12)then
		print*,'Etapa',k,'no necesita eliminación'
		deter=deter*a(k,k)
		cycle
	else
		if(a(k,k)>=1.e-12)then
			alfa=-sqrt(s2)
		else
			alfa=sqrt(s2)
		end if
		beta=alfa*(alfa-a(k,k))
		a(k,k)=a(k,k)-alfa
	end if
	do i=1,n
		w(i)=0
	end do
	print*,'Vector w'
	do j=k,n
		print*,a(j,k)
		w(j)=a(j,k)
	end do
	print*,'Norma de w:'
	print*,sqrt(dot_product(w,w))
	
	do j=k+1,n
		q=sum(a(k:n,k)*a(k:n,j))
		p=q/beta
		a(k:n,j)=a(k:n,j)-p*a(k:n,k)
	end do
	q=sum(a(k:n,k)*b(k:n))
	p=q/beta
	b(k:n)=b(k:n)-p*a(k:n,k)
	a(k,k)=alfa
	deter=-deter*a(k,k)
	print*,'·Etapa',k
	print*,'Alfa',alfa
	print*,'Beta',beta
	do j=1,n
		c(j)=a(j,k)
	end do
	c(1)=c(1)-alfa
end do
if(abs(a(n,n))<1.e-12)then
	print*,'Elemento diagonal ',n,' nulo'
	print*,'La matriz del sistema es singular'
	stop
end if
deter=deter*a(n,n)
print*,'Determinante:'
print*,deter
do i=1,n
	do j=1,n
		if(i>j)then
			a(i,j)=a(i,j)-a(i,j)
		end if
	end do
end do

print*,'Matriz R:'
do i=1,n
	print*,a(i,1:n)
end do
print*,'Vector b modificado:'
print*,b

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
auxx=0.
do j=1,n
	auxx=auxx+a(:,j)*u(j)
end do
r=auxx-b
print*,'Residuo:'
print*,r
print*,'Norma del residuo:'
print*,sqrt(dot_product(r,r))

!FINAL
deallocate(a,aa,u,b,bb,r,auxx,w)
end program householder
