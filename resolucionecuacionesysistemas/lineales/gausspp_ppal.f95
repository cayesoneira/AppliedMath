program gausspp
!DECLARACIONES
implicit none
integer,parameter::clreal=selected_real_kind(p=15,r=307)
real(kind=clreal),allocatable,dimension(:,:)::a,aa
real(kind=clreal),allocatable,dimension(:)::u,b,bb,r,aux
integer,allocatable,dimension(:)::ip
real(kind=clreal)::deter

!DECLARACIONES PARTICULARES
real::piv,cont,suma
integer::k,ipiv,ipi,ipk,i,j,n

!TOMA DE DATOS
print*,'GAUSS PIVOTE PARCIAL'
print*,'Introducir dimensión de la matriz:'
read*,n
print*,n
allocate(a(n,n),aa(n,n),u(n),b(n),bb(n),ip(n),r(n),aux(n))
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

aa=a
bb=b

!MÉTODO
deter=1.
cont=0.
ip=(/(i,i=1,n)/)
do k=1,n-1
	piv=a(ip(k),k)
	ipiv=k
	do i=k+1,n
		if(abs(piv)<abs(a(ip(i),k)))then
			piv=a(ip(i),k)
			ipiv=i
		end if
	end do
	if(abs(piv)<1.e-12) then
		print*,'pivote nulo na etapa: ',k
		print*,'A matriz do sistema e singular!'
		stop
	end if
	if(ipiv/=k) then
		ipk=ip(ipiv)
		ip(ipiv)=ip(k)
		ip(k)=ipk
		cont=cont+1
	else
		ipk=ip(k)
	end if
	deter=deter*piv
	do i=k+1,n
		ipi=ip(i)
		a(ipi,k)=a(ipi,k)/piv
			do j=k+1,n
				a(ipi,j)=a(ipi,j)-a(ipi,k)*a(ipk,j)
			end do
		b(ipi)=b(ipi)-a(ipi,k)*b(ipk)
	end do
	print*,'·Etapa',k
	print*,'Vector de permutaciones',ip
end do
piv=a(ip(n),n)
if(abs(piv)<1.e-12) then
	print*,'Pivote nulo, etapa',n
	print*,'Matriz singular'
	stop
end if
do i=1,n
	do j=1,n
		if(i>j)then
			a(i,j)=a(i,j)-a(i,j)
		end if
	end do
end do
deter=deter*piv*(-1)**cont
print*,'Determinante de A:',deter
print*,'Matriz U triangular superior:'
do i=1,n
	print*,a(i,1:n)
end do

!RESOLUCIÓN DEL SISTEMA
do i=n,1,-1
	suma=0.
	do j=i+1,n
		suma=suma+a(ip(i),j)*u(j)
	end do
	u(i)=(b(ip(i))-suma)/a(ip(i),i)
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
deallocate(a,aa,u,b,bb,ip,r,aux)
end program gausspp
