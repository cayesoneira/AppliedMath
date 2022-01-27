subroutine lecdic(a,b,del,eps)
implicit none
real,intent(inout)::a,b,del,eps

print*,'Extremo izquierdo del intervalo.'
read*,a
print*,a

print*,'Extremo derecho del intervalo.'
read*,b
print*,b

if(a>b)then
	print*,'Mal introducidos los extremos.'
	stop
end if

print*,'Delta de parada: distancia de F(x0) al 0.'
read*,del
print*,del

print*,'Épsilon: distancia de x0 a la raíz.'
read*,eps
print*,eps

return
end subroutine lecdic
