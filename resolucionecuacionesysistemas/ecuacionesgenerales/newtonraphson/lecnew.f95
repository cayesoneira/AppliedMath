subroutine lecnew(eps,del,nitmax,x0)
implicit none
real,intent(inout)::eps,del,nitmax,x0

print*,'Escribir el épsilon de parada:'
read*,eps
print*,eps

print*,'Escribir o delta de parada (residuo da imaxe)'
read*,del
print*,del

print*,'Escribir el número máximo de iteraciones:'
read*,nitmax
print*,nitmax

print*,'Escribir el punto inicial:'
read*,x0
print*,x0

end subroutine lecnew
