program prin_dicoto
implicit none
real::a,b,del,eps,z
external::f

call lecdic(a,b,del,eps)
call dicoto(f,a,b,del,eps,z)

stop
end program prin_dicoto
