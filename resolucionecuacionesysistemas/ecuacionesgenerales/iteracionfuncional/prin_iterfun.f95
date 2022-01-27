program iterfun
implicit none
real::eps,nitmax,x0
external::f

call leciter(eps,nitmax,x0)
call iter(f,eps,nitmax,x0)

end program iterfun
