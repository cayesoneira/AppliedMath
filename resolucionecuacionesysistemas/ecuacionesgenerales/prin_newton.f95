program prin_newton
implicit none
real::eps,del,nitmax,x0
external::f,g

call lecnew(eps,del,nitmax,x0)
call iter(f,g,eps,del,nitmax,x0)

end program prin_newton
