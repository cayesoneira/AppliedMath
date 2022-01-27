function f(z) result (res)
implicit none
real,intent(in)::z
real::res

res=exp(-z)-z

end function f
