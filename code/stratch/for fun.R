N = 1000000
x = runif(N)
y = runif(N)

idx <- which(
  (((x-0.5)^2+y^2)<=0.5^2)&
     ((x^2+(y-1)^2)<=1^2)
              )



16*length(idx)/N
4*atan(2)+16*atan(0.5)-8
