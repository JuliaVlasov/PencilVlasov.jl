@testset "Maxwell PSTD" begin

x = LinRange( 0, 4π, 1024)
y = LinRange( 0, 4π, 1024)

omega      =   sqrt(8*pi*pi)
bz         = - cos(2*pi*x)*cos(2*pi*y)*cos(omega*t)
ex         =   cos(2*pi*x)*sin(2*pi*y)*sin(omega*t)*2*pi/ omega
ey         = - sin(2*pi*x)*cos(2*pi*y)*sin(omega*t)*2*pi/ omega
diff(ex,t) =   diff(bz,y)
diff(ey,t) = - diff(bz,x)
diff(bz,t) =   diff(ex,y) - diff(ey,x)

omega = sqrt(2*pi*pi)
ez =  cos(pi*x)*cos(pi*y)*cos(omega*t)
hx =  cos(pi*x)*sin(pi*y)*sin(omega*t) * pi / omega
hy = -sin(pi*x)*cos(pi*y)*sin(omega*t) * pi / omega

diff(hx,t) = - diff(ez,y)
diff(hy,t) =   diff(ez,x)
diff(ez,t) =   diff(hy,x) - diff(hx,y)

  d(d2hz/dy2 + d2hz/dx2)/dy
- d(d2hz/dx2 + d2hz/dy2)/dx
  d3ex/dy3 - d3dey/dx3 + d3ex/dydx2 - d3ey/dxdy2 

end
