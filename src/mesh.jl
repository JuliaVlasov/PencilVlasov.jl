export Mesh

"""

    Mesh(nc_η1, nc_η2, nc_η3, nc_η4)

4D uniform mesh data.

"""
struct Mesh

   nx1
   nx2
   nx3
   nx4
   x1
   x2
   x3
   x4
   dx1
   dx2
   dx3
   dx4

   function Mesh( nx1, nx2, nx3, nx4)

       x1min, x1max =  0, 4π
       x2min, x2max =  0, 4π
       x3min, x3max = -6, 6
       x4min, x4max = -6, 6

       x1 = LinRange(  0, 4π, nx1+1)[1:end-1] |> collect
       x2 = LinRange(  0, 4π, nx2+1)[1:end-1] |> collect
       x3 = LinRange( -6,  6, nx3+1)[1:end-1] |> collect
       x4 = LinRange( -6,  6, nx4+1)[1:end-1] |> collect

       dx1 = (x1max - x1min) / nx1
       dx2 = (x2max - x2min) / nx2
       dx3 = (x3max - x3min) / nx3
       dx4 = (x4max - x4min) / nx4

       new( nx1, nx2, nx3, nx4, x1, x2, x3, x4, dx1, dx2, dx3, dx4)

   end

end
