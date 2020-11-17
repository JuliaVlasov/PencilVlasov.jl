export Mesh

"""

    Mesh(nc_η1, nc_η2, nc_η3, nc_η4)

4D uniform mesh data.

"""
struct Mesh

   x1
   x2
   x3
   x4

   function Mesh( nx1, nx2, nx3, nx4)

       x1 = LinRange(  0, 4π, nx1+1)[1:end-1]
       x2 = LinRange(  0, 4π, nx2+1)[1:end-1]
       x3 = LinRange( -6,  6, nx3+1)[1:end-1]
       x4 = LinRange( -6,  6, nx4+1)[1:end-1]

   end

end
