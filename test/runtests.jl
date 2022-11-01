using PencilVlasov
using FFTW
using Test

include("test_poisson.jl")

@testset "Vlasov-Maxwell 2d2v" begin

     ϵ = 0.01
     kx = 0.5
     ky = 0.5
     nx, ny = 32, 32
     nvx, nvy = 64, 64
     mesh = Mesh( nx, ny, nvx, nvy)
     fn = zeros(ComplexF64, (nx, ny, nvx, nvy))
     ft = zeros(ComplexF64, (nvx, nvy, nx, ny))

     for l=1:mesh.nx4, k=1:mesh.nx3, j=1:mesh.nx2, i=1:mesh.nx1
         x  = mesh.x1min+(i-1)*mesh.dx1
         y  = mesh.x2min+(j-1)*mesh.dx2
         vx = mesh.x3min+(k-1)*mesh.dx3
         vy = mesh.x4min+(l-1)*mesh.dx4
         v2 = vx*vx+vy*vy
         fn[i,j,k,l] = ( 1 + ϵ * cos(kx * x) * cos( ky * y )) * exp(-v2)
     end

     p = (3, 4, 1, 2) # permutation

     @test ft != permutedims(fn, p)

     permutedims!(ft, fn, p)

     @test ft == permutedims(fn, p)
     @test fn == permutedims(ft, p)

     vlasov = Vlasov(mesh)

     ρ = compute_charge(vlasov)

#    solve(poisson, ex, ey, rho)
#    bz = zeros(nx,ny)
#    
#    transposevx!(vlasov4d)
#    advection_x1!(vlasov4d,0.5dt)
#    advection_x2!(vlasov4d,0.5dt)
#    
#    for istep = 1:nstep
#    
#         transposexv!(vlasov4d)
#    
#         compute_current(vlasov4d)
#    
#         solve_ampere!(maxwell, dt)
#    
#         advection_x3x4(vlasov4d, dt)
#    
#         transposevx!(vlasov4d)
#    
#         advection_x1(vlasov4d,dt)
#         advection_x2(vlasov4d,dt)
#    
#    end
#
end
