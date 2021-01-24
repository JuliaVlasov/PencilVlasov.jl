using PencilVlasov
using FFTW
using Test

include("test_poisson.jl")
include("test_maxwell.jl")

@testset "Vlasov-Maxwell 2d2v" begin

     nx, ny = 32, 32
     nvx, nvy = 64, 64
     mesh = Mesh( nx, ny, nvx, nvy)
     maxwell = Maxwell( mesh )
     fn = zeros(ComplexF64, (nx, ny, nvx, nvy))
     ft = zeros(ComplexF64, (nvx, nvy, nx, ny))
#    for l=1:n_l, k=1:n_k, j=1:n_j, i=1:n_i
#        x  = η1_min+(i-1)*delta_η1
#        y  = η2_min+(j-1)*delta_η2
#        vx = η3_min+(k-1)*delta_η3
#        vy = η4_min+(l-1)*delta_η4
#        v2 = vx*vx+vy*vy
#        f[i,j,k,l] = ( 1 + ϵ * cos(kx * x) * cos( ky * y )) * exp(-v2)
#    end
#
#
#    compute_charge(vlasov4d)
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
