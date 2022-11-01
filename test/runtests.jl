using PencilVlasov
using FFTW
using Test

@testset "Poisson 2D on rectangular grid" begin

    mesh = Mesh( 64, 128, 128, 128)
    x = mesh.x1
    y = mesh.x2
    
    ex   = zeros(ComplexF64, (mesh.nx1, mesh.nx2))
    ey   = zeros(ComplexF64, (mesh.nx1, mesh.nx2))
    ρ    = zeros(ComplexF64, (mesh.nx1, mesh.nx2))
    
    ρ   .= - 2 * sin.(x) .* cos.(y')

    poisson!( ex, ey, mesh, ρ)

    @test maximum(abs.( ex .- (cos.(x) .* cos.(y')))) < 1e-14
    @test maximum(abs.( ey .+ (sin.(x) .* sin.(y')))) < 1e-14

end

@testset "Vlasov-Maxwell 2d2v" begin

     ϵ = 0.01
     kx = 0.5
     ky = 0.5
     nx, ny = 32, 32
     nvx, nvy = 64, 64
     mesh = Mesh( nx, ny, nvx, nvy)

     vlasov = Vlasov(mesh)

     for l=1:mesh.nx4, k=1:mesh.nx3, j=1:mesh.nx2, i=1:mesh.nx1
         x  = mesh.x1min+(i-1)*mesh.dx1
         y  = mesh.x2min+(j-1)*mesh.dx2
         vx = mesh.x3min+(k-1)*mesh.dx3
         vy = mesh.x4min+(l-1)*mesh.dx4
         v2 = vx*vx+vy*vy
         vlasov.fn[i,j,k,l] = ( 1 + ϵ * cos(kx * x) * cos( ky * y )) * exp(-v2)
     end


     p = (3, 4, 1, 2)
     @test vlasov.ft != permutedims(vlasov.fn, p)

     transposexv!(vlasov)

     @test vlasov.ft == permutedims(vlasov.fn, p)
     @test vlasov.fn == permutedims(vlasov.ft, p)


     ρ = complex.(compute_charge(vlasov))
     ex = zeros(ComplexF64, (mesh.nx1, mesh.nx2))
     ey = zeros(ComplexF64, (mesh.nx1, mesh.nx2))

     poisson!(ex, ey, mesh, ρ)
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
