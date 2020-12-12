using LinearAlgebra
using PencilVlasov

@testset "Maxwell 2D PSTD on rectangular grid" begin

function solution!( ex, ey, mesh :: Mesh, time :: Float64)

    nx, ny = mesh.nx1, mesh.nx2
    dx, dy = mesh.dx1, mesh.dx2
    ω = sqrt(2)

    x = mesh.x1
    y = mesh.x2

    ex .= + 1 ./ ω .* cos.(x) .* sin.(y) .* sin.(ω*time)
    ey .= - 1 ./ ω .* sin.(x) .* cos.(y) .* sin.(ω*time)

end

function solution!( bz, mesh :: Mesh, time :: Float64)

    nx, ny = mesh.nx1, mesh.nx2
    dx, dy = mesh.dx1, mesh.dx2
    ω = sqrt(2)

    x = mesh.x1
    y = mesh.x2

    bz .= - cos.(x) .* cos.(y) .* cos.(ω*time)

end

kx, ky =   2,   2
nx, ny = 128, 128

mesh = Mesh( nx, ny, nx, ny)

dx = mesh.dx1
dy = mesh.dx2

dt = (dx+dy) / (π * sqrt(2))

println(" dt = $dt ")

ex = zeros(ComplexF64, (nx,ny))
ey = zeros(ComplexF64, (nx,ny))
bz = zeros(ComplexF64, (nx,ny))

sol_ex = zeros(nx,ny)
sol_ey = zeros(nx,ny)
sol_bz = zeros(nx,ny)

solver = Maxwell( mesh )

time  = 0.

time = - 0.5dt
solution!(bz, mesh, time) # Initialize B^{n+1/2}
solution!(sol_bz, mesh, time)
err = maximum(abs.(real(bz) .- sol_bz))
println(" error bz = $err ")

ampere_maxwell!(ex, ey, solver, bz, dt) # E^{n} -> E^{n+1)

time += 0.5dt
solution!(sol_ex, sol_ey, mesh, time)

err = maximum(abs.(real(ex) .- sol_ex))
println(" error ex = $err ")

err = maximum(abs.(real(ey) .- sol_ey))
println(" error ey = $err ")

faraday!(bz, solver, ex, ey, dt)   # B^{n-1/2} -> B^{n+1/2}

time += 0.5dt
solution!(sol_bz, mesh, time)

err = maximum(abs.(real(bz) .- sol_bz))
println(" error bz = $err ")

@test true

end
