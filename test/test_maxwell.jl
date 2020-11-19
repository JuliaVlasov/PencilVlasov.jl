using LinearAlgebra
using PencilVlasov

@testset "Maxwell PSTD" begin

"""
    solution!( h, dt, time)

Set he electric fields at t=time and the magnetic field
at t=time+0.5dt

"""
function solution!( ex, ey, bz, 
                    mesh,
                    dt   :: Float64,
                    time :: Float64)

	nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    ω = sqrt(2)

    xn = range(0, step=dx, length=nx) 
    yn = range(0, step=dy, length=ny)

    xc = 0.5 .* (xn[1:end-1] .+ xn[ 2:end])
    yc = 0.5 .* (yn[1:end-1] .+ yn[ 2:end])

    for j=1:ny, i=1:nx
        ex[i,j] = + 1 / ω * cos(xc[i]) * sin(yn[j]) * sin(ω*time)
        ey[i,j] = - 1 / ω * sin(xn[i]) * cos(yc[j]) * sin(ω*time)
        bz[i,j] = - cos(xc[i]) * cos(yc[j]) * cos(ω*(time+0.5dt))
    end

end

function error( ex, ey, bz, sol_ex, sol_ey, sol_bz)

	nx, ny = size(ex)
    err_l2 = norm(ex .- sol_ex) / (nx*ny)
    println(" error ex = $err_l2 ")
    err_l2 = norm(ey .- sol_ey) / (nx*ny)
    println(" error ey = $err_l2 ")
    err_l2 = norm(bz .- sol_bz) / (nx*ny)
    println(" error bz = $err_l2 ")

end

function test_maxwell_pstd( )

    kx, ky =   2,   2
    nx, ny = 256, 256

    mesh = Mesh( nx, ny, nvx, nvy)

    dx = mesh.dx
    dy = mesh.dy

    cfl    = 0.1
    dt     = cfl / sqrt(1/(dx*dx)+1/(dy*dy))
    tfinal = 4π
    @show nstep  = trunc(Int64,tfinal/dt)

    time  = 0.

    ex = zeros(nx,ny)
    ey = zeros(nx,ny)
    bz = zeros(nx,ny)

	sol_ex = zeros(nx,ny)
    sol_ey = zeros(nx,ny)
    sol_bz = zeros(nx,ny)

    solver = Maxwell( mesh )

    solution!(f, dt, time)

    for istep = 1:nstep # Loop over time

        ampere_maxwell!(f, dt) # E^{n} -> E^{n+1)

        periodic_bc!(f, dt)

        time += 0.5dt

        faraday!(f, dt)   # B^{n+1/2} -> B^{n+3/2}

        time += 0.5dt

    end 

    solution!(h, dt, time)
    error( f, h)

end

@time test_maxwell_pstd()

end
