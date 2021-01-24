export Maxwell

struct Maxwell

    kx::Vector{Float64}
    ky::Vector{Float64}

    function Maxwell(mesh)

        nx = mesh.nx1
        Lx = mesh.x1max - mesh.x1min
        kx = 2π / Lx .* vcat(0:nx÷2-1, -nx÷2:-1)
        ny = mesh.nx2
        Ly = mesh.x2max - mesh.x2min
        ky = 2π / Ly .* vcat(0:ny÷2-1, -ny÷2:-1)

        new(kx, ky)

    end

end

export faraday!, ampere_maxwell!

"""
    ampere_maxwell!( ex, ey, solver, bz, dt)

```math
Ex^{t+dt} = Ex^{t} + dt \\frac{\\partial Bz}{\\partial y}
```
    
```math
Ey^{t+dt} = Ey^{t} - dt \\frac{\\partial Bz}{\\partial x}
```
"""
function ampere_maxwell!(ex, ey, s::Maxwell, bz::AbstractArray, dt)

    ex .= ex .- dt .* ifft(-1im .* s.ky' .* fft(bz, 2))
    ey .= ey .+ dt .* ifft(-1im .* s.kx .* fft(bz, 1))

end


"""
    faraday!( bz, solver, ex, ey, dt)

```math
Bz^{t+dt} = Bz^{t} + dt (  \\frac{\\partial Ey}{\\partial x}
-   \\frac{\\partial Ex}{\\partial y})
"""
function faraday!(bz, s::Maxwell, ex, ey, dt::Float64)

    dex_dy = ifft(-1im .* s.ky' .* fft(ex, 2))
    dey_dx = ifft(-1im .* s.kx .* fft(ey, 1))

    bz .= bz .- dt .* (dex_dy .- dey_dx)

end
