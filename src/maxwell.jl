#=
ez = ez + dt ( ifftx(-jkx fftx(hy)) - iffty(-jky ffty(hx)) 
	+ dt^2/24 ( 
	      2 iffty( -ky^2 ffty(ifftx(-jkx fftx(hy))))
	      + ifftx( jkx^3 fftx(hy)) 
	      - iffty( jky^3 ffty(hx))))
=#

export Maxwell

struct Maxwell

    kx :: Vector{Float64}
    ky :: Vector{Float64}

    function Maxwell( mesh)

        nx = mesh.nx1
        kx = 0.5 .* vcat(0:nx÷2-1,-nx÷2:-1)
        ny = mesh.nx2
        ky = 0.5 .* vcat(0:ny÷2-1,-ny÷2:-1)

        new( kx, ky )

    end

end

export faraday!, ampere_maxwell!

"""
    ampere_maxwell!( ex, ey, solver, bz, dt)

```math
ex .+= dt * ( ifft(- 1im * s.ky * ffty(hz), 2)
            + dt^2/24 ( ifftx(-kx^2 fftx( iffty ( -jky ffty(hz)))))
            + iffty(1im * ky.^3 * fft(hz,2))))
    
    ey = ey - dt ( ifftx(-jkx fftx(hz)) 
    	+ dt^2/24 *(iffty(-ky^2 ffty(ifftx(-jkx fftx(hz)))) 
    	+ ifftx(jkx^3 fftx(hz))))
    

```
"""
function ampere_maxwell!( ex, ey, s :: Maxwell, bz :: AbstractArray, dt :: Float64)

    ex .+= dt .* ifft(-1im .* s.ky .* fft(bz,2))
    ey .-= dt .* ifft(-1im .* s.kx .* fft(bz,1))

        #.+ dt^2/24 .* (ifftx(-s.kx.^2 .* fftx(iffty(-1im .* s.ky .* ffty(bz)))) 
        #.+ iffty(1im .* s.ky.^3 .* ffty(bz))))
    
        #.+ dt^2/24 .* (iffty(-s.ky.^2 .* ffty(ifftx(-1im .* s.kx .* fftx(bz)))) 
        #.+ ifftx(1im .* s.kx.^3 .* fftx(bz))))
end

    
"""
    hz = hz - dt (  ifftx(-jkx fftx(ey)) 
    	      - iffty(-jky ffty(ex)) 
    	      - dt^2/24 ( 
    	      - iffty(-ky^2 ffty(ifftx(-jkx fftx(ey))))
    	      + ifftx(-kx^2 fftx(iffty(-jky ffty(ex))))
    	      + ifftx(jkx^3 fftx(ey)) 
                  - iffty(jky^3 ffty(ex))))
"""
function faraday!( bz, s :: Maxwell, ex, ey, dt :: Float64)

    dex_dy = ifft(-1im .* s.ky .* fft(ex,2))
    dey_dx = ifft(-1im .* s.kx .* fft(ey,1))

    bz .+= dt .* ( dex_dy .- dey_dx )
    	        #.- dt^2/24  .* ( 
    	        #.- iffty(-s.ky.^2 .* ffty(ifftx(-1im .* s.kx .* fftx(ey))))
    	        #.+ ifftx(-s.kx.^2 .* fftx(iffty(-1im .* s.ky .* ffty(ex))))
    	        #.+ ifftx(1im .* s.kx.^3 .* fftx(ey)) 
                #.- iffty(1im .* s.ky.^3 .* ffty(ex))))

end
