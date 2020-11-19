#=
hx = hx - dt ( iffty(-jky ffty(ez) 
	+ dt^2/24( ifftx(-kx^2 fftx ( iffty ( -jky ffty(ez)))))
        + iffty(jky^3 ffty(ez))))

hy = hy + dt ( ifftx(-jkx fftx(ez)) 
	+ dt^2/24 *(iffty(-ky^2 ffty(ifftx(-jkx fftx(ez)))) 
	+ ifftx(jkx^3 fftx(ez))))

ez = ez + dt ( ifftx(-jkx fftx(hy)) - iffty(-jky ffty(hx)) 
	+ dt^2/24 ( 
	      2 iffty( -ky^2 ffty(ifftx(-jkx fftx(hy))))
	      + ifftx( jkx^3 fftx(hy)) 
	      - iffty( jky^3 ffty(hx))))
=#

struct MaxwellPSTD

    kx :: Vector{Float64}
    ky :: Vector{Float64}

    function PSTD( mesh)

        nx = mesh.nx
        kx = vcat(0:nx÷2,-nx÷2:-1)

        ny = mesh.ny
        ky = vcat(0:ny÷2,-ny÷2:-1)

        new( kx, ky )

    end

end

function faraday!( ex, ey, bz, s :: MaxwellPSTD, dt :: Float64)

    @. ex += dt * ( ifft(-1im * s.ky * fft(hz,1), 2)
          + dt^2/24 * (ifft(-kx^2 * fft(ifft(-1im * s.ky * fft(hz, 2), 2), 1), 1))
          + ifft(1im * ky^3 * fft(hz,2), 2))
    
    @. ey -= dt * ( ifft(-1im * s.kx * fft(hz,1), 1) 
    	  + dt^2/24 * (ifft(-ky^2 * fft(ifft(-1im * s.kx * fft(hz, 1), 1), 2), 2) 
    	  + ifft(1im * kx^3 * fft(hz, 1), 1) ) )
end
    
function ampere_maxwell!( ex, ey, bz, s :: MaxwellPSTD, dt :: Float64)

    @. hz -= dt * ( ifftx(-1im * s.kx * fft(ey, 1), 1) 
    	      - ifft(-1im * s.ky * fft(ex, 2), 2) 
    	      - dt^2/24  * ( 
    	      - ifft(-s.ky^2 * fft(ifft(- 1im * s.kx * fft(ey, 1), 1), 2), 2)
    	      + ifft(-s.kx^2 * fft(ifft(- 1im * s.ky * fft(ex, 2), 2), 1), 1)
    	      + ifft(1im * s.kx^3 * fft(ey, 1), 1) - ifft( 1im * s.ky^3 * fft(ex, 2), 2)))

end

