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

struct Maxwell

    kx :: Vector{Float64}
    ky :: Vector{Float64}

    function PSTD( meshx, meshy)

        nx = meshx.length
        kx = vcat(0:nx÷2,-nx÷2:-1)

        ny = meshy.length
        ky = vcat(0:ny÷2,-ny÷2:-1)

        new( kx, ky )

    end

end
#=

function ampere!( f :: MaxwellFieldsTE, s :: PSTD, dt :: Float64)


    f.ex .+= dt * ( ifft(- 1im * s.ky * ffty(hz), 2)
            + dt^2/24 ( ifftx(-kx^2 fftx ( iffty ( -jky ffty(hz)))))
            + iffty(1im * ky.^3 * fft(hz,2))))
    
    ey = ey - dt ( ifftx(-jkx fftx(hz)) 
    	+ dt^2/24 *(iffty(-ky^2 ffty(ifftx(-jkx fftx(hz)))) 
    	+ ifftx(jkx^3 fftx(hz))))
    
    hz = hz - dt (  ifftx(-jkx fftx(ey)) 
    	      - iffty(-jky ffty(ex)) 
    	      - dt^2/24 ( 
    	      - iffty(-ky^2 ffty(ifftx(-jkx fftx(ey))))
    	      + ifftx(-kx^2 fftx(iffty(-jky ffty(ex))))
    	      + ifftx(jkx^3 fftx(ey)) 
                  - iffty(jky^3 ffty(ex))))

end
=#

