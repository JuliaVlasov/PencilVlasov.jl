export poisson!

"""

   poisson!(ex, ey, mesh, ρ)

Solve the equation Δ Φ = - ρ

 ex = ∂ Φ / ∂ x
 ey = ∂ Φ / ∂ y 

WARNING: the ρ array is destroyed

"""
function poisson!( ex::Array{ComplexF64,2}, 
		   ey::Array{ComplexF64,2},
		   mesh::Mesh, ρ::Array{ComplexF64,2})

    nx1 = length(mesh.x1)
    nx2 = length(mesh.x2)
    kx0 = 2π / 4π
    ky0 = 2π / 4π

    fft!(ρ,[1,2])
    
    kx = kx0 * vcat(0:nx1÷2-1,-nx1÷2:-1)
    kx[1] = 1
    ky = ky0 * vcat(0:nx2÷2-1,-nx2÷2:-1)
    kx[1] = 1

    for i = 1:nx1
       kx2 = kx[i]*kx[i]
       for j =  1:nx2
          k2 = kx2 + ky[j]*ky[j]
          ex[i,j] = -1im * kx[i]/k2 * ρ[i,j]
	  ey[i,j] = -1im * ky[j]/k2 * ρ[i,j]
       end 
    end

    ifft!(ex,[1,2])
    ifft!(ey,[1,2])

end
