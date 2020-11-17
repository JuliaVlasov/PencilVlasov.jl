"""

   poisson!(ρ, mesh, ex, ey)

Solve the equation Δ Φ = - ρ

 ex = ∂ Φ / ∂ x
 ey = ∂ Φ / ∂ y 

WARNING: the ρ array is destroyed

"""
function poisson!( ρ::Array{ComplexF64,2}, 
		   mesh::Mesh, 
		   ex::Array{ComplexF64,2}, 
		   ey::Array{ComplexF64,2} )

    nx1 = length(meshx.x1)
    nx2 = length(meshx.x2)
    kx0 =  2π / (mesh.x1[end] - mesh.x1[1])
    ky0 =  2π / (mesh.x2[end] - mesh.x2[1])

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
