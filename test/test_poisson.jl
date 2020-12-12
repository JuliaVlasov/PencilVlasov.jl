using PencilVlasov

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
