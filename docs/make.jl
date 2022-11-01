using Documenter
using PencilVlasov

makedocs(
    sitename = "PencilVlasov.jl",
    authors = "Julia Vlasov",
    format = Documenter.HTML(),
    modules = [PencilVlasov]
)

deploydocs(
    repo="github.com/JuliaVlasov/PencilVlasov.jl",
    devbranch="master",
)
