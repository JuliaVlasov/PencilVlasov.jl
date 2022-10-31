export Mesh

"""

    Mesh(nc_η1, nc_η2, nc_η3, nc_η4)

4D uniform mesh data.

"""
struct Mesh

    nx1::Int
    nx2::Int
    nx3::Int
    nx4::Int
    x1::Vector{Float64}
    x2::Vector{Float64}
    x3::Vector{Float64}
    x4::Vector{Float64}
    dx1::Float64
    dx2::Float64
    dx3::Float64
    dx4::Float64
    x1min::Float64
    x1max::Float64
    x2min::Float64
    x2max::Float64
    x3min::Float64
    x3max::Float64
    x4min::Float64
    x4max::Float64

    function Mesh(nx1, nx2, nx3, nx4)

        x1min, x1max = 0, 4π
        x2min, x2max = 0, 4π
        x3min, x3max = -6, 6
        x4min, x4max = -6, 6

        x1 = LinRange(x1min, x1max, nx1 + 1)[1:end-1] |> collect
        x2 = LinRange(x2min, x2max, nx2 + 1)[1:end-1] |> collect
        x3 = LinRange(x3min, x3max, nx3 + 1)[1:end-1] |> collect
        x4 = LinRange(x4min, x4max, nx4 + 1)[1:end-1] |> collect

        dx1 = (x1max - x1min) / nx1
        dx2 = (x2max - x2min) / nx2
        dx3 = (x3max - x3min) / nx3
        dx4 = (x4max - x4min) / nx4

        new(
            nx1,
            nx2,
            nx3,
            nx4,
            x1,
            x2,
            x3,
            x4,
            dx1,
            dx2,
            dx3,
            dx4,
            x1min,
            x1max,
            x2min,
            x2max,
            x3min,
            x3max,
            x4min,
            x4max,
        )

    end

end
