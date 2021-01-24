export Mesh

"""

    Mesh(nc_η1, nc_η2, nc_η3, nc_η4)

4D uniform mesh data.

"""
struct Mesh

    nx1::Any
    nx2::Any
    nx3::Any
    nx4::Any
    x1::Any
    x2::Any
    x3::Any
    x4::Any
    dx1::Any
    dx2::Any
    dx3::Any
    dx4::Any
    x1min::Any
    x1max::Any
    x2min::Any
    x2max::Any
    x3min::Any
    x3max::Any
    x4min::Any
    x4max::Any

    function Mesh(nx1, nx2, nx3, nx4)

        x1min, x1max = 0, 4π
        x2min, x2max = 0, 4π
        x3min, x3max = -6, 6
        x4min, x4max = -6, 6

        x1 = LinRange(0, 4π, nx1 + 1)[1:end-1] |> collect
        x2 = LinRange(0, 4π, nx2 + 1)[1:end-1] |> collect
        x3 = LinRange(-6, 6, nx3 + 1)[1:end-1] |> collect
        x4 = LinRange(-6, 6, nx4 + 1)[1:end-1] |> collect

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
