using Test
using PencilArrays


@testset "PencilVlasov" begin

using MPI
using PencilArrays
using LinearAlgebra: transpose!

MPI.Init()
comm = MPI.COMM_WORLD       # we assume MPI.Comm_size(comm) == 12
rank = MPI.Comm_rank(comm)  # rank of local process, in 0:11

# Define MPI Cartesian topology: distribute processes on a 3×4 grid.
topology = MPITopology(comm, (2, 2))

# Let's decompose 3D arrays along dimensions (2, 3).
# This corresponds to the "x-pencil" configuration in the figure.
# This configuration is described by a Pencil object.
dims_global = (42, 31, 29)  # global dimensions of the array
decomp_dims = (2, 3)
pen_x = Pencil(topology, dims_global, decomp_dims)

# We can now allocate distributed arrays in the x-pencil configuration.
Ax = PencilArray{Float64}(undef, pen_x)
fill!(Ax, rank * π)  # each process locally fills its part of the array
parent(Ax)           # parent array holding the local data (here, an Array{Float64,3})
size(Ax)             # size of local part
size_global(Ax)      # total size of the array = (42, 31, 29)

# Create another pencil configuration, decomposing along dimensions (1, 3).
# We could use the same constructor as before, but it's recommended to reuse the
# previous Pencil instead.
pen_y = Pencil(pen_x, decomp_dims=(1, 3))

# Now transpose from the x-pencil to the y-pencil configuration, redistributing
# the data initially in Ax.
Ay = PencilArray{Float64}(undef, pen_y)
transpose!(Ay, Ax)


@test true

end
