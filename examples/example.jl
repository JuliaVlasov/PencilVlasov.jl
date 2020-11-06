using Test
using PencilArrays
using MPIClusterManagers
using Distributed
using LinearAlgebra: transpose!

manager = MPIManager(np=12)
addprocs(manager)

println("Added procs $(procs())")

@everywhere import MPI

println("Running hello as part of a Julia cluster")
@mpi_do manager (include("hello.jl"); do_hello())

println("Try the PencilArrays example")
@mpi_do manager (include("pencil_arrays.jl"); test_pencil_arrays())

println("Exiting")
exit()
