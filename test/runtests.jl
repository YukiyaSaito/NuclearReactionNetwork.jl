using NuclearReactionNetwork
using Test
using Profile
using StatProfilerHTML

@testset "NuclearReactionNetwork.jl" begin
    # Write your tests here.
    # read_test("/Users/yukiya/Documents/Physics/PhD/ReactionNetwork.jl/moller2003.dat")
    @time networksize, timestep = read_test("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/extent/default.dat")
    # @profilehtml networksize, timestep = read_test("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/extent/default.dat")
    # Profile.print()
    # @test networksize == 13396
    # @test typeof(timestep)==Float64
end

