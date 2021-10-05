using NuclearReactionNetwork
using Test
using Profile
using StatProfilerHTML

@testset "NuclearReactionNetwork.jl" begin
    # Write your tests here.
    # read_test("/Users/yukiya/Documents/Physics/PhD/ReactionNetwork.jl/moller2003.dat")
    @time read_test("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/ControlFiles/wind-beta+ncap+alpha+photo.json")
    # Profile.init(n=5*10^7, delay=0.005)
    # @profilehtml read_test("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/ControlFiles/wind-beta+ncap+alpha+photo.json")
    # Profile.print()
    # @test networksize == 13396
    # @test typeof(timestep)==Float64
end

