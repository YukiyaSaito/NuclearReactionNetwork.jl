using ..ReactionTypes
using ..Network
using ..NetworkSolve
using ..Astro
using LinearAlgebra
using SparseArrays
using Printf
using InteractiveUtils
using Profile
# using ..InOut

export main
export read_test

function main()
    reactant = [[2,8]]
    product = [[3,7],[3,6],[0,1]]
    rate = 1.0
    average_number = [0.8, 0.2, 0.2]
    somedecay = ProbDecay(reactant, product, rate, average_number)
    println(somedecay)
end

function read_test(path::String)
    BLAS.set_num_threads(4)
    reaction_data::ReactionData = initialize_reactions()
    read_probdecay!("/Users/yukiya/Documents/Physics/PhD/prism/input/nuclear/Nubase/betam_nubase2016_moller.dat", reaction_data)
    read_ncap!("/Users/yukiya/Documents/Physics/PhD/ReactionNetwork.jl/Reaclib_ng.dat", reaction_data)
    trajectory = read_trajectory("/Users/yukiya/Documents/Physics/PhD/prism/input/examples/conditions/rprocess-dynamical-merger_trajectory")
    # display(reaction_data.neutroncapture[[[0,1],[55,110]]])
    extent = read_boundary(path)
    # println(extent)
    networksize = get_networksize(extent)
    zn_to_index_dict, mass_vector = zn_to_index(extent)
    # display(mass_vector)
    # println(zn_to_index_dict[[0,1]])
    # println(typeof(zn_to_index_dict))
    # println(zn_to_index_dict)
    abundance = initialize_abundance(networksize)
    abundance = read_initial_abundance("/Users/yukiya/Documents/Physics/PhD/prism/input/examples/conditions/rprocess-dynamical-merger_initx",abundance,zn_to_index_dict)
    display(check_mass_fraction_unity(abundance,mass_vector))
    display(abundance)
    ydot = initialize_ydot(networksize)
    fill_probdecay_ydot!(ydot,abundance,reaction_data,zn_to_index_dict)
    # display(ydot)
    timestep, current_time = initialize_timestep()
    jacobian = initialize_and_fill_sparse_jacobian(networksize,abundance,reaction_data,zn_to_index_dict,timestep)
    # display(jacobian)
    time_limit::Float64 = 20.0
    # display(abundance[zn_to_index_dict[[54,93]]])
    # display(reaction_data.probdecay[[[0,1]]].rate)
    SolveNetwork!(abundance, jacobian, reaction_data, ydot, timestep, current_time, mass_vector, time_limit, zn_to_index_dict)
    # Profile.print()
    # ydelta = Vector{Float64}(undef,size(abundance)[1])
    # yproposed = Vector{Float64}(undef,size(abundance)[1])
    # F = lu(jacobian)
    # while current_time < 1e-2
    #     current_time = newton_raphson_iteration!(abundance,yproposed,jacobian, F ,ydot, ydelta,timestep, current_time, mass_vector)
    #     # @printf "%e\n" current_time
    # end
    
    # display(check_mass_fraction_unity(abundance,mass_vector))
    display(abundance)
    # display(yproposed)
    # display(yproposed - abundance)
    # display(ydelta)
    # println(sum(abundance))
    # println(timestep)
    # println(typeof(jacobian))
    # println(typeof(reaction_data))
    return networksize, timestep
    # println(reaction_data.probdecay[[[57,109]]])
end