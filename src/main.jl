using ..ReactionTypes
using ..Network
using ..NetworkSolve
using ..Astro
using ..InOut
using LinearAlgebra
using SparseArrays
using Printf
using InteractiveUtils
using Profile

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
    read_probdecay!("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/nuclear/Nubase/betam_nubase2016_moller.dat", reaction_data)
    read_ncap!("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Reaclib_ng.dat", reaction_data)
    read_alphadecay!("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/nuclear/Nubase/alpha_nubase2016.dat", reaction_data)
    # trajectory = read_trajectory("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-dynamical-merger_trajectory")
    trajectory = read_trajectory("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-wind_trajectory")
    # display(reaction_data.neutroncapture[[[0,1],[55,110]]])
    extent = read_boundary(path)
    # println(extent)
    networksize = get_networksize(extent)
    net_idx = NetworkIndex(extent)
    # display(mass_vector)
    # println(zn_to_index_dict[[0,1]])
    # println(typeof(zn_to_index_dict))
    # println(zn_to_index_dict)
    abundance = initialize_abundance(networksize)
    # abundance = read_initial_abundance("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-dynamical-merger_initx", abundance, net_idx)
    abundance = read_initial_abundance("/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-wind_initx", abundance, net_idx)
    # display(check_mass_fraction_unity(abundance, net_idx.mass_vector))
    # display(abundance)
    ydot = initialize_ydot(networksize)
    fill_probdecay_ydot!(ydot, abundance, reaction_data, net_idx)
    # display(ydot)
    # time = Time()
    time = Time(0.0, 1e-15, 10.0)
    jacobian = initialize_and_fill_sparse_jacobian(networksize, abundance, reaction_data, trajectory, net_idx, time)
    # display(jacobian)
    # display(abundance[zn_to_index_dict[[54,93]]])
    # display(reaction_data.probdecay[[[0,1]]].rate)
    dump_each_iteration = false
    try
        println("Solving network")
        SolveNetwork!(abundance, jacobian, reaction_data, ydot, time, net_idx, trajectory, dump_each_iteration)
    finally
        res = Result(abundance, net_idx);
        dump_result(res, "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/output/wind-beta+ncap+alpha-Y.txt")
        println("Wrote results to disk")
    end 
    # Profile.print()
    # ydelta = Vector{Float64}(undef,size(abundance)[1])
    # yproposed = Vector{Float64}(undef,size(abundance)[1])
    # F = lu(jacobian)
    # while current_time < 1e-2
    #     current_time = newton_raphson_iteration!(abundance,yproposed,jacobian, F ,ydot, ydelta,timestep, current_time, mass_vector)
    #     # @printf "%e\n" current_time
    # end

    # display(check_mass_fraction_unity(abundance,mass_vector))
    # display(abundaneare)
    # display(yproposed)
    # display(yproposed - abundance)
    # display(ydelta)
    # println(sum(abundance))
    # println(timestep)
    # println(typeof(jacobian))
    # println(typeof(reaction_data))
    return networksize, time.step
    # println(reaction_data.probdecay[[[57,109]]])
end