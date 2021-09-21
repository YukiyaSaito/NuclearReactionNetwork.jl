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
    probdecay_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/nuclear/Nubase/betam_nubase2016_moller.dat"
    ncap_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Reaclib_ng.dat"
    alphadecay_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/nuclear/Nubase/alpha_nubase2016.dat"
    photodissociation_files = ["/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Converted/NDI/s1n_ame2016.dat", "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Converted/NDI/s1n_frdm2012.dat"]
    trajectory_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-wind_trajectory"
    # trajectory_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-dynamical-merger_trajectory"
    initial_abundance_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-wind_initx"
    # initial_abundance_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-dynamical-merger_initx"

    BLAS.set_num_threads(4)
    reaction_data::ReactionData = initialize_reactions()
    read_probdecay!(probdecay_file, reaction_data)
    read_ncap!(ncap_file, reaction_data)
    read_alphadecay!(alphadecay_file, reaction_data)
    for file in photodissociation_files
        read_photodissociation!(file, reaction_data)
    end
    trajectory = read_trajectory(trajectory_file)
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
    abundance = read_initial_abundance(initial_abundance_file, abundance, net_idx)
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
    dump_each_iteration = true
    try
        println("Solving network")
        SolveNetwork!(abundance, jacobian, reaction_data, ydot, time, net_idx, trajectory, dump_each_iteration)
    finally
        res = Result(abundance, net_idx);
        # dump_result(res, "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/output/wind-beta+ncap+alpha+photo-Y.txt")
        dump_result(res, "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/output/wind-beta-Y.txt")
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