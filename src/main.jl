using ..ReactionTypes
using ..Network
using ..NetworkSolve
using ..Astro
using ..InOut
using ..NetworkDatas
using LinearAlgebra
using SparseArrays
using Printf
using InteractiveUtils
using Profile
using Interpolations

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
    probdecay_files = ["/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/nuclear/Nubase/betam_nubase2016_moller.dat", "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Converted/Reaclib/moller2003.dat"]
    ncap_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Reaclib_ng.dat"
    alphadecay_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/nuclear/Nubase/alpha_nubase2016.dat"
    photodissociation_files = ["/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Converted/NDI/s1n_ame2016.dat", "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Data/Converted/NDI/s1n_frdm2012.dat"]
    trajectory_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-wind_trajectory"
    # trajectory_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-dynamical-merger_trajectory"
    initial_abundance_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-wind_initx"
    # initial_abundance_file = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/examples/conditions/rprocess-dynamical-merger_initx"

    BLAS.set_num_threads(4)
    println("Reading network boundary...")
    extent = read_boundary(path)
    networksize = get_networksize(extent)
    net_idx = NetworkIndex(extent)

    println("Reading reaction data...")
    reaction_data::ReactionData = initialize_reactions()
    included_reactions = IncludedReactions(false, false, false, false)
    for file in probdecay_files
        read_probdecay!(file, reaction_data, net_idx)
    end
    read_ncap!(ncap_file, reaction_data, net_idx)
    read_alphadecay!(alphadecay_file, reaction_data, net_idx)
    for file in photodissociation_files
        read_photodissociation!(file, reaction_data, net_idx)
    end
    trajectory = read_trajectory(trajectory_file)
    abundance = initialize_abundance(networksize)
    abundance = read_initial_abundance(initial_abundance_file, abundance, net_idx)
    yproposed = Vector{Float64}(undef, length(abundance))
    ydot = initialize_ydot(networksize)
    times = knots(trajectory.temperatures)
    start_time = times[1]
    final_time = times[length(times)]

    # time = Time(start_time, 1e-15, final_time)
    time = Time(start_time, 1e-15, 1e-4)
    if time.current != start_time
        println("WARNING: Using start time of $(time.current)s, instead of start of trajectory ($(start_time)s)")
    end
    if time.stop != final_time
        println("WARNING: Using final time of $(time.stop)s, instead of end of trajectory ($(final_time)s)")
    end
    println("Starting simulation from $(time.current)s to $(time.stop)s")

    jacobian = initialize_jacobian(networksize)

    nd = NetworkData(net_idx, reaction_data, trajectory, abundance, yproposed, ydot, time, jacobian, OutputInfo(false, missing, false, missing), included_reactions)
    initialize_and_fill_sparse_jacobian(networksize, nd)
    fill_probdecay_ydot!(nd)

    dump_each_iteration = false
    try
        SolveNetwork!(nd, dump_each_iteration)
    finally
        res = Result(nd);
        dump_result(res, "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/output/wind-all-NEW-Y.txt")
        println("Wrote results to disk")
    end 
    return networksize, nd.time.step
end