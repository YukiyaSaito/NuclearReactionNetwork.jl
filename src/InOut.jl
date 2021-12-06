"""
    InOut

Handles everything to do with input/output for the simulation.
"""
module InOut

using DelimitedFiles
using DataStructures
using JSON
using JLD2
using SparseArrays
using StaticArrays
using Pardiso
using ..LinearInterpolations
using ..Astro
using ..ReactionTypes
using ..Network
using ..NetworkDatas
using ..LinearSolvers

export Result
export dump_result
export dump_iteration
export save_checkpoint
export initialize_network_data

# TODO: Should this be a matrix rather than a series of vectors of the same size (we would lose type inhomogeneity)?
"""
    Result

Holds the results of the network calculations.

# Fields:
- `proton_nums::Vector{Int}`: The proton number of each species.
- `neutron_nums::Vector{Int}`: The neutron number of each species.
- `abundance::Vector{Float64}`: The abundance of each species.
"""
struct Result
    """The proton number of each species."""
    proton_nums::Vector{Int}
    """The neutron number of each species."""
    neutron_nums::Vector{Int}
    """The abundance of each species."""
    abundance::Vector{Float64} # TODO: Add more things that we can output, like PRISM does
end

"""
    Result(nd::NetworkData)::Result

Constructs a [`Result`](@ref) from the data of the network.
"""
function Result(nd::NetworkData)::Result
    boundary::Matrix{Int} = nd.rd.net_idx.networkboundary.matrix

    zs::Vector{Int} = zeros(Int, length(nd.wd.abundance))
    ns::Vector{Int} = zeros(Int, length(nd.wd.abundance))

    index::Int = 1
    for (z::Int, n_low::Int, n_high::Int) in eachrow(boundary)
        dn::Int = n_high - n_low
        zs[index:index+dn] = repeat([z], dn+1)
        ns[index:index+dn] = range(n_low, n_high, length=dn+1)
        index += dn + 1
    end
    return Result(zs, ns, nd.wd.abundance)
end

"""
    dump_y(nd::NetworkData)::Nothing

Outputs the abundances of each species to the disk.
"""
function dump_y(nd::NetworkData)::Nothing
    result::Result = Result(nd)
    open(nd.rd.output_info.final_y_path, "w") do out_file
        # TODO: This should really write [Int Int Float64], but right now it gets promoted to [Float64 Float64 Float64]
        idxs::BitVector = .!iszero.(result.abundance)
        Zs::Vector{Int} = result.proton_nums[idxs]
        As::Vector{Int} = Zs .+ result.neutron_nums[idxs]
        ys::Vector{Float64} = result.abundance[idxs]
        writedlm(out_file, [Zs As ys])
    end
    return nothing
end

"""
    dump_y(nd::NetworkData)::Nothing

Outputs the abundances of each mass number to the disk.
"""
function dump_ya(nd::NetworkData)::Nothing
    result::Result = Result(nd)
    open(nd.rd.output_info.final_ya_path, "w") do out_file
        # TODO: This should really write [Int Int Float64], but right now it gets promoted to [Float64 Float64 Float64]
        idxs::BitVector = .!iszero.(result.abundance)
        Zs::Vector{Int} = result.proton_nums[idxs]
        As::Vector{Int} = Zs .+ result.neutron_nums[idxs]
        ys::Vector{Float64} = result.abundance[idxs]
        As_unique::SortedSet{Int} = SortedSet(As)
        yas::Vector{Float64} = zeros(Float64,length(As_unique))
        for (i, a) in zip(1:length(As_unique), As_unique)
            yas[i] = sum(ys[As .== a])
        end
        writedlm(out_file, [collect(As_unique) yas])
    end
    return nothing
end

"""
    dump_result(nd::NetworkData)::Nothing

Outputs the results of the network calculations to the disk. This is usually called at the
end of network calculation.
"""
function dump_result(nd::NetworkData)::Nothing
    if nd.rd.output_info.dump_final_y
        dump_y(nd)
    end
    if nd.rd.output_info.dump_final_ya
        dump_ya(nd)
    end
    return nothing
end

"""
    dump_iteration(nd::NetworkData)::Nothing

Outputs the current state of the network calculations to the disk. This is usually called
after any iteration.
"""
function dump_iteration(nd::NetworkData, iteration::Int)::Nothing
    if !nd.rd.output_info.dump_each_iteration
        return
    end

    result::Result = Result(nd)
    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)
    mode::String = iteration == 0 ? "w" : "a"
    open(nd.rd.output_info.iteration_output_path, mode) do out_file
        write(out_file, "$(iteration)\t$(nd.wd.time.current)\t$(curr_traj.temperature)\t$(curr_traj.density)\n")
        idxs::BitVector = result.abundance .> 1e-15 #.!iszero.(result.abundance)
        Zs::Vector{Int} = result.proton_nums[idxs]
        As::Vector{Int} = Zs .+ result.neutron_nums[idxs]
        ys::Vector{Float64} = result.abundance[idxs]
        writedlm(out_file, [Zs As ys])
        write(out_file, "\n")
    end
    return nothing
end

function save_checkpoint(nd::NetworkData)::Nothing
    if !ismissing(nd.rd.output_info.checkpoint.path)
        for (i, (time, done)) in enumerate(zip(nd.rd.output_info.checkpoint.times, nd.rd.output_info.checkpoint.completed))
            if done || nd.wd.time.current < time
                continue
            end
            nd.rd.output_info.checkpoint.completed[i] = true
            println("Saving checkpoint at time $(nd.wd.time.current)")
            save_object("$(nd.rd.output_info.checkpoint.path)_$(time).jld", nd)
        end
    end
end

function load_checkpoint(path::String)::NetworkData
    return load_object(path)
end

function read_boundary(path::String)::Nothing
    # Read Fortran formatted extent file. It determines the limit on the nuclear chart.
    raw_boundary::Matrix{Int} = readdlm(path, ' ', Int)
    fill_boundary(raw_boundary)
end

function read_ncap!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    ncaps::Vector{NeutronCaptureIO} = load_object(path)
    num_not_in_network::Int = 0
    for ncap::NeutronCaptureIO in ncaps
        # Make sure every species involved in the reaction is in the network
        out_of_network::Bool = false
        for species::Tuple{Int, Int} in [ncap.product; ncap.reactant]
            z::Int, n::Int = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            num_not_in_network += 1
            continue
        end

        # Add the reaction to the dictionary
        z_r::Int, n_r::Int = ncap.reactant[2]
        reactant_idx::Int = zn_to_index(z_r, n_r, net_idx)
        reaction_data_io.ncap_dict[reactant_idx] = ncap
    end
    if num_not_in_network > 0
        println("$(num_not_in_network) reactions outside the network in file $(path)")
    end
end

function read_probdecay!(reaction_data_io::ReactionDataIO, path::String, id::String, net_idx::NetworkIndex)::Nothing
    probdecay::Vector{ProbDecayIO} = load_object(path)
    num_not_in_network::Int = 0
    for decay::ProbDecayIO in probdecay
        # Make sure every species involved in the reaction is in the network
        out_of_network::Bool = false
        for species::Tuple{Int, Int} in [decay.product; decay.reactant]
            z::Int, n::Int = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            num_not_in_network += 1
            continue
        end

        # Check if we already have this reaction in the network
        # idx::Union{Nothing, Int} = findfirst(other -> check_eq_reaction(decay, other[2]), reaction_data_io.probdecay)
        idx::Union{Nothing, Int} = findfirst(other -> decay.reactant == other[2].reactant && id == other[1], reaction_data_io.probdecay)
        if !isnothing(idx)
            # Replace the old data
            reaction_data_io.probdecay[idx] = (id, decay)
        else
            # Add the reaction to the array
            push!(reaction_data_io.probdecay, (id, decay))
        end
    end
    if num_not_in_network > 0
        println("$(num_not_in_network) reactions outside the network in file $(path)")
    end
end

function read_alphadecay!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    alphadecay::Vector{AlphaDecayIO} = load_object(path)
    num_not_in_network::Int = 0
    for decay::AlphaDecayIO in alphadecay
        # Make sure every species involved in the reaction is in the network
        out_of_network::Bool = false
        for species::Tuple{Int, Int} in [decay.product; decay.reactant]
            z::Int, n::Int = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            num_not_in_network += 1
            continue
        end

        # Check if we already have this reaction in the network
        # idx::Union{Nothing, Int} = findfirst(other -> check_eq_reaction(decay, other), reaction_data_io.alphadecay)
        idx::Union{Nothing, Int} = findfirst(other -> decay.reactant == other.reactant, reaction_data_io.alphadecay)
        if !isnothing(idx) # It already was in the network
            # Replace the old data
            reaction_data_io.alphadecay[idx] = decay
        else # It was not in the network
            # Add the reaction to the array
            push!(reaction_data_io.alphadecay, decay)
        end
    end
    if num_not_in_network > 0
        println("$(num_not_in_network) reactions outside the network in file $(path)")
    end
end

function read_photodissociation!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    photodissociation_dict::Dict{Tuple{Int, Int}, Photodissociation} = load_object(path)
    for (reactant::Tuple{Int, Int}, photodissociation::Photodissociation) in photodissociation_dict
        z_r::Int, n_r::Int = reactant
        # Make sure the product is in the network
        if !zn_in_network(z_r, n_r-1, net_idx)
            continue
        end
        product_idx::Int = zn_to_index(z_r, n_r-1, net_idx)

        # Make sure we have the forward rate associated with this reverse rate
        if !haskey(reaction_data_io.ncap_dict, product_idx)
            continue
        end
        ncap::NeutronCaptureIO = reaction_data_io.ncap_dict[product_idx]

        # Add the q value to the neutroncapture
        ncap.q = photodissociation.q
    end
end

function read_probrxn!(reaction_data_io::ReactionDataIO, path::String, id::String, net_idx::NetworkIndex)::Nothing
    probrxns::Vector{ProbRxnIO} = load_object(path)
    num_not_in_network::Int = 0
    for rxn::ProbRxnIO in probrxns
        # Make sure every species involved in the reaction is in the network
        out_of_network::Bool = false
        for (z::Int, n::Int) in [rxn.products; rxn.reactants]
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            num_not_in_network += 1
            continue
        end

        # Check if we already have this reaction in the network
        idx::Union{Nothing, Int} = findfirst(other -> rxn.reactants == other[2].reactants && id == other[1], reaction_data_io.probrxn)
        if !isnothing(idx)
            # Replace the old data
            reaction_data_io.probrxn[idx] = (id, rxn)
        else
            # Add the reaction to the array
            push!(reaction_data_io.probrxn, (id, rxn))
        end
    end
    if num_not_in_network > 0
        println("$(num_not_in_network) reactions outside the network in file $(path)")
    end
end

function read_rxn!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    rxns::Vector{RxnIO} = load_object(path)
    num_not_in_network::Int = 0
    for rxn::RxnIO in rxns
        # Make sure every species involved in the reaction is in the network
        out_of_network::Bool = false
        for (z::Int, n::Int) in [rxn.products; rxn.reactants]
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            num_not_in_network += 1
            continue
        end

        # Check if we already have this reaction in the network
        idx::Union{Nothing, Int} = findfirst(other -> rxn.reactants == other.reactants && rxn.products == other.products, reaction_data_io.rxn)
        if !isnothing(idx)
            # Replace the old data
            reaction_data_io.rxn[idx] = rxn
        else
            # Add the reaction to the array
            push!(reaction_data_io.rxn, rxn)
        end
    end
    if num_not_in_network > 0
        println("$(num_not_in_network) reactions outside the network in file $(path)")
    end
end

function read_decay!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    decays::Vector{DecayIO} = load_object(path)
    num_not_in_network::Int = 0
    for decay::DecayIO in decays
        # Make sure every species involved in the reaction is in the network
        out_of_network::Bool = false
        for (z::Int, n::Int) in [decay.products; decay.reactants]
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            num_not_in_network += 1
            continue
        end

        # Check if we already have this reaction in the network
        idx::Union{Nothing, Int} = findfirst(other -> decay.reactants == other.reactants && decay.products == other.products, reaction_data_io.decay)
        if !isnothing(idx)
            # Replace the old data
            reaction_data_io.decay[idx] = decay
        else
            # Add the reaction to the array
            push!(reaction_data_io.decay, decay)
        end
    end
    if num_not_in_network > 0
        println("$(num_not_in_network) reactions outside the network in file $(path)")
    end
end

function read_dataset!(reaction_data_io::ReactionDataIO, included_reactions::IncludedReactions, dataset::DataStructures.OrderedDict{String, Any}, net_idx::NetworkIndex)
    if !get(dataset, "active", false)
        return
    end
    id = get(dataset, "id", "")
    if dataset["rxn_type"] == "ncap"
        println("Reading neutron capture...")
        read_ncap!(reaction_data_io, dataset["path"], net_idx)
        included_reactions.ncap = true
    elseif dataset["rxn_type"] == "probdecay"
        println("Reading probdecay...")
        read_probdecay!(reaction_data_io, dataset["path"], id, net_idx)
        included_reactions.probdecay = true
    elseif dataset["rxn_type"] == "alphadecay"
        println("Reading alpha decay...")
        read_alphadecay!(reaction_data_io, dataset["path"], net_idx)
        included_reactions.alphadecay = true
    elseif dataset["rxn_type"] == "photodissociation"
        println("Reading photodissociation...")
        read_photodissociation!(reaction_data_io, dataset["path"], net_idx)
        included_reactions.photodissociation = true
    elseif dataset["rxn_type"] == "probrxn"
        println("Reading probreaction...")
        read_probrxn!(reaction_data_io, dataset["path"], id, net_idx)
        included_reactions.probrxn = true
    elseif dataset["rxn_type"] == "rxn"
        println("Reading rxn...")
        read_rxn!(reaction_data_io, dataset["path"], net_idx)
        included_reactions.rxn = true
    elseif dataset["rxn_type"] == "decay"
        println("Reading decay...")
        read_decay!(reaction_data_io, dataset["path"], net_idx)
        included_reactions.decay = true
    else
        error("Unknown reaction type: $(dataset["rxn_type"])")
    end
end

function read_initial_abundance(path::String, net_idx::NetworkIndex)::Vector{Float64}
    data::Vector{Tuple{Int, Int, Float64}} = load_object(path)
    total::Float64 = sum(zay->zay[3], data)
    data = [(z, a, y/total) for (z, a, y) in data]
    abundance::Vector{Float64} = zeros(Float64, get_networksize(net_idx))
    for (z::Int, a::Int, y::Float64) in data
        n::Int = a - z
        while !zn_in_network(z, n, net_idx)
            n -= 1
        end
        idx::Int = zn_to_index(z, n, net_idx)
        abundance[idx] += y/(z+n)
    end
    return abundance
end

function get_solver(type::String)
    if uppercase(type) == "UMFPACK"
        println("Using UMFPACK")
        return LS_UMFPACK()
    elseif uppercase(type) == "PARDISO"
        println("Using Pardiso")
        return LS_MKLPardisoSolver(MKLPardisoSolver())
    else
        println("Using UMFPACK")
        return LS_UMFPACK()
    end
end

function post_process_probdecay(probdecay::Vector{Tuple{String, ProbDecayIO}}, net_idx::NetworkIndex)
    probdecay_data = Vector{ProbDecay}()
    for (_, decay) in probdecay
        # Remove reactions with zero rate
        if iszero(decay.rate) || iszero(decay.average_number)
            continue
        end

        reactant_idxs = SVector{1, Int}(zn_to_index(decay.reactant[1][1], decay.reactant[1][2], net_idx))
        product_idxs = Vector{Int}([zn_to_index(product[1], product[2], net_idx) for product in decay.product])
        real_decay = ProbDecay(reactant_idxs, product_idxs, decay.average_number, decay.rate)
        push!(probdecay_data, real_decay)
    end
    return probdecay_data
end

function post_process_ncap(ncap_dict::Dict{Int, NeutronCaptureIO}, net_idx::NetworkIndex)
    if length(ncap_dict) == 0
        return []
    end
    vec_size = maximum(keys(ncap_dict))
    ncap_data = Vector{Union{Nothing, NeutronCapture}}(nothing, vec_size)
    for (idx, ncap) in ncap_dict
        if isnothing(ncap) || idx > vec_size
            continue
        end

        # Remove reactions that have zero rate
        if iszero(ncap.rates_pfuncs_lerp.rates)
            continue
        end

        reactant_idxs = SVector{2, Int}([zn_to_index(reactant[1], reactant[2], net_idx) for reactant in ncap.reactant])
        product_idxs = SVector{1, Int}([zn_to_index(product[1], product[2], net_idx) for product in ncap.product])
        A_r = ncap.product[1][1] + ncap.product[1][2]
        A_p1 = ncap.reactant[1][1] + ncap.reactant[1][2]
        A_p2 = ncap.reactant[2][1] + ncap.reactant[2][2]
        A_factor = (A_p1*A_p2/A_r)^(3/2)

        real_ncap = NeutronCapture(reactant_idxs, product_idxs, ncap.rates_pfuncs_lerp, ncap.q, A_factor)
        ncap_data[idx] = real_ncap
    end
    return ncap_data
end

function post_process_alpha(alphadecay::Vector{AlphaDecayIO}, net_idx::NetworkIndex)
    alpha_data = Vector{AlphaDecay}()
    for decay in alphadecay
        # Get rid of reactions with zero rate
        if iszero(decay.rate)
            continue
        end

        reactant_idx = zn_to_index(decay.reactant[1], decay.reactant[2], net_idx)
        product_idxs = SVector{2, Int}([zn_to_index(product[1], product[2], net_idx) for product in decay.product])
        real_decay = AlphaDecay(reactant_idx, product_idxs, decay.rate)
        push!(alpha_data, real_decay)
    end
    return alpha_data
end

function post_process_probrxn(probrxn::Vector{Tuple{String, ProbRxnIO}}, net_idx::NetworkIndex)
    probrxn_data = Vector{ProbRxn}()
    for (_, rxn) in probrxn
        # Remove reactions with zero rate
        if iszero(rxn.rates_lerp.rates) || iszero(rxn.average_numbers)
            continue
        end

        reactant_idxs = Vector{Int}([zn_to_index(reactant[1], reactant[2], net_idx) for reactant in rxn.reactants])
        product_idxs = Vector{Int}([zn_to_index(product[1], product[2], net_idx) for product in rxn.products])
        real_rxn = ProbRxn(reactant_idxs, product_idxs, rxn.average_numbers, rxn.rates_lerp)
        push!(probrxn_data, real_rxn)
    end
    return probrxn_data
end

function post_process_rxn(rxns::Vector{RxnIO}, net_idx::NetworkIndex)
    rxn_data = Vector{Rxn}()
    for rxn in rxns
        # Remove reactions with zero rate
        if iszero(rxn.rates_lerp.rates)
            continue
        end

        reactant_idxs = Vector{Int}([zn_to_index(reactant[1], reactant[2], net_idx) for reactant in rxn.reactants])
        product_idxs = Vector{Int}([zn_to_index(product[1], product[2], net_idx) for product in rxn.products])
        real_rxn = Rxn(reactant_idxs, product_idxs, rxn.rates_lerp)
        push!(rxn_data, real_rxn)
    end
    return rxn_data
end

function post_process_decay(decays::Vector{DecayIO}, net_idx::NetworkIndex)
    decay_data = Vector{Decay}()
    for decay in decays
        # Remove reactions with zero rate
        if iszero(decay.rate)
            continue
        end

        reactant_idxs = Vector{Int}([zn_to_index(reactant[1], reactant[2], net_idx) for reactant in decay.reactants])
        product_idxs = Vector{Int}([zn_to_index(product[1], product[2], net_idx) for product in decay.products])
        real_decay = Decay(reactant_idxs, product_idxs, decay.rate)
        push!(decay_data, real_decay)
    end
    return decay_data
end

function post_process_reaction_data(reaction_data_io::ReactionDataIO, net_idx::NetworkIndex)
    probdecay_data = post_process_probdecay(reaction_data_io.probdecay, net_idx)
    ncap_data = post_process_ncap(reaction_data_io.ncap_dict, net_idx)
    alpha_data = post_process_alpha(reaction_data_io.alphadecay, net_idx)
    probrxn_data = post_process_probrxn(reaction_data_io.probrxn, net_idx)
    rxn_data = post_process_rxn(reaction_data_io.rxn, net_idx)
    decay_data = post_process_decay(reaction_data_io.decay, net_idx)

    reaction_data::ReactionData = ReactionData(probdecay_data, ncap_data, alpha_data, probrxn_data, rxn_data, decay_data)
    return reaction_data
end

"""
    initialize_network_data(path::String)

Reads in the control file at `path` and constructs a [`NetworkData`](@ref) from the information provided in the control file.
"""
function initialize_network_data(path::String)
    # Parse the JSON control file
    println("Parsing JSON...")
    j = JSON.parsefile(path, dicttype=DataStructures.OrderedDict)

    if get(get(get(j, "network", Dict()), "checkpoint", Dict()), "use", false)
        println("Loading checkpoint $(j["network"]["checkpoint"]["path"])...")
        nd = load_checkpoint(j["network"]["checkpoint"]["path"])
        return nd
    end

    # Get the network index
    println("Reading index...")
    net_idx::NetworkIndex = load_object(j["network"]["extent"]["path"])

    # Get the reaction data
    reaction_data_io::ReactionDataIO = ReactionDataIO([], Dict(), [], [], [], [])
    included_reactions::IncludedReactions = IncludedReactions(false, false, false, false, false, false, false)
    for dataset in j["reactions"]
        read_dataset!(reaction_data_io, included_reactions, dataset, net_idx)
    end
    reaction_data = post_process_reaction_data(reaction_data_io, net_idx)

    # Get the trajectory
    println("Reading trajectory...")
    trajectory::TrajectoryLerp = load_object(j["conditions"]["trajectory"]["path"])

    # Get the initial abundance
    println("Reading initial abundance...")
    abundance::Vector{Float64} = read_initial_abundance(j["conditions"]["initial_composition"]["path"], net_idx)
    yproposed::Vector{Float64} = Vector{Float64}(undef, length(abundance))

    # The size of the network
    networksize::Int = get_networksize(net_idx)

    # Get the time
    println("Reading time data...")
    time::Time = Time(0.0, 1e-15, 0.0)
    if get(j["network"]["start"], "use_trajectory", false)
        time.current = trajectory.times[1]
    elseif haskey(j["network"]["start"], "time")
        time.current = j["network"]["start"]["time"]
    end
    if get(j["network"]["stop"], "use_trajectory", false)
        time.stop = trajectory.times[end]
    elseif haskey(j["network"]["stop"], "time")
        time.stop = j["network"]["stop"]["time"]
    end

    # Initialize ydot
    println("Creating Ẏ...")
    ydot::Vector{Float64} = zeros(Float64, networksize)

    # Create the Jacobian
    println("Creating the Jacobian...")
    jacobian::SparseMatrixCSC{Float64, Int} = spzeros(Float64, Int, (networksize, networksize))

    # Grab the output info
    println("Reading the output information...")
    dump_final_y::Bool = get(get(get(j, "output", Dict()), "y", Dict()), "active", false)
    final_y_path::Union{Missing, String} = dump_final_y ? j["output"]["y"]["path"] : missing
    dump_final_ya::Bool = get(get(get(j, "output", Dict()), "ya", Dict()), "active", false)
    final_ya_path::Union{Missing, String} = dump_final_y ? j["output"]["ya"]["path"] : missing
    dump_each_iteration::Bool = get(get(get(j, "output", Dict()), "ytime", Dict()), "active", false)
    iteration_output_path::Union{Missing, String} = dump_final_y ? j["output"]["ytime"]["path"] : missing

    # Checkpoints
    checkpoint_path::Union{Missing, String} = get(get(get(j, "output", Dict()), "checkpoints", Dict()), "path", missing)
    checkpoint_times::Vector{Float64} = Vector{Float64}()
    if !ismissing(checkpoint_path)
        checkpoint_times = j["output"]["checkpoints"]["times"]
    end
    checkpoint_completed = BitVector(repeat([false], length(checkpoint_times)))
    checkpoint = Checkpoint(checkpoint_path, checkpoint_times, checkpoint_completed)

    output_info::OutputInfo = OutputInfo(dump_final_y, final_y_path, dump_final_ya, final_ya_path, dump_each_iteration, iteration_output_path, checkpoint)

    # Get the linear solver
    println("Initializing solver...")
    solver = get_solver(j["computational"]["solver"]["type"])

    # Create the network data
    ro_data = ROData(net_idx, reaction_data, trajectory, output_info, included_reactions)
    rw_data = RWData(abundance, yproposed, ydot, time, jacobian, solver)
    nd::NetworkData = NetworkData(ro_data, rw_data)

    # Update both the Jacobian and ydot
    println("Initializing the network...")
    fill_jacobian!(nd)
    update_ydot!(nd)

    return nd
end

end