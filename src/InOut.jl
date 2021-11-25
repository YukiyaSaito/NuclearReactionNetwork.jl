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
    boundary::Matrix{Int} = nd.net_idx.networkboundary.matrix

    zs::Vector{Int} = zeros(Int, length(nd.abundance))
    ns::Vector{Int} = zeros(Int, length(nd.abundance))

    index::Int = 1
    for (z::Int, n_low::Int, n_high::Int) in eachrow(boundary)
        dn::Int = n_high - n_low
        zs[index:index+dn] = repeat([z], dn+1)
        ns[index:index+dn] = range(n_low, n_high, length=dn+1)
        index += dn + 1
    end
    return Result(zs, ns, nd.abundance)
end

"""
    dump_y(nd::NetworkData)::Nothing

Outputs the abundances of each species to the disk.
"""
function dump_y(nd::NetworkData)::Nothing
    result::Result = Result(nd)
    open(nd.output_info.final_y_path, "w") do out_file
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
    open(nd.output_info.final_ya_path, "w") do out_file
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
    if nd.output_info.dump_final_y
        dump_y(nd)
    end
    if nd.output_info.dump_final_ya
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
    if !nd.output_info.dump_each_iteration
        return
    end

    result::Result = Result(nd)
    curr_traj::CurrentTrajectory = get_current_trajectory(nd.trajectory, nd.time.current)
    mode::String = iteration == 0 ? "w" : "a"
    open(nd.output_info.iteration_output_path, mode) do out_file
        write(out_file, "$(iteration)\t$(nd.time.current)\t$(curr_traj.temperature)\t$(curr_traj.density)\n")
        idxs::BitVector = .!iszero.(result.abundance)
        Zs::Vecotr{Int} = result.proton_nums[idxs]
        As::Vector{Int} = Zs .+ result.neutron_nums[idxs]
        ys::Vector{Int} = result.abundance[idxs]
        writedlm(out_file, [Zs As ys])
        write(out_file, "\n")
    end
end

function read_boundary(path::String)::Nothing
    # Read Fortran formatted extent file. It determines the limit on the nuclear chart.
    raw_boundary::Matrix{Int} = readdlm(path, ' ', Int)
    fill_boundary(raw_boundary)
end

function read_ncap!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    ncaps::Vector{NeutronCaptureIO} = load_object(path)
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
            continue
        end

        # Add the reaction to the dictionary
        z_r::Int, n_r::Int = ncap.reactant[2]
        reactant_idx::Int = zn_to_index(z_r, n_r, net_idx)
        reaction_data_io.ncap_dict[reactant_idx] = ncap
    end
end

function read_probdecay!(reaction_data_io::ReactionDataIO, path::String, id::String, net_idx::NetworkIndex)::Nothing
    probdecay::Vector{ProbDecayIO} = load_object(path)
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
            continue
        end

        # Remove reactions with zero rate
        if iszero(decay.rate)
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
end

function read_alphadecay!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    alphadecay::Vector{AlphaDecayIO} = load_object(path)
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
            continue
        end

        if iszero(decay.rate)
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
end

function read_photodissociation!(reaction_data_io::ReactionDataIO, path::String, net_idx::NetworkIndex)::Nothing
    photodissociation_dict::Dict{Tuple{Int, Int}, Photodissociation} = load_object(path)
    for (reactant::Tuple{Int, Int}, photodissociation::Photodissociation) in photodissociation_dict
        z_r::Int, n_r::Int = reactant
        # Make sure the reactant is in the network
        if !zn_in_network(z_r, n_r-1, net_idx)
            continue
        end
        reactant_idx::Int = zn_to_index(z_r, n_r-1, net_idx)

        # Make sure we have the forward rate associated with this reverse rate
        if !haskey(reaction_data_io.ncap_dict, reactant_idx)
            continue
        end
        ncap::NeutronCaptureIO = reaction_data_io.ncap_dict[reactant_idx]

        # Add the q value to the neutroncapture
        ncap.q = photodissociation.q
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
        println("Reading beta decay...")
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
        reactant_idxs = SVector{1, Int}(zn_to_index(decay.reactant[1][1], decay.reactant[1][2], net_idx))
        product_idxs = Vector{Int}([zn_to_index(product[1], product[2], net_idx) for product in decay.product])
        real_decay = ProbDecay(reactant_idxs, product_idxs, decay.average_number, decay.rate)
        push!(probdecay_data, real_decay)
    end
    return probdecay_data
end

function post_process_ncap(ncap_dict::Dict{Int, NeutronCaptureIO}, net_idx::NetworkIndex)
    vec_size = maximum(keys(ncap_dict))
    ncap_data = Vector{Union{Nothing, NeutronCapture}}(nothing, vec_size)
    for (idx, ncap) in ncap_dict
        if isnothing(ncap)
            real_ncap = nothing
        else
            reactant_idxs = SVector{2, Int}([zn_to_index(reactant[1], reactant[2], net_idx) for reactant in ncap.reactant])
            product_idxs = SVector{1, Int}([zn_to_index(product[1], product[2], net_idx) for product in ncap.product])
            A_r = ncap.product[1][1] + ncap.product[1][2]
            A_p1 = ncap.reactant[1][1] + ncap.reactant[1][2]
            A_p2 = ncap.reactant[2][1] + ncap.reactant[2][2]
            A_factor = (A_p1*A_p2/A_r)^(3/2)
            real_ncap = NeutronCapture(reactant_idxs, product_idxs, ncap.rates_pfuncs_lerp, ncap.q, A_factor)
        end
        ncap_data[idx] = real_ncap
    end
    return ncap_data
end

function post_process_alpha(alphadecay::Vector{AlphaDecayIO}, net_idx::NetworkIndex)
    alpha_data = Vector{AlphaDecay}()
    for decay in alphadecay
        reactant_idx = zn_to_index(decay.reactant[1], decay.reactant[2], net_idx)
        product_idxs = SVector{2, Int}([zn_to_index(product[1], product[2], net_idx) for product in decay.product])
        real_decay = AlphaDecay(reactant_idx, product_idxs, decay.rate)
        push!(alpha_data, real_decay)
    end
    return alpha_data
end

function post_process_reaction_data(reaction_data_io::ReactionDataIO, net_idx::NetworkIndex)
    probdecay_data = post_process_probdecay(reaction_data_io.probdecay, net_idx)
    ncap_data = post_process_ncap(reaction_data_io.ncap_dict, net_idx)
    alpha_data = post_process_alpha(reaction_data_io.alphadecay, net_idx)

    reaction_data::ReactionData = ReactionData(probdecay_data, ncap_data, alpha_data)
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

    # Get the network index
    println("Reading index...")
    net_idx::NetworkIndex = load_object(j["network"]["extent"]["path"])

    # Get the reaction data
    reaction_data_io::ReactionDataIO = ReactionDataIO([], Dict(), [])
    included_reactions::IncludedReactions = IncludedReactions(false, false, false, false)
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
    output_info::OutputInfo = OutputInfo(dump_final_y, final_y_path, dump_final_ya, final_ya_path, dump_each_iteration, iteration_output_path)

    # Get the linear solver
    println("Initializing solver...")
    solver = get_solver(j["computational"]["solver"]["type"])

    # Create the network data
    nd::NetworkData = NetworkData(net_idx, reaction_data, trajectory, abundance, yproposed, ydot, time, jacobian, output_info, included_reactions, solver)

    # Update both the Jacobian and ydot
    println("Initializing the network...")
    fill_jacobian!(nd)
    update_ydot!(nd)

    return nd
end

end