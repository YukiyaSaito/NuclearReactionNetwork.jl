module InOut

using DelimitedFiles
using DataStructures
using JSON
using JLD2
using SparseArrays
using Interpolations
using Pardiso
using ..Astro
using ..ReactionTypes
using ..Network
using ..NetworkDatas

export Result
export dump_result
export dump_iteration
export initialize_network_data

# TODO: Should this be a matrix rather than a series of vectors of the same size (we would lose type inhomogeneity)?
struct Result
    proton_nums::Vector{Int}
    neutron_nums::Vector{Int}
    abundance::Vector{Float64} # TODO: Add more things that we can output, like PRISM does
end

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

function dump_result(nd::NetworkData)::Nothing
    if nd.output_info.dump_final_y
        dump_y(nd)
    end
    if nd.output_info.dump_final_ya
        dump_ya(nd)
    end
    return nothing
end

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

function read_ncap!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)::Nothing
    ncaps::Vector{NeutronCapture} = load_object(path)
    for ncap::NeutronCapture in ncaps
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
        z_p::Int, n_p::Int = ncap.product[1]
        product_idx::Int = zn_to_index(z_p, n_p, net_idx)
        reaction_data.neutroncapture[product_idx] = ncap
    end
end

function read_probdecay!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)::Nothing
    probdecay::Vector{ProbDecay} = load_object(path)
    for decay::ProbDecay in probdecay
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

        # Check if we already have this reaction in the network
        idx::Union{Nothing, Int} = findfirst(other -> check_eq_reaction(decay, other), reaction_data.probdecay)
        if !isnothing(idx)
            # Replace the old data
            reaction_data.probdecay[idx] = decay
        else
            # Add the reaction to the array
            push!(reaction_data.probdecay, decay)
        end
    end
end

function read_alphadecay!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)::Nothing
    alphadecay::Vector{AlphaDecay} = load_object(path)
    for decay::AlphaDecay in alphadecay
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

        # Check if we already have this reaction in the network
        idx::Union{Nothing, Int} = findfirst(other -> check_eq_reaction(decay, other), reaction_data.alphadecay)
        if !isnothing(idx)
            # Replace the old data
            reaction_data.prodecay[idx] = decay
        else
            # Add the reaction to the array
            push!(reaction_data.alphadecay, decay)
        end
    end
end

function read_photodissociation!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)::Nothing
    photodissociation_dict::Dict{Tuple{Int, Int}, Photodissociation} = load_object(path)
    for (reactant::Tuple{Int, Int}, photodissociation::Photodissociation) in photodissociation_dict
        z_r::Int, n_r::Int = reactant
        # Make sure the reactant is in the network
        if !zn_in_network(z_r, n_r, net_idx)
            continue
        end
        reactant_idx::Int = zn_to_index(z_r, n_r, net_idx)

        # Make sure we have the forward rate associated with this reverse rate
        if !haskey(reaction_data.neutroncapture, reactant_idx)
            continue
        end
        ncap::NeutronCapture = reaction_data.neutroncapture[reactant_idx]

        # Add the q value to the neutroncapture
        ncap.q = photodissociation.q
    end
end

function read_dataset!(reaction_data::ReactionData, included_reactions::IncludedReactions, dataset::DataStructures.OrderedDict{String, Any}, net_idx::NetworkIndex)
    if !get(dataset, "active", false)
        return
    end
    if dataset["rxn_type"] == "ncap"
        println("Reading neutron capture...")
        read_ncap!(reaction_data, dataset["path"], net_idx)
        included_reactions.ncap = true
    elseif dataset["rxn_type"] == "probdecay"
        println("Reading beta decay...")
        read_probdecay!(reaction_data, dataset["path"], net_idx)
        included_reactions.probdecay = true
    elseif dataset["rxn_type"] == "alphadecay"
        println("Reading alpha decay...")
        read_alphadecay!(reaction_data, dataset["path"], net_idx)
        included_reactions.alphadecay = true
    elseif dataset["rxn_type"] == "photodissociation"
        println("Reading photodissociation...")
        read_photodissociation!(reaction_data, dataset["path"], net_idx)
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

function initialize_network_data(path::String)
    # Parse the JSON control file
    println("Parsing JSON...")
    j = JSON.parsefile(path, dicttype=DataStructures.OrderedDict)

    # Get the network index
    println("Reading index...")
    net_idx::NetworkIndex = load_object(j["network"]["extent"]["path"])

    # Get the reaction data
    reaction_data::ReactionData = initialize_reactions()
    included_reactions::IncludedReactions = IncludedReactions(false, false, false, false)
    for dataset in j["reactions"]
        read_dataset!(reaction_data, included_reactions, dataset, net_idx)
    end

    # Get the trajectory
    println("Reading trajectory...")
    trajectory::Trajectory = load_object(j["conditions"]["trajectory"]["path"])

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
        lerped_times = knots(trajectory.temperatures)
        time.current = lerped_times[1]
    elseif haskey(j["network"]["start"], "time")
        time.current = j["network"]["start"]["time"]
    end
    if get(j["network"]["stop"], "use_trajectory", false)
        lerped_times = knots(trajectory.temperatures)
        time.stop = lerped_times[length(lerped_times)]
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

    println("Initializing MKLPARDISO...")
    pardiso = MKLPardisoSolver()

    # Create the network data
    nd::NetworkData = NetworkData(net_idx, reaction_data, trajectory, abundance, yproposed, ydot, time, jacobian, output_info, included_reactions, pardiso)

    # Update both the Jacobian and ydot
    println("Initializing the network...")
    fill_jacobian!(nd)
    update_ydot!(nd)

    return nd
end

end