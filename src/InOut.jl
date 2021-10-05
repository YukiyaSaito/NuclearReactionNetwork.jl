module InOut

using DelimitedFiles
using DataStructures
using JSON
using JLD2
using SparseArrays
using Interpolations
using ..Astro
using ..ReactionTypes
using ..Network
using ..NetworkDatas

export Result
export NetworkData
export dump_result
export dump_iteration

# TODO: Should this be a matrix rather than a series of vectors of the same size (we would lose type inhomogeneity)?
struct Result
    proton_nums::Vector{Int64}
    neutron_nums::Vector{Int64}
    abundance::Vector{Float64} # TODO: Add more things that we can output, like PRISM does
    function Result(proton_nums::Vector{Int64}, neutron_nums::Vector{Int64}, abundance::Vector{Float64})
        return new(proton_nums, neutron_nums, abundance)
    end
end

function Result(nd::NetworkData)
    boundary = nd.net_idx.networkboundary.matrix

    zs = zeros(Int64, length(nd.abundance))
    ns = zeros(Int64, length(nd.abundance))

    index = 1
    for (z, n_low, n_high) in eachrow(boundary)
        dn = n_high - n_low
        zs[index:index+dn] = repeat([z], dn+1)
        ns[index:index+dn] = range(n_low, n_high, length=dn+1)
        index += dn + 1
    end
    return Result(zs, ns, nd.abundance)
end

function dump_result(nd::NetworkData)
    if !nd.output_info.dump_final_output
        return
    end

    result = Result(nd)
    open(nd.output_info.final_output_path, "w") do out_file
        # TODO: This should really write [Int64 Int64 Float64], but right now it gets promoted to [Float64 Float64 Float64]
        idxs = .!iszero.(result.abundance)
        Zs = result.proton_nums[idxs]
        As = Zs .+ result.neutron_nums[idxs]
        ys = result.abundance[idxs]
        writedlm(out_file, [Zs As ys])
    end
end

function dump_iteration(nd::NetworkData, iteration::Int64)
    if !nd.output_info.dump_each_iteration
        return
    end

    result = Result(nd)
    curr_traj = get_current_trajectory(nd.trajectory, nd.time.current)
    mode = iteration == 0 ? "w" : "a"
    open(nd.output_info.iteration_output_path, mode) do out_file
        write(out_file, "$(iteration)\t$(nd.time.current)\t$(curr_traj.temperature)\t$(curr_traj.density)\n")
        idxs = .!iszero.(result.abundance)
        Zs = result.proton_nums[idxs]
        As = Zs .+ result.neutron_nums[idxs]
        ys = result.abundance[idxs]
        writedlm(out_file, [Zs As ys])
        write(out_file, "\n")
    end
end

function read_boundary(path::String)
    # Read Fortran formatted extent file. It determines the limit on the nuclear chart.
    raw_boundary::Matrix{Int64} = readdlm(path, ' ', Int)
    fill_boundary(raw_boundary)
end

function read_ncap!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)
    ncap_dict::Dict{Vector{Vector{Int64}},NeutronCapture} = load_object(path)
    for ncap in values(ncap_dict)
        # Make sure every species involved in the reaction is in the network
        out_of_network = false
        for species in [ncap.product; ncap.reactant]
            z, n = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            continue
        end

        # Add the reaction to the dictionary
        reaction_data.neutroncapture[ncap.reactant] = ncap
    end
end

function read_probdecay!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)
    probdecay_dict::Dict{Vector{Vector{Int64}}, ProbDecay} = load_object(path)
    for decay in values(probdecay_dict)
        # Make sure every species involved in the reaction is in the network
        out_of_network = false
        for species in [decay.product; decay.reactant]
            z, n = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            continue
        end

        # Add the reaction to the dictionary
        reaction_data.probdecay[decay.reactant] = decay
    end
end

function read_alphadecay!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)
    alphadecay::Vector{AlphaDecay} = load_object(path)
    for decay in alphadecay
        # Make sure every species involved in the reaction is in the network
        out_of_network = false
        for species in [decay.product; decay.reactant]
            z, n = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            continue
        end

        # Add the reaction to the vector
        push!(reaction_data.alphadecay, decay)
    end
end

function read_photodissociation!(reaction_data::ReactionData, path::String, net_idx::NetworkIndex)
    photodissociation_dict::Dict{Vector{Vector{Int64}}, Photodissociation} = load_object(path)
    for (products, photodissociation) in photodissociation_dict
        key = products
        # Make sure we have the forward rate associated with this reverse rate
        if !haskey(reaction_data.neutroncapture, key)
            reverse!(key) # Maybe the key is backwards?
            if !haskey(reaction_data.neutroncapture, key)
                continue
            end
        end
        ncap = reaction_data.neutroncapture[key]

        # Make sure every species involved in the reaction is in the network
        out_of_network = false
        for species in [ncap.product; ncap.reactant]
            z, n = species
            if !zn_in_network(z, n, net_idx)
                out_of_network = true
                break
            end
        end
        if out_of_network
            continue
        end

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

function read_initial_abundance(path::String, net_idx::NetworkIndex)
    data = load_object(path)
    total = sum(zay->zay[3], data)
    data = [(z, a, y/total) for (z, a, y) in data]
    abundance = zeros(Float64, get_networksize(net_idx))
    for (z, a, y) in data
        n = a - z
        if zn_in_network(z, n, net_idx)
            idx = zn_to_index(z, n, net_idx)
            abundance[idx] += y/a
        else
            subtract = 1
            while !zn_in_network(z, a-subtract, net_idx) # FIXME: Are we sure this is correct?
                subtract += 1
            end
            idx = zn_to_index(z, n, net_idx)
            abundance[idx] += y/(z+n)
        end
    end
    return abundance
end

function NetworkData(path::String)
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
    networksize::Int64 = get_networksize(net_idx)

    # Get the time
    println("Reading time data...")
    time = Time()
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
    jacobian::Matrix{Float64} = zeros(Float64, (networksize, networksize))

    # Grab the output info
    println("Reading the output information...")
    dump_final_output::Bool = get(get(get(j, "output", Dict()), "y", Dict()), "active", false)
    final_output_path = dump_final_output ? j["output"]["y"]["path"] : missing
    dump_each_iteration::Bool = get(get(get(j, "output", Dict()), "ytime", Dict()), "active", false)
    iteration_output_path = dump_final_output ? j["output"]["ytime"]["path"] : missing
    output_info = OutputInfo(dump_final_output, final_output_path, dump_each_iteration, iteration_output_path)

    # Create the network data
    nd::NetworkData = NetworkData(net_idx, reaction_data, trajectory, abundance, yproposed, ydot, time, jacobian, output_info, included_reactions)

    println("Initializing the network...")
    # Update both the Jacobian and ydot
    fill_jacobian!(nd)
    other_thing(nd)
    # update_ydot!(nd)

    return nd
end

end