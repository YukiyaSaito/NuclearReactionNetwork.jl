"""
    Network

Provides basic functionality for the network.
"""
module Network

export NetworkBoundary
export NetworkIndex
export get_networksize
export zn_to_index
export zn_in_network
export Time
export step_time!

"""
    Time

Keeps track of the current time of the simulation, when it should stop, and the step size.

# Fields
- `current::Float64`: The current time of the simulation.
- `step::Float64`: The current step size of the simulation.
- `stop::Float64`: The time the simulation should end on.
"""
mutable struct Time
    """The current time of the simulation."""
    current::Float64
    """The current step size of the simulation."""
    step::Float64
    """The time the simulation should end on."""
    stop::Float64
end

"""
    step_time!(time::Time)

Advance the current time by `time.step`.
"""
function step_time!(time::Time)::Union{Nothing, Float64}
    time.current += time.step

    # Cap the time
    if time.current >= time.stop
        # time.step = time.stop - (time.current - time.step) TODO: Add this?
        time.current = time.stop
    end
end

"""
    NetworkBoundary{T<:Matrix{Int}}

Keeps track of what nuclear species to include in the network calculations.

# Fields:
- `matrix::T`: A matrix that tells us which species are in the network.
"""
struct NetworkBoundary{T<:Matrix{Int}}
    matrix::T
end

# TODO: Do we need a better name for this?
# TODO: Does it make sense to hold the mass vector with these other fields?
"""
   NetworkIndex
   
Useful data structure to use to check if a species is in the network.

# Fields:
- `networkboundary::NetworkBoundary{Matrix{Int}}`: Tells us which species are in the network.
- `cum_isotopes::Vector{Int}`: A helper piece of data to quickly identify if a species is in the network.
- `mass_vector::Vector{Float64}:` A vector holding the atomic masses of all the species in the network
"""
struct NetworkIndex
    networkboundary::NetworkBoundary{Matrix{Int}}
    cum_isotopes::Vector{Int}
    mass_vector::Vector{Float64}
end

"""
    get_networksize(net_idx::NetworkIndex)

Return the number of species in the network.
"""
function get_networksize(net_idx::NetworkIndex)::Int
    return net_idx.cum_isotopes[end]
end

"""
    get_networksize(networkboundary::NetworkBoundary)

Return the number of species in the network.
"""
function get_networksize(networkboundary::NetworkBoundary)::Int
    boundary::Matrix{Int} = networkboundary.matrix
    networksize::Int = 0
    for z::Int in 1:size(boundary, 1)
        n_low::Int = boundary[z,2]
        n_high::Int = boundary[z,3]
        networksize::Int += n_high - n_low + 1
    end
    return networksize
    #get boundary shape and loop over the first (second) dimension. Make the inner loop over column.
end

"""
    NetworkIndex(networkboundary::NetworkBoundary)

Construct a `NetworkIndex` from a `NetworkBoundary`.
"""
function NetworkIndex(networkboundary::NetworkBoundary)::NetworkIndex
    boundary::Matrix{Int} = networkboundary.matrix
    num_isotopes::Vector{Int} = @views boundary[:, 3] - boundary[:, 2] .+ 1
    cum_isotopes::Vector{Int} = cumsum(num_isotopes)

    networksize::Int = get_networksize(networkboundary)
    mass_vector::Vector{Float64} = zeros(Float64, networksize)
    index::Int = 1
    for (z::Int, n_low::Int, n_high::Int) in eachrow(boundary)
        dn::Int = n_high - n_low
        mass_vector[index:index+dn] = range(n_low, n_high, length=dn+1) .+ z
        index += dn + 1
    end

    return NetworkIndex(networkboundary, cum_isotopes, mass_vector)
end

"""
    zn_to_index(z::Int, n::Int, net_idx::NetworkIndex)

Convert the pair (`z`, `n`) to an index that can be used to index into other arrays.
"""
function zn_to_index(z::Int, n::Int, net_idx::NetworkIndex)::Int
    idx_so_far::Int = z > 0 ? net_idx.cum_isotopes[z] : 0
    n_low::Int = net_idx.networkboundary.matrix[z+1, 2]
    return idx_so_far + n - n_low + 1
end

"""
    zn_in_network(z::Int, n::Int, net_idx::NetworkIndex)

Check if the pair (`z`, `n`) is in the network.
"""
function zn_in_network(z::Int, n::Int, net_idx::NetworkIndex)::Bool
    # TODO: Maybe only accept UInts and not Ints?
    if z+1 > size(net_idx.networkboundary.matrix, 1) || z < 0
        return false
    end
    return n <= net_idx.networkboundary.matrix[z+1, 3] 
end

end