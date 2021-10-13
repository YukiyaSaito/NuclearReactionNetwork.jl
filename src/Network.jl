module Network

export NetworkBoundary
export NetworkIndex
export get_networksize
export zn_to_index
export zn_in_network
export Time
export step_time!

mutable struct Time
    current::Float64
    step::Float64
    stop::Float64
end

function step_time!(time::Time)::Union{Nothing, Float64}
    time.current += time.step

    # Cap the time
    if time.current >= time.stop
        # time.step = time.stop - (time.current - time.step) TODO: Add this?
        time.current = time.stop
    end
end

struct NetworkBoundary{T<:Matrix{Int}}
    matrix::T
end

# TODO: Do we need a better name for this?
# TODO: Does it make sense to hold the mass vector with these other fields?
struct NetworkIndex
    networkboundary::NetworkBoundary{Matrix{Int}}
    cum_isotopes::Vector{Int}
    mass_vector::Vector{Float64}
end

function get_networksize(net_idx::NetworkIndex)::Int
    return net_idx.cum_isotopes[end]
end

function get_networksize(networkboundary::NetworkBoundary)::Int # Get the number of species included in the network
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

function zn_to_index(z::Int, n::Int, net_idx::NetworkIndex)::Int
    idx_so_far::Int = z > 0 ? net_idx.cum_isotopes[z] : 0
    n_low::Int = net_idx.networkboundary.matrix[z+1, 2]
    return idx_so_far + n - n_low + 1
end

function zn_in_network(z::Int, n::Int, net_idx::NetworkIndex)::Bool
    # TODO: Maybe only accept UInts and not Ints?
    if z+1 > size(net_idx.networkboundary.matrix, 1) || z < 0
        return false
    end
    return n <= net_idx.networkboundary.matrix[z+1, 3] 
end

end