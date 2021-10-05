module Network
# using ..InOut
using DelimitedFiles

export NetworkBoundary
export NetworkIndex
export read_boundary
export get_networksize
export zn_to_index
export zn_in_network
export Time
export step_time!

mutable struct Time
    current::Float64
    step::Float64
    stop::Float64
    function Time(curr_time::Float64, time_step::Float64, stop_time::Float64)
        return new(curr_time, time_step, stop_time)
    end
    function Time()
        return new(0.0, 1e-15, 20.0)
    end
end

function step_time!(time::Time)
    time.current += time.step

    # Cap the time
    if time.current >= time.stop
        # time.step = time.stop - (time.current - time.step) TODO: Add this?
        time.current = time.stop
    end
end

struct NetworkBoundary{T<:Matrix{Int64}}
    matrix::T
end

# TODO: Do we need a better name for this?
# TODO: Does it make sense to hold the mass vector with these other fields?
struct NetworkIndex
    networkboundary::NetworkBoundary{Matrix{Int64}}
    cum_isotopes::Vector{Int64}
    mass_vector::Vector{Float64}
    function NetworkIndex(networkboundary::NetworkBoundary{Matrix{Int64}}, cum_isotopes::Vector{Int64}, mass_vector::Vector{Float64})
        return new(networkboundary, cum_isotopes, mass_vector)
    end
end

function read_boundary(path::String)
    # Read Fortran formatted extent file. It determines the limit on the nuclear chart.
    raw_boundary::Matrix{Int64} = readdlm(path,' ',Int)
    fill_boundary(raw_boundary)
end

function fill_boundary(raw_boundary::Matrix{Int64}) # Function berriering
    return NetworkBoundary(raw_boundary)
end

function get_networksize(net_idx::NetworkIndex)
    return net_idx.cum_isotopes[end]
end

function get_networksize(networkboundary::NetworkBoundary) # Get the number of species included in the network
    boundary = networkboundary.matrix
    networksize::Int64 = 0
    for z in 1:size(boundary,1)
        n_low::Int64 = boundary[z,2]
        n_high::Int64 = boundary[z,3]
        networksize += n_high - n_low + 1
    end
    return networksize
    #get boundary shape and loop over the first (second) dimension. Make the inner loop over column.
end

function NetworkIndex(networkboundary::NetworkBoundary)
    boundary = networkboundary.matrix
    num_isotopes = @views boundary[:, 3] - boundary[:, 2] .+ 1
    cum_isotopes = cumsum(num_isotopes)

    networksize = get_networksize(networkboundary)
    mass_vector = zeros(Float64, networksize)
    index = 1
    for (z, n_low, n_high) in eachrow(boundary)
        dn = n_high - n_low
        mass_vector[index:index+dn] = range(n_low, n_high, length=dn+1) .+ z
        index += dn + 1
    end

    return NetworkIndex(networkboundary, cum_isotopes, mass_vector)
end

function zn_to_index(z::Int64, n::Int64, net_idx::NetworkIndex)
    idx_so_far = z > 0 ? net_idx.cum_isotopes[z] : 0
    n_low = net_idx.networkboundary.matrix[z+1, 2]
    return idx_so_far + n - n_low + 1
end

function zn_in_network(z::Int64, n::Int64, net_idx::NetworkIndex)
    # TODO: Maybe only accept uints and not ints?
    if z+1 > size(net_idx.networkboundary.matrix, 1) || z < 0
        return false
    end
    return n <= net_idx.networkboundary.matrix[z+1, 3] 
end

end