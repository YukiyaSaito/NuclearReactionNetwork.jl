module Network
# using ..InOut
using DelimitedFiles

export NetworkBoundary
export read_boundary
export get_networksize
export zn_to_index

struct NetworkBoundary{T<:Matrix{Int64}}
    matrix::T
end

# to do: organize zn_to_index_dict and mass_vector into struct

function read_boundary(path::String)
    raw_boundary::Matrix{Int64} = readdlm(path,' ',Int)
    fill_boundary(raw_boundary)
end

function fill_boundary(raw_boundary::Matrix{Int64}) # Function berriering
    return NetworkBoundary(raw_boundary)
end

function get_networksize(networkboundary::NetworkBoundary)
    boundary = networkboundary.matrix
    networksize::Int64 = 0
    # println(size(boundary,1))
    for z in 1:size(boundary,1)
        n_low::Int64 = boundary[z,2]
        n_high::Int64 = boundary[z,3]
        networksize += n_high - n_low + 1
    end
    return networksize
    #get boundary shape and loop over the first (second) dimension. Make the inner loop over column.
end

function zn_to_index(networkboundary::NetworkBoundary)
    boundary::Matrix{Int64} = networkboundary.matrix
    zn_to_index_dict = Dict{Vector{Int64},Int64}()
    index::Int64 = 1
    networksize::Int64 = get_networksize(networkboundary)
    mass_vector::Vector{Float64} = zeros(Float64,networksize)
    for nline in 1:size(boundary,1)
        z::Int64 = boundary[nline,1]
        n_low::Int64 = boundary[nline,2]
        n_high::Int64 = boundary[nline,3]
        for n in n_low:n_high
            zn_to_index_dict[[z,n]] = index
            mass_vector[index] = z+n
            index += 1            
        end        
    end
    return zn_to_index_dict, mass_vector
end


end