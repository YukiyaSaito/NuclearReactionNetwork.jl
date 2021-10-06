module ReactionTypes

using Interpolations
using ..Astro
using ..Network

export ProbDecay
export NeutronCapture
export AlphaDecay
export Photodissociation
export ReactionData
export initialize_reactions

abstract type AbstractReaction end

# mutable struct NBodyReaction <: AbstractReaction
#     reactant::Vector{Vector{Int32}}
#     product::Vector{Vector{Int32}}
# end

mutable struct ProbDecay <: AbstractReaction
    reactant::Vector{Tuple{Int64, Int64}} # [(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]
    product::Vector{Tuple{Int64, Int64}}
    rate::Float64
    average_number::Vector{Float64}
end

mutable struct NeutronCapture <: AbstractReaction # Temperature Dependent Reactions
    reactant::Vector{Vector{Int64}} # [[Z_0, N_0], [Z_1, N_1], ... [Z_n, N_n]]
    product::Vector{Vector{Int64}}
    rate::Interpolations.Extrapolation
    pfunc::Interpolations.Extrapolation
    current_rate::Float64 
    q::Union{Missing, Float64}
end

# TODO: Make every decay struct have the same structures (i.e. tuples, or vectors, or something els)
mutable struct AlphaDecay <: AbstractReaction
    reactant::Tuple{Int64, Int64} # (Z, N)
    product::Vector{Tuple{Int64, Int64}} # [(Z_0, N_0), (Z_1, N_1)]
    rate::Float64
end

struct Photodissociation <: AbstractReaction
    q::Float64
end

# TODO: Replace all dicts to simple vectors
mutable struct ReactionData
    probdecay::Vector{ProbDecay}
    neutroncapture::Dict{Vector{Vector{Int64}},NeutronCapture}
    alphadecay::Vector{AlphaDecay}
end

function initialize_reactions()
    probdecay = Vector{ProbDecay}()
    ncap_dict = Dict{Vector{Vector{Int64}},NeutronCapture}()
    alphadecay = Vector{AlphaDecay}()
    return ReactionData(probdecay, ncap_dict, alphadecay)
end

end