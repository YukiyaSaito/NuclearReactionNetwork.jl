module ReactionTypes

using ..LinearInterpolations
using ..Astro
using ..Network

export ProbDecay
export NeutronCapture
export AlphaDecay
export Photodissociation
export ReactionData
export initialize_reactions
export check_eq_reaction

abstract type AbstractReaction end

# mutable struct NBodyReaction <: AbstractReaction
#     reactant::Vector{Vector{Int}}
#     product::Vector{Vector{Int}}
# end

mutable struct ProbDecay <: AbstractReaction
    reactant::Vector{Tuple{Int, Int}} # [(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]
    product::Vector{Tuple{Int, Int}}
    rate::Float64
    average_number::Vector{Float64}
end

mutable struct NeutronCapture <: AbstractReaction # Temperature Dependent Reactions
    reactant::Vector{Tuple{Int, Int}} # [(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]
    product::Vector{Tuple{Int, Int}}
    rates_pfuncs_lerp::NcapLerp
    current_rate::Float64 
    q::Union{Missing, Float64}
end

mutable struct AlphaDecay <: AbstractReaction
    reactant::Tuple{Int, Int} # (Z, N)
    product::Vector{Tuple{Int, Int}} # [(Z_0, N_0), (Z_1, N_1)]
    rate::Float64
end

# TODO: Should this be == instead?
function check_eq_reaction(lhs::AbstractReaction, rhs::AbstractReaction)::Bool
    return (typeof(lhs) == typeof(rhs)) && (lhs.reactant == rhs.reactant) && (lhs.product == rhs.product)
end

struct Photodissociation
    q::Float64
end

mutable struct ReactionData
    probdecay::Vector{ProbDecay}
    neutroncapture::Dict{Int, NeutronCapture}
    alphadecay::Vector{AlphaDecay}
end

function initialize_reactions()
    probdecay = Vector{ProbDecay}()
    ncap_dict = Dict{Int, NeutronCapture}()
    alphadecay = Vector{AlphaDecay}()
    return ReactionData(probdecay, ncap_dict, alphadecay)
end

end