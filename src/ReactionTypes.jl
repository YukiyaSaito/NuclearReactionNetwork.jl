"""
    ReactionTypes

Provides a series of reaction types as well as some functionality for them.
"""
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

"""A generic reaction type to refer to any reaction."""
abstract type AbstractReaction end

# mutable struct NBodyReaction <: AbstractReaction
#     reactant::Vector{Vector{Int}}
#     product::Vector{Vector{Int}}
# end

"""
    ProbDecay <: AbstractReaction

Probabilistic decay (like ``\\beta^-``decay)

# Fields:
- `reactant::Vector{Tuple{Int, Int}}`: A vector of reactants of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``
- `product::Vector{Tuple{Int, Int}}`: A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``
- `rate::Float64`: The rate of the reaction.
- `average_number::Vector{Float64}`: A vector of average numbers of species produced.
"""
mutable struct ProbDecay <: AbstractReaction
    """A vector of reactants of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    reactant::Vector{Tuple{Int, Int}}
    """A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    product::Vector{Tuple{Int, Int}}
    """The rate of the reaction."""
    rate::Float64
    """A vector of average numbers of species produced."""
    average_number::Vector{Float64}
end

"""
    NeutronCapture

Neutron capture reaction

# Fields:
- `reactant::Vector{Tuple{Int, Int}}`: A vector of reactants of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``
- `product::Vector{Tuple{Int, Int}}`: A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``
- `rates_pfuncs_lerp::NcapLerp`: The rates and partition function data to be interpolated.
- `current_rate::Float64`: Unused field.
- `q::Union{Missinng, Float64}`: The Q value of the reverse reactionn. Used for the calculation of the reverse rate.
"""
mutable struct NeutronCapture <: AbstractReaction # Temperature Dependent Reactions
    """A vector of reactants of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    reactant::Vector{Tuple{Int, Int}}
    """A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    product::Vector{Tuple{Int, Int}}
    """The rates and partition function data to be interpolated."""
    rates_pfuncs_lerp::NcapLerp
    """Unused field."""
    current_rate::Float64 
    """The Q value of the reverse reactionn. Used for the calculation of the reverse rate."""
    q::Union{Missing, Float64}
end

"""
    AlphaDecay

``\\alpha`` decay

# Fields:
- `reactant::Vector{Tuple{Int, Int}}`: A tuple of reactants of the form ``(Z, N)``
- `product::Vector{Tuple{Int, Int}}`: A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1)]``
- `rate::Float64`: The rate of the reaction.
"""
mutable struct AlphaDecay <: AbstractReaction
    """A tuple of reactants of the form ``(Z, N)``"""
    reactant::Tuple{Int, Int}
    """A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1)]``"""
    product::Vector{Tuple{Int, Int}}
    """The rate of the reaction."""
    rate::Float64
end

# TODO: Should this be == instead?
"""
    check_eq_reaction(lhs::AbstractReaction, rhs::AbstractReaction)::Bool

Check if `lhs == rhs`, where `lhs` and `rhs` are reactions.

This only checks if the reactants and the products are the same.
"""
function check_eq_reaction(lhs::AbstractReaction, rhs::AbstractReaction)::Bool
    return (typeof(lhs) == typeof(rhs)) && (lhs.reactant == rhs.reactant) && (lhs.product == rhs.product)
end

"""
    Photodissociation

The reverse reaction of a neutronn capture

# Fields:
- `q::Float64`: The Q value of the reverse reaction.
"""
struct Photodissociation
    q::Float64
end

"""
    ReactionData

A collection of all the reaction data in the network

# Fields
- `probdecay::Vector{ProbDecay}`: The probabilistic decays in the network
- `neutroncapture::Dict{Int, NeutronCapture}`: The neutron captures in the network, indexed by `zn_to_index()`
- `alphadecay::Vector{AlphaDecay}`: The alpha decays in the network.

See also: [`zn_to_index`](@ref)
"""
mutable struct ReactionData
    probdecay::Vector{ProbDecay}
    neutroncapture::Dict{Int, NeutronCapture}
    alphadecay::Vector{AlphaDecay}
end

"""
    initialize_reactions()

Create an empty Reactiondata data structure

See also: [`ReactionData`](@ref)
"""
function initialize_reactions()
    probdecay = Vector{ProbDecay}()
    ncap_dict = Dict{Int, NeutronCapture}()
    alphadecay = Vector{AlphaDecay}()
    return ReactionData(probdecay, ncap_dict, alphadecay)
end

end