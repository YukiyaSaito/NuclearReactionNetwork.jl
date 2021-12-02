"""
    ReactionTypes

Provides a series of reaction types as well as some functionality for them.
"""
module ReactionTypes

using StaticArrays
using ..LinearInterpolations
using ..Astro
using ..Network

export ProbDecay
export ProbDecayIO
export NeutronCapture
export NeutronCaptureIO
export AlphaDecay
export AlphaDecayIO
export Photodissociation
export PhotodissociationIO
export ReactionData
export ReactionDataIO
export ProbRxn
export ProbRxnIO
export Rxn
export RxnIO
export Decay
export DecayIO
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
mutable struct ProbDecayIO <: AbstractReaction
    """A vector of reactants of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    reactant::SVector{1, Tuple{Int, Int}}
    """A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    product::Vector{Tuple{Int, Int}}
    """The rate of the reaction."""
    rate::Float64
    """A vector of average numbers of species produced."""
    average_number::Vector{Float64}
end

struct ProbDecay <: AbstractReaction
    reactant_idxs::SVector{1, Int}
    product_idxs::Vector{Int}
    average_number::Vector{Float64}
    rate::Float64
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
mutable struct NeutronCaptureIO <: AbstractReaction # Temperature Dependent Reactions
    """A vector of reactants of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    reactant::SVector{2, Tuple{Int, Int}}
    """A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1), ... (Z_n, N_n)]``"""
    product::SVector{1, Tuple{Int, Int}}
    """The rates and partition function data to be interpolated."""
    rates_pfuncs_lerp::NcapLerp
    """Unused field."""
    current_rate::Float64 
    """The Q value of the reverse reactionn. Used for the calculation of the reverse rate."""
    q::Union{Nothing, Float64}
end

struct NeutronCapture <: AbstractReaction
    reactant_idxs::SVector{2, Int}
    product_idxs::SVector{1, Int}
    rates_pfuncs_lerp::NcapLerp
    q::Union{Nothing, Float64}
    A_factor::Float64
end

"""
    AlphaDecay

``\\alpha`` decay

# Fields:
- `reactant::Vector{Tuple{Int, Int}}`: A tuple of reactants of the form ``(Z, N)``
- `product::Vector{Tuple{Int, Int}}`: A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1)]``
- `rate::Float64`: The rate of the reaction.
"""
mutable struct AlphaDecayIO <: AbstractReaction
    """A tuple of reactants of the form ``(Z, N)``"""
    reactant::Tuple{Int, Int}
    """A vector of products of the form ``[(Z_0, N_0), (Z_1, N_1)]``"""
    product::SVector{2, Tuple{Int, Int}}
    """The rate of the reaction."""
    rate::Float64
end

struct AlphaDecay <: AbstractReaction
    reactant_idx::Int
    product_idxs::SVector{2, Int}
    rate::Float64
end

mutable struct ProbRxnIO <: AbstractReaction
    reactants::Vector{Tuple{Int, Int}}
    products::Vector{Tuple{Int, Int}}
    average_numbers::Vector{Float64}
    rates_lerp::RxnLerp
end

struct ProbRxn <: AbstractReaction
    reactant_idxs::Vector{Int}
    product_idxs::Vector{Int}
    average_numbers::Vector{Float64}
    rates_lerp::RxnLerp
end

struct RxnIO <: AbstractReaction
    reactants::Vector{Tuple{Int, Int}}
    products::Vector{Tuple{Int, Int}}
    rates_lerp::RxnLerp
end

struct Rxn <: AbstractReaction
    reactant_idxs::Vector{Int}
    product_idxs::Vector{Int}
    rates_lerp::RxnLerp
end

struct DecayIO <: AbstractReaction
    reactants::Vector{Tuple{Int, Int}}
    products::Vector{Tuple{Int, Int}}
    rate::Float64
end

struct Decay <: AbstractReaction
    reactant_idxs::Vector{Int}
    product_idxs::Vector{Int}
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

The reverse reaction of a neutron capture

# Fields:
- `q::Float64`: The Q value of the reverse reaction.
"""
struct Photodissociation
    q::Float64
end

mutable struct ReactionDataIO
    probdecay::Vector{Tuple{String, ProbDecayIO}}
    ncap_dict::Dict{Int, NeutronCaptureIO}
    alphadecay::Vector{AlphaDecayIO}
    probrxn::Vector{Tuple{String, ProbRxnIO}}
    rxn::Vector{RxnIO}
    decay::Vector{DecayIO}
end

"""
    ReactionData

A collection of all the reaction data in the network

# Fields
- `probdecay::Vector{ProbDecay}`: The probabilistic decays in the network, indexed byb `zn_to_index()`.
- `neutroncapture::Vec{Union{NeutronCapture}}`: The neutron captures in the network, indexed by `zn_to_index()`.
- `alphadecay::Vector{AlphaDecay}`: The alpha decays in the network, indexed by `zn_to_index()`.

See also: [`zn_to_index`](@ref)
"""
mutable struct ReactionData
    """The probabilistic decays in the network, indexed byb `zn_to_index()`."""
    probdecay::Vector{ProbDecay}
    """The neutron captures in the network, indexed by `zn_to_index()`."""
    neutroncapture::Vector{Union{Nothing, NeutronCapture}}
    """The alpha decays in the network, indexed by `zn_to_index()`."""
    alphadecay::Vector{AlphaDecay}
    probrxn::Vector{ProbRxn}
    rxn::Vector{Rxn}
    decay::Vector{Decay}
end

"""
    initialize_reactions()

Create an empty Reactiondata data structure

See also: [`ReactionData`](@ref)
"""
function initialize_reactions()
    probdecay = Vector{ProbDecay}()
    ncap = Vector{Union{Nothing, NeutronCapture}}()
    alphadecay = Vector{AlphaDecay}()
    probrxn = Vector{ProbRxn}()
    rxn = Vector{Rxn}()
    decay = Vector{Decay}()
    return ReactionData(probdecay, ncap, alphadecay, probrxn, rxn, decay)
end

end