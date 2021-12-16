"""
    NetworkDatas

This module handles everything to do with preparing for the Newton-Raphson iteration. In
here, we fill up the Jacobian, ``\\mathcal{J}``, as well as ``\\vec{\\dot{Y}}``. The
main data structure of this module is [`NetworkData`](@ref) which holds all of the
information needed to solve the network calculations.
"""
module NetworkDatas

using ..LinearInterpolations
using ..Astro
using ..ReactionTypes
using ..Network
using ..LinearSolvers
using SparseArrays
using LinearAlgebra
using DelimitedFiles

export RWData
export ROData
export NetworkData
export OutputInfo
export Checkpoint
export Modulator
export IncludedReactions
export fill_jacobian!
export update_ydot!

struct Modulator
    id::String
    rxn_type::Type
    value::Float64
end

struct Checkpoint
    path::Union{Missing, String}
    times::Vector{Float64}
    completed::BitVector
end

"""
    OutputInfo

Holds the data about where we should dump the output of the simulation.

# Fields:
- `dump_final_y::Bool`: Whether the simulation should output the final abundances for each species
- `final_y_path::Union{Missing, String}`: The file to which the simulation should output the final abundances for each species
- `dump_final_ya::Bool`: Whether the simulation should output the final abundances as a function of mass number
- `final_ya_path::Union{Missing, String}`: The file to which the simulation should output the final abundances as a fucntion of mass number
- `dump_each_iteration::Bool`: Whether the simulationn should output the abundances of each species at each time step
- `iteration_output_path::Union{Missing, String}`: The file to which the simulation shoould the abundances of each species at each time step
"""
struct OutputInfo
    """Whether the simulation should output the final abundances for each species"""
    dump_final_y::Bool
    """The file to which the simulation should output the final abundances for each species"""
    final_y_path::Union{Missing, String}
    """Whether the simulation should output the final abundances as a function of mass number"""
    dump_final_ya::Bool
    """The file to which the simulation should output the final abundances as a fucntion of mass number"""
    final_ya_path::Union{Missing, String}
    """Whether the simulationn should output the abundances of each species at each time step"""
    dump_each_iteration::Bool
    """The file to which the simulation shoould the abundances of each species at each time step"""
    iteration_output_path::Union{Missing, String}
    checkpoint::Checkpoint
end

"""
    IncludedReactions

Holds the data about which reactions are used in the network.

This structure exists for performance reasons.

# Fields:
- `ncap::Bool`: Whether neutron captures are included in the network calculations
- `probdecay::Bool`: Whether probabilistic decays are included in the network calculations
- `alphadecay::Bool`: Whether ``\\alpha`` decays are included in the network calculations
- `photodissociation::Bool`: Whether the reverse reaction of neutron captures are included in the network calculations
"""
mutable struct IncludedReactions
    """Whether neutron captures are included in the network calculations"""
    ncap::Bool
    """Whether probabilistic decays are included in the network calculations"""
    probdecay::Bool
    """Whether ``\\alpha`` decays are included in the network calculations"""
    alphadecay::Bool
    """Whether the reverse reaction of neutron captures are included in the network calculations"""
    photodissociation::Bool
    probrxn::Bool
    rxn::Bool
    decay::Bool
end

struct ROData
    """Used to check is a species is in a network."""
    net_idx::NetworkIndex
    """Holds all of the data used for the reactions (like rates, reactants, products, etc.)."""
    reaction_data::ReactionData
    """Holds the data for the astrophysical trajectory to interpolate (like density and temperature)."""
    trajectory::TrajectoryLerp
    """Used to know if there are reactions that we can skip for performance."""
    included_reactions::IncludedReactions
    """The solver used to invert the Jacobian."""
    solver::LinearSolver
end

struct RWData
    """The abundances of each species, indexed using `zn_to_index()`."""
    abundance::Vector{Float64}
    """The proposed abundances of each species if a Newton-Raphson iteration fails"""
    yproposed::Vector{Float64}
    """The time derivative of the abundances. The goal of the network solver is to integrate this over time."""
    ydot::Vector{Float64}
    """The current time, time step, and time at which the simulation should end."""
    time::Time
    """The Jacobian that we have to invert in order to integrate ``\\vec{\\dot{Y}}``."""
    jacobian::SparseMatrixCSC{Float64, Int}
    """Tells us where we should dump the output of the program to."""
    output_info::OutputInfo
    modulator::Modulator
end

"""
    NetworkData

A conglomerate of all the information needed for the network to run.

This struct should only ever be instantiated as a singleton.

# Fields:
- `net_idx::NetworkIndex`: Used to check is a species is in a network.
- `reaction_data::ReactionData`: Holds all of the data used for the reactions (like rates, reactants, products, etc.).
- `trajectory::TrajectoryLerp`: Holds the data for the astrophysical trajectory to interpolate (like density and temperature).
- `abundance::Vector{Float64}`: The abundances of each species, indexed using `zn_to_index()`.
- `yproposed::Vector{Float64}`: The proposed abundances of each species if a Newton-Raphson iteration fails
- `ydot::Vector{Float64}`: The time derivative of the abundances. The goal of the network solver is to integrate this over time.
- `time::Time`: The current time, time step, and time at which the simulation should end.
- `jacobian::SparseMatrixCSC{Float64, Int}`: The Jacobian that we have to invert in order to integrate ``\\vec{\\dot{Y}}``.
- `output_info::OutputInfo`: Tells us where we should dump the output of the program to.
- `included_reactions::IncludedReactions`: Used to know if there are reactions that we can skip for performance.
- `solver::LinearSolver`: The solver used to invert the Jacobian.

See also: [`zn_to_index`](@ref), [`NetworkIndex`](@ref), [`ReactionData`](@ref), [`TrajectoryLerp`](@ref),
[`Time`](@ref), [`OutputInfo`](@ref), [`IncludedReactions`](@ref),
[`LinearSolver`](@ref)
"""
struct NetworkData
    rd::ROData
    wd::RWData
end

@inbounds function initialize_ydot!(nd::NetworkData)::Vector{Float64}
    fill!(nd.wd.ydot, 0.0)
    # Threads.@threads for i in eachindex(nd.wd.ydot)
    #     nd.wd.ydot[i] = 0.0
    # end
    return nd.wd.ydot
end

@inbounds function fill_probdecay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == ProbDecay

    for decay::ProbDecay in nd.rd.reaction_data.probdecay
        rate::Float64 = decay.rate
        if iszero(rate)
            continue
        end

        if use_modulation && decay.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        reactant_idx::Int = decay.reactant_idxs[1]

        if iszero(abundance[reactant_idx])
            continue
        end

        nd.wd.ydot[reactant_idx] += -1.0 * rate * abundance[reactant_idx]

        for (product_idx::Int, average_number::Float64) in zip(decay.product_idxs, decay.average_number)
            if iszero(average_number)
                continue
            end
            nd.wd.ydot[product_idx] += average_number * rate * abundance[reactant_idx]
        end
    end
end

@inbounds function fill_neutroncapture_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == NeutronCapture

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)
    for capture::Union{Nothing, NeutronCapture} in nd.rd.reaction_data.neutroncapture
        if isnothing(capture)
            continue
        end

        rate::Float64 = get_rate(capture.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        if use_modulation && capture.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        # Grab the product of all the abundances
        abundance_factor::Float64 = 1.0
        for reactant_idx::Int in capture.reactant_idxs
            abundance_factor *= abundance[reactant_idx]
        end
        if iszero(abundance_factor)
            continue
        end

        # Double counting factor
        dc_factor::Float64 = capture.reactant_idxs[1] == capture.reactant_idxs[2] ? 0.5 : 1.0

        # Update ydot for the reactants
        for reactant_idx::Int in capture.reactant_idxs
            nd.wd.ydot[reactant_idx] += -1.0 * dc_factor * curr_traj.density * rate * abundance_factor
        end

        # Update ydot for the products
        for product_idx::Int in capture.product_idxs
            nd.wd.ydot[product_idx] += rate * dc_factor * curr_traj.density * abundance_factor
        end
    end
end

@inbounds function fill_alphadecay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == AlphaDecay

    for decay::AlphaDecay in nd.rd.reaction_data.alphadecay
        rate::Float64 = decay.rate

        if iszero(rate)
            continue
        end

        if use_modulation && decay.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        # Grab the abundace of the reactant
        reactant_idx::Int = decay.reactant_idx
        abundance_factor::Float64 = abundance[reactant_idx]
        if iszero(abundance_factor)
            continue
        end

        # Update ydot for the reactant
        nd.wd.ydot[reactant_idx] += -1.0 * rate * abundance_factor

        # Update ydot for the products
        for product_idx::Int in decay.product_idxs
            nd.wd.ydot[product_idx] += 1.0 * rate * abundance_factor
        end
    end
end

@inbounds function fill_photodissociation_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)
    for reaction::Union{Nothing, NeutronCapture} in nd.rd.reaction_data.neutroncapture
        if isnothing(reaction)
            continue
        end

        q::Union{Nothing, Float64} = reaction.q
        if isnothing(reaction.q)
            continue
        end

        forward_rate::Float64, pfunc_p::Float64 = get_rate_pfunc(reaction.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(forward_rate)
            continue
        end

        # Grab the abundace of the reactant
        reactant_idx::Int = reaction.product_idxs[1]
        abundance_factor::Float64 = abundance[reactant_idx]
        if iszero(abundance_factor)
            continue
        end

        # Lookup partition function for the reactant
        if reactant_idx > length(nd.rd.reaction_data.neutroncapture) || isnothing(nd.rd.reaction_data.neutroncapture[reactant_idx])
            pfunc_r::Float64 = 1.0
        else
            reactant = nd.rd.reaction_data.neutroncapture[reactant_idx]
            pfunc_r = get_pfunc(reactant.rates_pfuncs_lerp, curr_traj.temperature)
        end
        pfunc_n::Float64 = 2.0

        pfunc::Float64 = pfunc_n * pfunc_p / pfunc_r # FIXME: What about division by zero?
        if iszero(pfunc)
            continue
        end

        # Compute the rate of the reaction
        if q < 0
            reverse_rate::Float64 = 1e16
        else
            reverse_rate = forward_rate * 9.8678e9 * pfunc * reaction.A_factor * curr_traj.temperature^(3/2) * exp(-11.605 * q / curr_traj.temperature)
        end

        # Update ydot for the reactant
        nd.wd.ydot[reactant_idx] += -1.0 * reverse_rate * abundance_factor

        # Update ydot for the products
        for product_idx::Int in reaction.reactant_idxs
            nd.wd.ydot[product_idx] += 1.0 * reverse_rate * abundance_factor
        end
    end
end

@inbounds function fill_probrxn_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == ProbRxn

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)

    for rxn::ProbRxn in nd.rd.reaction_data.probrxn
        # Grab the rate of the reaction
        rate::Float64 = get_rate(rxn.rates_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        if use_modulation && rxn.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        density_factor::Float64 = curr_traj.density^(length(rxn.reactant_idxs) - 1)

        # Grab the abundance of the reactants
        abundance_factor::Float64 = prod(abundance[rxn.reactant_idxs])
        if iszero(abundance_factor)
            continue
        end

        # Double counting factor
        # num_same = length(rxn.reactant_idxs) - length(unique(rxn.reactant_idxs))
        # dc_factor::Float64 = 1/factorial(num_same)
        dc_factor::Float64 = 1.0 #rxn.reactant_idxs[1] == rxn.reactant_idxs[2] ? 0.5 : 1.0 # FIXME: This does not generalize

        for reactant_idx::Int in rxn.reactant_idxs
            nd.wd.ydot[reactant_idx] += -1.0 * dc_factor * rate * density_factor * abundance_factor
        end
        for (product_idx::Int, average_number::Float64) in zip(rxn.product_idxs, rxn.average_numbers)
            if iszero(average_number)
                continue
            end
            nd.wd.ydot[product_idx] += dc_factor * average_number * rate * density_factor * abundance_factor
        end
    end
end

@inbounds function fill_rxn_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == Rxn

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)

    for rxn::Rxn in nd.rd.reaction_data.rxn
        # Grab the rate of the reaction
        rate::Float64 = get_rate(rxn.rates_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        if use_modulation && rxn.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        density_factor::Float64 = curr_traj.density^(length(rxn.reactant_idxs) - 1)

        # Grab the abundance of the reactants
        abundance_factor::Float64 = prod(abundance[rxn.reactant_idxs])
        if iszero(abundance_factor)
            continue
        end

        # Double counting factor
        # num_same = length(rxn.reactant_idxs) - length(unique(rxn.reactant_idxs))
        # dc_factor::Float64 = 1/factorial(num_same)
        dc_factor::Float64  = 1.0

        for reactant_idx::Int in rxn.reactant_idxs
            nd.wd.ydot[reactant_idx] += -1.0 * dc_factor * rate * density_factor * abundance_factor
        end
        for product_idx::Int in rxn.product_idxs
            nd.wd.ydot[product_idx] += dc_factor * rate * density_factor * abundance_factor
        end
    end
end

@inbounds function fill_decay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == Decay

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)

    for decay::Decay in nd.rd.reaction_data.decay
        # Grab the rate of the reaction
        rate::Float64 = decay.rate
        if iszero(rate)
            continue
        end

        if use_modulation && decay.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        density_factor::Float64 = curr_traj.density^(length(decay.reactant_idxs) - 1)

        # Grab the abundance of the reactants
        abundance_factor::Float64 = prod(abundance[decay.reactant_idxs])
        if iszero(abundance_factor)
            continue
        end

        # Double counting factor
        # num_same = length(decay.reactant_idxs) - length(unique(decay.reactant_idxs))
        # dc_factor::Float64 = 1/factorial(num_same)
        dc_factor::Float64  = 1.0

        for reactant_idx::Int in decay.reactant_idxs
            nd.wd.ydot[reactant_idx] += -1.0 * dc_factor * rate * density_factor * abundance_factor
        end
        for product_idx::Int in decay.product_idxs
            nd.wd.ydot[product_idx] += dc_factor * rate * density_factor * abundance_factor
        end
    end
end

"""
    update_ydot!(nd::NetworkData; use_yproposed::Bool=false)::Nothing

Fill the vector ``\\vec{\\dot{Y}}`` with the data from the current time-step.

The flag `use_yproposed` is used to determine whether the network should use the abundances
held in `nd.wd.abundance` or `nd.wd.yproposed`. The latter should be used for an iteration that
failed to converge after a Newton-Raphson iteration.
"""
function update_ydot!(nd::NetworkData; use_yproposed::Bool=false)::Nothing
    initialize_ydot!(nd)
    if nd.rd.included_reactions.photodissociation
        fill_photodissociation_ydot!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.ncap
        fill_neutroncapture_ydot!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.alphadecay
        fill_alphadecay_ydot!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.probdecay
        fill_probdecay_ydot!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.probrxn
        fill_probrxn_ydot!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.rxn
        fill_rxn_ydot!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.decay
        fill_decay_ydot!(nd, use_yproposed)
    end
end

@inbounds function fill_jacobian_probdecay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == ProbDecay

    for decay::ProbDecay in nd.rd.reaction_data.probdecay
        rate::Float64 = decay.rate
        if iszero(rate)
            continue
        end

        if use_modulation && decay.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        reactant_idx::Int = decay.reactant_idxs[1]
        nd.wd.jacobian[reactant_idx, reactant_idx] += -1.0 * rate

        for (product_idx::Int, average_number::Float64) in zip(decay.product_idxs, decay.average_number)
            if iszero(average_number)
                continue
            end
            nd.wd.jacobian[product_idx, reactant_idx] += average_number * rate
        end
    end
end

@inbounds function fill_jacobian_neutroncapture!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == NeutronCapture

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)
    if iszero(curr_traj.density)
        return
    end

    # Neutron index and abundance
    neutron_idx::Int = zn_to_index(Int(0), Int(1), nd.rd.net_idx)
    neutron_abundance::Float64 = abundance[neutron_idx]

    # TODO: Convert this to a loop instead of 6 hard coded changes to the jacobian?
    for capture::Union{Nothing, NeutronCapture} in nd.rd.reaction_data.neutroncapture
        if isnothing(capture)
            continue
        end

        rate::Float64 = get_rate(capture.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        if use_modulation && capture.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        # Reactant index and abundance
        reactant_idx::Int = capture.reactant_idxs[2]
        reactant_abundance::Float64 = abundance[reactant_idx]

        # Product index
        product_idx::Int = capture.product_idxs[1]

        # Double counting factor TODO: Is this always 1.0?
        dc_factor::Float64 = reactant_idx == neutron_idx ? 0.5 : 1.0

        # Reactants
        nd.wd.jacobian[reactant_idx, reactant_idx] -= dc_factor * curr_traj.density * rate * neutron_abundance  # J_RR
        nd.wd.jacobian[reactant_idx, neutron_idx]  -= dc_factor * curr_traj.density * rate * reactant_abundance # J_Rn
        nd.wd.jacobian[neutron_idx, reactant_idx]  -= dc_factor * curr_traj.density * rate * neutron_abundance  # J_nR
        nd.wd.jacobian[neutron_idx, neutron_idx]   -= dc_factor * curr_traj.density * rate * reactant_abundance # J_nn

        # Products
        nd.wd.jacobian[product_idx, reactant_idx]  += dc_factor * curr_traj.density * rate * neutron_abundance  # J_PR
        nd.wd.jacobian[product_idx, neutron_idx]   += dc_factor * curr_traj.density * rate * reactant_abundance # J_Pn
    end
end

@inbounds function fill_jacobian_alphadecay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == AlphaDecay

    for decay::AlphaDecay in nd.rd.reaction_data.alphadecay
        rate::Float64 = decay.rate

        if use_modulation && decay.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        reactant_idx::Int = decay.reactant_idx
        nd.wd.jacobian[reactant_idx, reactant_idx] += -1.0 * rate

        for product_idx::Int in decay.product_idxs
            nd.wd.jacobian[product_idx, reactant_idx] += rate
        end
    end
end

@inbounds function fill_jacobian_photodissociation!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)
    if iszero(curr_traj.density)
        return
    end

    # TODO: Convert this to a loop instead of 6 hard coded changes to the jacobian?
    for reaction::Union{Nothing, NeutronCapture} in nd.rd.reaction_data.neutroncapture
        if isnothing(reaction)
            continue
        end

        forward_rate::Float64 = get_rate(reaction.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(forward_rate)
            continue
        end

        q::Union{Nothing, Float64} = reaction.q
        if isnothing(q)
            continue
        end

        # Reactant index
        reactant_idx::Int = reaction.product_idxs[1]

        # Neutron (product) index
        neutron_idx::Int = zn_to_index(Int(0), Int(1), nd.rd.net_idx)

        # Product index
        product_idx::Int = reaction.reactant_idxs[2]

        # Lookup partition function for the reactant (product of the forward reaction)
        if reactant_idx > length(nd.rd.reaction_data.neutroncapture) || isnothing(nd.rd.reaction_data.neutroncapture[reactant_idx])
            pfunc_r::Float64 = 1.0
        else
            reactant = nd.rd.reaction_data.neutroncapture[reactant_idx]
            pfunc_r = get_pfunc(reactant.rates_pfuncs_lerp, curr_traj.temperature)
        end
        pfunc_n::Float64 = 2.0
        pfunc_p::Float64 = get_pfunc(reaction.rates_pfuncs_lerp, curr_traj.temperature)

        pfunc = pfunc_n * pfunc_p / pfunc_r # FIXME: What about division by zero?
        if iszero(pfunc)
            continue
        end

        # Reverse rate
        if q < 0
            reverse_rate::Float64 = 1e16
        else
            reverse_rate = forward_rate * 9.8678e9 * pfunc * reaction.A_factor * curr_traj.temperature^(3/2) * exp(-11.605 * q / curr_traj.temperature)
        end

        # Reactant
        nd.wd.jacobian[reactant_idx, reactant_idx] -= reverse_rate # J_RR

        # Products
        nd.wd.jacobian[product_idx, reactant_idx]  += reverse_rate # J_PR
        nd.wd.jacobian[neutron_idx, reactant_idx]  += reverse_rate # J_nR
    end
end

@inbounds function fill_jacobian_probrxn!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == ProbRxn

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)

    for rxn::ProbRxn in nd.rd.reaction_data.probrxn
        # Grab the rate of the reaction
        rate::Float64 = get_rate(rxn.rates_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        if use_modulation && rxn.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        density_factor::Float64 = curr_traj.density^(length(rxn.reactant_idxs) - 1)

        # Double counting factor
        # num_same = length(rxn.reactant_idxs) - length(unique(rxn.reactant_idxs))
        # dc_factor::Float64 = 1/factorial(num_same)
        dc_factor::Float64 = 1.0 #rxn.reactant_idxs[1] == rxn.reactant_idxs[2] ? 0.5 : 1.0 # FIXME: This does not generalize

        for (i::Int, r1_idx::Int) in enumerate(rxn.reactant_idxs)
            abundance_factor::Float64 = 1.0
            for (j::Int, r2_idx::Int) in enumerate(rxn.reactant_idxs)
                if i != j
                    abundance_factor *= abundance[r2_idx]
                end
            end
            if iszero(abundance_factor)
                continue
            end

            for r2_idx::Int in rxn.reactant_idxs
                nd.wd.jacobian[r2_idx, r1_idx] += -1.0 * dc_factor * rate * density_factor * abundance_factor
            end

            for (product_idx::Int, average_number::Float64) in zip(rxn.product_idxs, rxn.average_numbers)
                if iszero(average_number)
                    continue
                end
                nd.wd.jacobian[product_idx, r1_idx] += dc_factor * average_number * rate * density_factor * abundance_factor
            end
        end
    end
end

@inbounds function fill_jacobian_rxn!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == Rxn

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)

    for rxn::Rxn in nd.rd.reaction_data.rxn
        # Grab the rate of the reaction
        rate::Float64 = get_rate(rxn.rates_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        if use_modulation && rxn.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        density_factor::Float64 = curr_traj.density^(length(rxn.reactant_idxs) - 1)

        # Double counting factor
        # num_same = length(rxn.reactant_idxs) - length(unique(rxn.reactant_idxs))
        # dc_factor::Float64 = 1/factorial(num_same)
        dc_factor::Float64 = 1.0

        for (i::Int, r1_idx::Int) in enumerate(rxn.reactant_idxs)
            abundance_factor::Float64 = 1.0
            for (j::Int, r2_idx::Int) in enumerate(rxn.reactant_idxs)
                if i != j
                    abundance_factor *= abundance[r2_idx]
                end
            end
            if iszero(abundance_factor)
                continue
            end

            for r2_idx::Int in rxn.reactant_idxs
                nd.wd.jacobian[r2_idx, r1_idx] += -1.0 * dc_factor * rate * density_factor * abundance_factor
            end

            for product_idx::Int in rxn.product_idxs
                nd.wd.jacobian[product_idx, r1_idx] += dc_factor * rate * density_factor * abundance_factor
            end
        end
    end
end

@inbounds function fill_jacobian_decay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.wd.abundance
    if use_yproposed
        abundance = nd.wd.yproposed
    end

    use_modulation = nd.wd.modulator.rxn_type == Decay

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.rd.trajectory, nd.wd.time.current)

    for decay::Decay in nd.rd.reaction_data.decay
        # Grab the rate of the reaction
        rate::Float64 = decay.rate
        if iszero(rate)
            continue
        end

        if use_modulation && decay.id == nd.wd.modulator.id
            rate *= nd.wd.modulator.value
        end

        density_factor::Float64 = curr_traj.density^(length(decay.reactant_idxs) - 1)

        # Double counting factor
        # num_same = length(decay.reactant_idxs) - length(unique(decay.reactant_idxs))
        # dc_factor::Float64 = 1/factorial(num_same)
        dc_factor::Float64 = 1.0

        for (i::Int, r1_idx::Int) in enumerate(decay.reactant_idxs)
            abundance_factor::Float64 = 1.0
            for (j::Int, r2_idx::Int) in enumerate(decay.reactant_idxs)
                if i != j
                    abundance_factor *= abundance[r2_idx]
                end
            end
            if iszero(abundance_factor)
                continue
            end

            for r2_idx::Int in decay.reactant_idxs
                nd.wd.jacobian[r2_idx, r1_idx] += -1.0 * dc_factor * rate * density_factor * abundance_factor
            end

            for product_idx::Int in decay.product_idxs
                nd.wd.jacobian[product_idx, r1_idx] += dc_factor * rate * density_factor * abundance_factor
            end
        end
    end
end

"""
    fill_jacobian!(nd::NetworkData; use_yproposed::Bool=false)::Vector{Float64}

Fill Jacobian with the data from the current time-step.

The flag `use_yproposed` is used to determine whether the network should use the abundances
held in `nd.wd.abundance` or `nd.wd.yproposed`. The latter should be used for an iteration that
failed to converge after a Newton-Raphson iteration.
"""
function fill_jacobian!(nd::NetworkData; use_yproposed::Bool=false)::Vector{Float64}
    # Jacobian coordinate: (reactant, product)

    mul!(nd.wd.jacobian, nd.wd.jacobian, 0.0)

    if nd.rd.included_reactions.photodissociation
        fill_jacobian_photodissociation!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.ncap
        fill_jacobian_neutroncapture!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.alphadecay
        fill_jacobian_alphadecay!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.probdecay
        fill_jacobian_probdecay!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.probrxn
        fill_jacobian_probrxn!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.rxn
        fill_jacobian_rxn!(nd, use_yproposed)
    end
    if nd.rd.included_reactions.decay
        fill_jacobian_decay!(nd, use_yproposed)
    end

    mul!(nd.wd.jacobian, nd.wd.jacobian, -1.0)
    nd.wd.jacobian[diagind(nd.wd.jacobian)] .+= 1.0/nd.wd.time.step
end

end