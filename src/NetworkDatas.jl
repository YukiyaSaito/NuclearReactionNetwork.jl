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

export NetworkData
export OutputInfo
export IncludedReactions
export fill_jacobian!
export update_ydot!

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
    """Used to check is a species is in a network."""
    net_idx::NetworkIndex
    """Holds all of the data used for the reactions (like rates, reactants, products, etc.)."""
    reaction_data::ReactionData
    """Holds the data for the astrophysical trajectory to interpolate (like density and temperature)."""
    trajectory::TrajectoryLerp
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
    """Used to know if there are reactions that we can skip for performance."""
    included_reactions::IncludedReactions
    """The solver used to invert the Jacobian."""
    solver::LinearSolver
end

@inbounds function initialize_ydot!(nd::NetworkData)::Vector{Float64}
    # fill!(nd.ydot, 0.0)
    Threads.@threads for i in eachindex(nd.ydot)
        nd.ydot[i] = 0.0
    end
    return nd.ydot
end

@inbounds function fill_probdecay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    for decay::ProbDecay in nd.reaction_data.probdecay
        z_r::Int, n_r::Int = decay.reactant[1]
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)

        if iszero(decay.rate) || iszero(abundance[reactant_idx])
            continue
        end

        nd.ydot[reactant_idx] += -1.0 * decay.rate * abundance[reactant_idx]

        for (product::Tuple{Int, Int}, average_number::Float64) in zip(decay.product, decay.average_number)
            z_p::Int, n_p::Int = product
            product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)
            nd.ydot[product_idx] += average_number * decay.rate * abundance[reactant_idx]
        end
    end
end

@inbounds function fill_neutroncapture_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.trajectory, nd.time.current)
    for capture::NeutronCapture in values(nd.reaction_data.neutroncapture)
        rate::Float64 = get_rate(capture.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        # Grab the product of all the abundances
        abundance_factor::Float64 = 1.0
        for reactant::Tuple{Int, Int} in capture.reactant
            z_r::Int, n_r::Int = reactant
            reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
            abundance_factor *= abundance[reactant_idx]
        end
        if iszero(abundance_factor)
            continue
        end

        # Update ydot for the reactants
        for reactant::Tuple{Int, Int} in capture.reactant
            z_r::Int, n_r::Int = reactant
            reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
            nd.ydot[reactant_idx] += -1.0 * curr_traj.density * rate * abundance_factor
        end

        # Update ydot for the products
        for product::Tuple{Int, Int} in capture.product
            z_p::Int, n_p::Int = product
            product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)
            nd.ydot[product_idx] += rate * curr_traj.density * abundance_factor
        end
    end
end

@inbounds function fill_alphadecay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    for decay::AlphaDecay in nd.reaction_data.alphadecay
        rate::Float64 = decay.rate
        if iszero(rate)
            continue
        end

        # Grab the abundace of the reactant
        z_r::Int, n_r::Int = decay.reactant
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
        abundance_factor::Float64 = abundance[reactant_idx]
        if iszero(abundance_factor)
            continue
        end

        # Update ydot for the reactant
        nd.ydot[reactant_idx] += -1.0 * rate * abundance_factor

        # Update ydot for the products
        for product::Tuple{Int, Int} in decay.product
            z_p::Int, n_p::Int = product
            product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)
            nd.ydot[product_idx] += 1.0 * rate * abundance_factor
        end
    end
end

@inbounds function fill_photodissociation_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.trajectory, nd.time.current)
    for reaction::NeutronCapture in values(nd.reaction_data.neutroncapture)
        q::Union{Missing, Float64} = reaction.q
        if ismissing(reaction.q)
            continue
        end

        forward_rate::Float64 = get_rate(reaction.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(forward_rate)
            continue
        end

        # Grab the abundace of the reactant
        z_r::Int, n_r::Int = reaction.product[1]
        A_r::Int = z_r + n_r
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
        abundance_factor::Float64 = abundance[reactant_idx]
        if iszero(abundance_factor)
            continue
        end

        # Lookup partition function for the reactant
        if !haskey(nd.reaction_data.neutroncapture, reactant_idx)
            pfunc_r::Float64 = 1.0
        else
            reactant = nd.reaction_data.neutroncapture[reactant_idx]
            pfunc_r = get_pfunc(reactant.rates_pfuncs_lerp, curr_traj.temperature)
        end
        pfunc_n::Float64 = 2.0
        pfunc_p::Float64 = get_pfunc(reaction.rates_pfuncs_lerp, curr_traj.temperature)

        pfunc::Float64 = pfunc_n * pfunc_p / pfunc_r # FIXME: What about division by zero?
        if iszero(pfunc)
            continue
        end

        # Compute the rate of the reaction
        A_p::Int = 0
        for product::Tuple{Int, Int} in reaction.reactant
            z_p::Int, n_p::Int = product
            A_p += z_p + n_p
        end
        A_n::Int = 1
        if q < 0
            reverse_rate::Float64 = 1e16
        else
            reverse_rate = forward_rate * 9.8678e9 * pfunc * (A_n*A_p/A_r)^(3/2) * curr_traj.temperature^(3/2) * exp(-11.605 * q / curr_traj.temperature)
        end

        # Update ydot for the reactant
        nd.ydot[reactant_idx] += -1.0 * reverse_rate * abundance_factor

        # Update ydot for the products
        for product::Tuple{Int, Int} in reaction.reactant
            z_p::Int, n_p::Int = product
            product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)
            nd.ydot[product_idx] += 1.0 * reverse_rate * abundance_factor
        end
    end
end

"""
    update_ydot!(nd::NetworkData; use_yproposed::Bool=false)::Nothing

Fill the vector ``\\vec{\\dot{Y}}`` with the data from the current time-step.

The flag `use_yproposed` is used to determine whether the network should use the abundances
held in `nd.abundance` or `nd.yproposed`. The latter should be used for an iteration that
failed to converge after a Newton-Raphson iteration.
"""
function update_ydot!(nd::NetworkData; use_yproposed::Bool=false)::Nothing
    initialize_ydot!(nd)
    if nd.included_reactions.photodissociation
        fill_photodissociation_ydot!(nd, use_yproposed)
    end
    if nd.included_reactions.ncap
        fill_neutroncapture_ydot!(nd, use_yproposed)
    end
    if nd.included_reactions.alphadecay
        fill_alphadecay_ydot!(nd, use_yproposed)
    end
    if nd.included_reactions.probdecay
        fill_probdecay_ydot!(nd, use_yproposed)
    end
end

@inbounds function fill_jacobian_probdecay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    for decay::ProbDecay in nd.reaction_data.probdecay
        if iszero(decay.rate)
            continue
        end

        z_r::Int, n_r::Int = decay.reactant[1]
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
        nd.jacobian[reactant_idx, reactant_idx] += -1.0 * decay.rate

        for (product::Tuple{Int, Int}, average_number::Float64) in zip(decay.product, decay.average_number)
            z_p::Int, n_p::Int = product
            product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)
            nd.jacobian[product_idx, reactant_idx] += average_number * decay.rate
        end
    end
end

@inbounds function fill_jacobian_neutroncapture!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.trajectory, nd.time.current)
    if iszero(curr_traj.density)
        return
    end

    # TODO: Convert this to a loop instead of 6 hard coded changes to the jacobian?
    for capture::NeutronCapture in values(nd.reaction_data.neutroncapture)
        rate::Float64 = get_rate(capture.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(rate)
            continue
        end

        # Neutron index and abundance
        neutron_idx::Int = zn_to_index(Int(0), Int(1), nd.net_idx)
        neutron_abundance::Float64 = abundance[neutron_idx]

        # Reactant index and abundance
        reactant::Tuple{Int, Int} = capture.reactant[2]
        z_r::Int, n_r::Int = reactant
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
        reactant_abundance::Float64 = abundance[reactant_idx]

        # Product index
        product::Tuple{Int, Int} = capture.product[1]
        z_p::Int, n_p::Int = product
        product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)

        # Double counting factor TODO: Is this always 1.0?
        dc_factor::Float64 = reactant_idx == product_idx ? 0.5 : 1.0

        # Reactants
        nd.jacobian[reactant_idx, reactant_idx] -= dc_factor * curr_traj.density * rate * neutron_abundance  # J_RR
        nd.jacobian[reactant_idx, neutron_idx]  -= dc_factor * curr_traj.density * rate * reactant_abundance # J_Rn
        nd.jacobian[neutron_idx, reactant_idx]  -= dc_factor * curr_traj.density * rate * neutron_abundance  # J_nR
        nd.jacobian[neutron_idx, neutron_idx]   -= dc_factor * curr_traj.density * rate * reactant_abundance # J_nn

        # Products
        nd.jacobian[product_idx, reactant_idx]  += dc_factor * curr_traj.density * rate * neutron_abundance  # J_PR
        nd.jacobian[product_idx, neutron_idx]   += dc_factor * curr_traj.density * rate * reactant_abundance # J_Pn
    end
end

@inbounds function fill_jacobian_alphadecay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    for decay::AlphaDecay in nd.reaction_data.alphadecay
        z_r::Int, n_r::Int = decay.reactant
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)
        nd.jacobian[reactant_idx, reactant_idx] += -1.0 * decay.rate

        for product::Tuple{Int, Int} in decay.product
            z_p::Int, n_p::Int = product
            product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)
            nd.jacobian[product_idx, reactant_idx] += decay.rate
        end
    end
end

@inbounds function fill_jacobian_photodissociation!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.trajectory, nd.time.current)
    if iszero(curr_traj.density)
        return
    end

    # TODO: Convert this to a loop instead of 6 hard coded changes to the jacobian?
    for reaction::NeutronCapture in values(nd.reaction_data.neutroncapture)
        forward_rate::Float64 = get_rate(reaction.rates_pfuncs_lerp, curr_traj.temperature)
        if iszero(forward_rate)
            continue
        end

        q::Union{Missing, Float64} = reaction.q
        if ismissing(q)
            continue
        end

        # Reactant index
        z_r::Int, n_r::Int = reaction.product[1]
        reactant_idx::Int = zn_to_index(z_r, n_r, nd.net_idx)

        # Neutron (product) index
        neutron_idx::Int = zn_to_index(Int(0), Int(1), nd.net_idx)

        # Product index
        z_p::Int, n_p::Int = reaction.reactant[2]
        product_idx::Int = zn_to_index(z_p, n_p, nd.net_idx)

        # Lookup partition function for the reactant (product of the forward reaction)
        if !haskey(nd.reaction_data.neutroncapture, reactant_idx)
            pfunc_r::Float64 = 1.0
        else
            reactant = nd.reaction_data.neutroncapture[reactant_idx]
            pfunc_r = get_pfunc(reactant.rates_pfuncs_lerp, curr_traj.temperature)
        end
        pfunc_n::Float64 = 2.0
        pfunc_p::Float64 = get_pfunc(reaction.rates_pfuncs_lerp, curr_traj.temperature)

        pfunc = pfunc_n * pfunc_p / pfunc_r # FIXME: What about division by zero?
        if iszero(pfunc)
            continue
        end

        A_n::Int = 1
        A_p::Int = 0 + 1 + z_p + n_p
        A_r::Int = z_r + n_r

        # Reverse rate
        if q < 0
            reverse_rate::Float64 = 1e16
        else
            reverse_rate = forward_rate * 9.8678e9 * pfunc * (A_n*A_p/A_r)^(3/2) * curr_traj.temperature^(3/2) * exp(-11.605 * q / curr_traj.temperature)
        end

        # Reactant
        nd.jacobian[reactant_idx, reactant_idx] -= reverse_rate # J_RR

        # Products
        nd.jacobian[product_idx, reactant_idx]  += reverse_rate # J_PR
        nd.jacobian[neutron_idx, reactant_idx]  += reverse_rate # J_nR
    end
end

"""
    fill_jacobian!(nd::NetworkData; use_yproposed::Bool=false)::Vector{Float64}

Fill Jacobian with the data from the current time-step.

The flag `use_yproposed` is used to determine whether the network should use the abundances
held in `nd.abundance` or `nd.yproposed`. The latter should be used for an iteration that
failed to converge after a Newton-Raphson iteration.
"""
function fill_jacobian!(nd::NetworkData; use_yproposed::Bool=false)::Vector{Float64}
    # Jacobian coordinate: (reactant, product)

    mul!(nd.jacobian, nd.jacobian, 0.0)

    if nd.included_reactions.photodissociation
        fill_jacobian_photodissociation!(nd, use_yproposed)
    end
    if nd.included_reactions.ncap
        fill_jacobian_neutroncapture!(nd, use_yproposed)
    end
    if nd.included_reactions.alphadecay
        fill_jacobian_alphadecay!(nd, use_yproposed)
    end
    if nd.included_reactions.probdecay
        fill_jacobian_probdecay!(nd, use_yproposed)
    end

    mul!(nd.jacobian, nd.jacobian, -1.0)
    nd.jacobian[diagind(nd.jacobian)] .+= 1.0/nd.time.step
end

end