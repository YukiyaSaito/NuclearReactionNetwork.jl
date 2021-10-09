module NetworkDatas

using ..Astro
using ..ReactionTypes
using ..Network
using SparseArrays
using LinearAlgebra
using DelimitedFiles

export NetworkData
export OutputInfo
export IncludedReactions
export fill_jacobian!
export update_ydot!

struct OutputInfo
    dump_final_y::Bool
    final_y_path::Union{Missing, String}
    dump_final_ya::Bool
    final_ya_path::Union{Missing, String}
    dump_each_iteration::Bool
    iteration_output_path::Union{Missing, String}
end

mutable struct IncludedReactions
    ncap::Bool
    probdecay::Bool
    alphadecay::Bool
    photodissociation::Bool
end

struct NetworkData
    net_idx::NetworkIndex
    reaction_data::ReactionData
    trajectory::Trajectory
    abundance::Vector{Float64}
    yproposed::Vector{Float64}
    ydot::Vector{Float64}
    time::Time
    jacobian::SparseMatrixCSC{Float64, Int}
    output_info::OutputInfo
    included_reactions::IncludedReactions
end

function initialize_ydot!(nd::NetworkData)::Vector{Float64}
    fill!(nd.ydot, 0.0)
end

function fill_probdecay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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

function fill_neutroncapture_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
    abundance::Vector{Float64} = nd.abundance
    if use_yproposed
        abundance = nd.yproposed
    end

    curr_traj::CurrentTrajectory = get_current_trajectory(nd.trajectory, nd.time.current)
    for capture::NeutronCapture in values(nd.reaction_data.neutroncapture)
        rate::Float64 = capture.rate(curr_traj.temperature)
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

function fill_alphadecay_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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

function fill_photodissociation_ydot!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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

        forward_rate::Float64 = reaction.rate(curr_traj.temperature)
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
            pfunc_r = nd.reaction_data.neutroncapture[reactant_idx].pfunc(curr_traj.temperature)
        end
        pfunc_n::Float64 = 2.0
        pfunc_p::Float64 = reaction.pfunc(curr_traj.temperature)

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

function fill_jacobian_probdecay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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

function fill_jacobian_neutroncapture!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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
        rate::Float64 = capture.rate(curr_traj.temperature)
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

function fill_jacobian_alphadecay!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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

function fill_jacobian_photodissociation!(nd::NetworkData, use_yproposed::Bool=false)::Nothing
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
        forward_rate::Float64 = reaction.rate(curr_traj.temperature)
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
            pfunc_r = nd.reaction_data.neutroncapture[reactant_idx].pfunc(curr_traj.temperature)
        end
        pfunc_n::Float64 = 2.0
        pfunc_p::Float64 = reaction.pfunc(curr_traj.temperature)

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

function fill_jacobian!(nd::NetworkData; use_yproposed::Bool=false)::Vector{Float64}
    # Jacobian coordinate: (reactant, product)

    fill!(nd.jacobian, 0.0)

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

    mul!(nd.jacobian, nd.jacobian, -1)
    nd.jacobian[diagind(nd.jacobian)] .+= 1/nd.time.step
end

end