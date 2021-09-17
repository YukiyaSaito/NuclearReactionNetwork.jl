module NetworkSolve
using ..Astro
using ..Network
using ..ReactionTypes
using ..InOut
using LinearAlgebra
using DelimitedFiles
using Printf
using SparseArrays
using InteractiveUtils
using SuiteSparse.UMFPACK
using SparseArrays: getcolptr
using SuiteSparse: decrement
using SuiteSparse.UMFPACK: UmfpackLU, UMFVTypes, UMFITypes, umfpack_numeric!
using Pardiso

export initialize_jacobian
export initialize_abundance
export read_initial_abundance
export fill_jacobian!
export initialize_ydot
export fill_probdecay_ydot!
export newton_raphson_iteration!
export initialize_and_fill_sparse_jacobian
export check_mass_fraction_unity
export SolveNetwork!

function initialize_jacobian(networksize::Int64)
    # println(networksize)
    jacobian::Matrix{Float64} = zeros(Float64,networksize,networksize)
    return jacobian
    # println(size(jacobian))
end

function initialize_abundance(networksize::Int64)
    abundance::Vector{Float64} = zeros(Float64,networksize)
    return abundance
end

function initialize_ydot(networksize::Int64)
    ydot::Vector{Float64} = zeros(Float64,networksize)
    return ydot
end

function initialize_ydot!(ydot::Vector{Float64})
    lmul!(0,ydot)    
end


function fill_probdecay_ydot!(ydot::Vector{Float64}, abundance::Vector{Float64}, reaction_data::ReactionData, net_idx::NetworkIndex)
    for (_, decay) in reaction_data.probdecay
        if zn_in_network(decay.reactant[1]..., net_idx)
            reactant_idx = zn_to_index(decay.reactant[1]..., net_idx)
            ydot[reactant_idx] += -1.0 * decay.rate * abundance[reactant_idx]
            # @printf "%d: %.10f += -1.0 * %f * %f\n" zn_to_index_dict[value.reactant[1]] ydot[zn_to_index_dict[value.reactant[1]]] value.rate abundance[zn_to_index_dict[value.reactant[1]]]
            # println("yes")
            # display(ydot[zn_to_index_dict[value.reactant[1]]])
        else 
            # throw(DomainError(value.reactant[1],"The species outside the reaction network."))
            continue
        end
        for product_num in 1:size(decay.product)[1]
            if zn_in_network(decay.product[product_num]..., net_idx)
                # println(value.reactant[1])
                # println(value.product[product_num])
                ydot[zn_to_index(decay.product[product_num]..., net_idx)] += decay.average_number[product_num] * decay.rate * abundance[zn_to_index(decay.reactant[1]..., net_idx)]
            elseif decay.rate==0.0 || decay.average_number[product_num]==0.0
                # throw(DomainError(value.product[product_num],"The species outside the reaction network."))
                continue
            else
                throw(DomainError(decay.product[product_num],"The species outside the reaction network."))
            end
        end
        # println(value.product[1])
        # println(zn_to_index_dict[key[1]])
        # println(key)
        # println(value.rate)
        # println(value)
    end
end

function fill_neutroncapture_ydot!(ydot::Vector{Float64}, abundance::Vector{Float64}, reaction_data::ReactionData, net_idx::NetworkIndex, trajectory::Trajectory, time::Time)
    curr_traj = get_current_trajectory(trajectory, time.current)
    for capture in values(reaction_data.neutroncapture)
        rate = capture.rate(curr_traj.temperature)
        if iszero(rate)
            println("Zero rate in fill_neutroncapture_ydot!")
            continue
        end
        # Grab the product of all the abundances
        abundance_factor = 1.0
        for reactant in capture.reactant
            z, n = reactant
            if zn_in_network(z, n, net_idx)
                reactant_idx = zn_to_index(z, n, net_idx)
                abundance_factor *= abundance[reactant_idx]
            end
        end
        if iszero(abundance_factor)
            continue
        end

        # Update ydot for the reactants
        for reactant in capture.reactant
            z, n = reactant
            if zn_in_network(z, n, net_idx)
                reactant_idx = zn_to_index(z, n, net_idx)
                ydot[reactant_idx] += -1.0 * curr_traj.density * rate * abundance_factor
            end
        end

        # Update ydot for the products
        for product in capture.product
            z, n = product
            if zn_in_network(z, n, net_idx)
                product_idx = zn_to_index(z, n, net_idx)
                ydot[product_idx] += rate * curr_traj.density * abundance_factor
            end
        end
    end
end

function fill_alphadecay_ydot!(ydot::Vector{Float64}, abundance::Vector{Float64}, reaction_data::ReactionData, net_idx::NetworkIndex)
    for decay in reaction_data.alphadecay
        rate = decay.rate
        if iszero(rate)
            continue
        end

        # Grab the abundace of the reactant
        z_r, n_r = decay.reactant
        if !zn_in_network(z_r, n_r, net_idx)
            continue
        end
        reactant_idx = zn_to_index(z_r, n_r, net_idx)
        abundance_factor = abundance[reactant_idx]
        if iszero(abundance_factor)
            continue
        end

        # Update ydot for the reactant
        ydot[reactant_idx] += -1.0 * rate * abundance_factor

        # Update ydot for the products
        for product in decay.product
            z_p, n_p = product
            if !zn_in_network(z_p, n_p, net_idx)
                continue
            end
            product_idx = zn_to_index(z_p, n_p, net_idx)
            ydot[product_idx] += 1.0 * rate * abundance_factor
        end
    end
end

function update_ydot!(ydot::Vector{Float64}, abundance::Vector{Float64}, reaction_data::ReactionData, net_idx::NetworkIndex, trajectory::Trajectory, time::Time)
    initialize_ydot!(ydot)
    fill_probdecay_ydot!(ydot, abundance, reaction_data, net_idx)
    fill_neutroncapture_ydot!(ydot, abundance, reaction_data, net_idx, trajectory, time)
    fill_alphadecay_ydot!(ydot, abundance, reaction_data, net_idx)
end

function fill_initial_abundance!(abundance_index::Matrix{Int64}, abundance_vector::Vector{Float64}, abundance::Vector{Float64}, net_idx::NetworkIndex)
    neutron_num::Int64 = 0
    proton_num::Int64 = 0
    A_massnum::Float64 = 0.0
    for i in 1:size(abundance_index)[1]
        A_massnum = abundance_index[i,2]
        proton_num = abundance_index[i,1]
        neutron_num = abundance_index[i,2] - abundance_index[i,1]
        # println(abundance_index[i,:])
        if zn_in_network(proton_num, neutron_num, net_idx)
            # println(sum(abundance_index[i,:]))
            abundance[zn_to_index(proton_num, neutron_num, net_idx)] = abundance_vector[i]/A_massnum
        else
            neutron_subtract::Int64 = 1
            # println(haskey(zn_to_index_dict, [abundance_index[i,1],abundance_index[i,2]-neutron_subtract]))
            while !zn_in_network(abundance_index[i, 1], abundance_index[i, 2] - neutron_subtract, net_idx)
                # println(haskey(zn_to_index_dict, [abundance_index[i,1],abundance_index[i,2]-neutron_subtract]))
                neutron_subtract += 1
                # println(neutron_subtract)
            end				
            # println(typeof([abundance_index[i,1],abundance_index[i,2]-neutron_subtract]))
            
            abundance[zn_to_index(abundance_index[i, 1], abundance_index[i, 2] - neutron_subtract, net_idx)] += abundance_vector[i]/A_massnum
            @printf "Abundance for [%d, %d] added to [%d, %d]" abundance_index[i,1] abundance_index[i,2] abundance_index[i,1] abundance_index[i,2]-neutron_subtract
        end
        # @printf "[%d, %d]: %e\n" proton_num neutron_num abundance[zn_to_index_dict[[proton_num,neutron_num]]]
    end        
    return abundance
end

function read_initial_abundance(path::String, abundance::Vector{Float64}, net_idx::NetworkIndex)
    raw_abundance::Matrix{Float64} = readdlm(path)
    abundance_index::Matrix{Int64} = round.(Int,raw_abundance[:,1:2])
    abundance_vector::Vector{Float64} = raw_abundance[:,3]/sum(raw_abundance[:,3])
    # println(abundance_index)
    # println(size(abundance_index)[1])
    # for i in 1:size(abundance_index)[1]
    # 	println(typeof(abundance_index[i,:]))
    # end
    fill_initial_abundance!(abundance_index, abundance_vector, abundance, net_idx)
end

function fill_jacobian_probdecay!(jacobian::Union{Matrix{Float64},SparseMatrixCSC{Float64, Int64}}, reaction_data::ReactionData, net_idx::NetworkIndex)
    for (_, decay) in reaction_data.probdecay
        z_r, n_r = decay.reactant[1]
        if zn_in_network(z_r, n_r, net_idx)
            reactant_idx = zn_to_index(z_r, n_r, net_idx)
            jacobian[reactant_idx, reactant_idx] += -1.0 * decay.rate
            # println("yes")
        else 
            # throw(DomainError(value.reactant[1],"The species outside the reaction network."))
            continue
        end
        for product_num in 1:size(decay.product)[1]
            if zn_in_network(decay.product[product_num]..., net_idx)
                # println(value.reactant[1])
                # println(value.product[product_num])
                z_p, n_p = decay.product[product_num]
                product_idx = zn_to_index(z_p, n_p, net_idx)
                jacobian[product_idx, reactant_idx] += decay.average_number[product_num] * decay.rate
            elseif decay.rate==0.0 || decay.average_number[product_num]==0.0
                # throw(DomainError(value.product[product_num],"The species outside the reaction network."))
                continue
            else
                throw(DomainError(decay.product[product_num],"The species outside the reaction network."))
            end
        end
    end
end

function fill_jacobian_neutroncapture!(jacobian::Union{Matrix{Float64},SparseMatrixCSC{Float64, Int64}}, abundance::Vector{Float64}, reaction_data::ReactionData, trajectory::Trajectory, net_idx::NetworkIndex, time::Time)
    curr_traj = get_current_trajectory(trajectory, time.current)
    if iszero(curr_traj.density)
        println("Zero density")
        return
    end

    # TODO: Convert this to a loop instead of 6 hard coded changes to the jacobian?
    for capture in values(reaction_data.neutroncapture)
        rate = capture.rate(curr_traj.temperature)
        if iszero(rate)
            println("Zero rate in fill_jacobian_neutroncapture!")
            continue
        end

        # Neutron index and abundance
        if !zn_in_network(0, 1, net_idx) # TODO: Are neutrons always in the network making this check redundant?
            continue
        end
        neutron_idx = zn_to_index(0, 1, net_idx)
        neutron_abundance = abundance[neutron_idx]

        # Reactant index and abundance
        reactant = capture.reactant[2]
        z_r, n_r = reactant
        if !zn_in_network(z_r, n_r, net_idx)
            continue
        end
        reactant_idx = zn_to_index(z_r, n_r, net_idx)
        reactant_abundance = abundance[reactant_idx]

        # Product index
        product = capture.product[1]
        z_p, n_p = product
        if !zn_in_network(z_p, n_p, net_idx)
            continue
        end
        product_idx = zn_to_index(z_p, n_p, net_idx)

        # Doule counting factor TODO: Is this always 1.0?
        dc_factor = reactant_idx == product_idx ? 0.5 : 1.0

        # Reactants
        jacobian[reactant_idx, reactant_idx] -= dc_factor * curr_traj.density * rate * neutron_abundance  # J_RR
        jacobian[reactant_idx, neutron_idx]  -= dc_factor * curr_traj.density * rate * reactant_abundance # J_Rn
        jacobian[neutron_idx, reactant_idx]  -= dc_factor * curr_traj.density * rate * neutron_abundance  # J_nR
        jacobian[neutron_idx, neutron_idx]   -= dc_factor * curr_traj.density * rate * reactant_abundance # J_nn

        # Products
        jacobian[product_idx, reactant_idx]  += dc_factor * curr_traj.density * rate * neutron_abundance  # J_PR
        jacobian[product_idx, neutron_idx]   += dc_factor * curr_traj.density * rate * reactant_abundance # J_Pn
    end
end

function fill_jacobian_alphadecay!(jacobian::Union{Matrix{Float64},SparseMatrixCSC{Float64, Int64}}, reaction_data::ReactionData, net_idx::NetworkIndex)
    for decay in reaction_data.alphadecay
        z_r, n_r = decay.reactant
        if !zn_in_network(z_r, n_r, net_idx)
            continue
        end
        reactant_idx = zn_to_index(z_r, n_r, net_idx)
        jacobian[reactant_idx, reactant_idx] += -1.0 * decay.rate

        for product in decay.product
            z_p, n_p = product
            if !zn_in_network(z_p, n_p, net_idx)
                continue
            end
            product_idx = zn_to_index(z_p, n_p, net_idx)
            jacobian[product_idx, reactant_idx] += decay.rate
        end
    end
end

function fill_jacobian!(jacobian::Union{Matrix{Float64},SparseMatrixCSC{Float64, Int64}}, abundance::Vector{Float64}, reaction_data::ReactionData, trajectory::Trajectory, net_idx::NetworkIndex, time::Time)
    # jacobian = Matrix{Float64}(0*I,size(jacobian)) #Jacobian coordinate: (reactant, product)
    if typeof(jacobian)==SparseMatrixCSC{Float64, Int64}
        mul!(jacobian,jacobian,0)
    end

    fill_jacobian_probdecay!(jacobian, reaction_data, net_idx)
    fill_jacobian_neutroncapture!(jacobian, abundance, reaction_data, trajectory, net_idx, time)
    fill_jacobian_alphadecay!(jacobian, reaction_data, net_idx)

    if typeof(jacobian)==Matrix{Float64}
        jacobian = sparse(jacobian)
    end
    # display(jacobian)
    mul!(jacobian,jacobian,-1)
    # display(jacobian)
    jacobian[diagind(jacobian)] .+= 1/time.step
    # jacobian .+= Matrix((1.0/timestep)I, size(jacobian))
    # @printf "%.5f\n" jacobian[1,1]
    # display(jacobian)
    # jacobian = jacobian*-1
    # jacobian .+= Matrix{Float64}(I, size(jacobian)) / timestep
    # num_entry::Int64 = 0
    # num_nonzero_entry::Int64 = 0
    # display(jacobian)
    # display(typeof(jacobian)==SparseMatrixCSC{Float64, Int64})
    # for i in 1:size(jacobian)[1]
    #     for j in 1:size(jacobian)[1]
    #         # num_entry += 1
    #         if jacobian[j,i] != 0.0 && j==i && (j<100 || i<100)
    #             @printf "[%d, %d]: %f\n" j i jacobian[j,i]
    #         # if jacobian[j,i] != 0.0
    #         #     num_nonzero_entry+=1
    #         end
    #     end
    # end
    # @printf "%f, %f\n" num_nonzero_entry num_entry
    return jacobian
end

function initialize_and_fill_sparse_jacobian(networksize::Int64, abundance::Vector{Float64}, reaction_data::ReactionData, trajectory::Trajectory, net_idx::NetworkIndex, time::Time)
    jacobian::Matrix{Float64} = initialize_jacobian(networksize)
    fill_jacobian!(jacobian, abundance, reaction_data, trajectory, net_idx, time)
end

function newton_raphson_step(yproposed::Vector{Float64},jacobian::Matrix{Float64},jacobian_inv::Matrix{Float64},ydot::Vector{Float64})
    jacobian_inv = inv(jacobian) 
    yproposed = yproposed + jacobian_inv * ydot 
    return yproposed  
end

function check_mass_fraction_unity(yproposed::Vector{Float64},mass_vector::Vector{Float64})
    # mass_fraction_sum::Float64 = 0
    # mass_fraction_sum = dot(yproposed,mass_vector)
    # display(mass_fraction_sum)
    # println(abs(1-dot(yproposed,mass_vector)))
    return abs(1 - dot(yproposed, mass_vector)) < 1e-10
end

function lu_dot!(F::UmfpackLU, S::SparseMatrixCSC{<:UMFVTypes,<:UMFITypes}; check::Bool=true)
    zerobased = getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)
    if zerobased
        F.colptr .= getcolptr(S)
        F.rowval .= rowvals(S)
    else
        F.colptr .= getcolptr(S) .- oneunit(eltype(S))
        F.rowval .= rowvals(S) .- oneunit(eltype(S))
    end
    F.nzval .= nonzeros(S)

    umfpack_numeric!(F, reuse_numeric = false)
    check && (issuccess(F) || throw(LinearAlgebra.SingularException(0)))
    return F
end


function newton_raphson_iteration!(abundance::Vector{Float64}, yproposed::Vector{Float64}, jacobian::SparseMatrixCSC{Float64, Int64}, F::UmfpackLU, ydot::Vector{Float64}, ydelta::Vector{Float64}, time::Time, reaction_data::ReactionData, net_idx::NetworkIndex, trajectory::Trajectory)
    yproposed .= abundance
    # println("ok so far")
    # lu_dot!(F,jacobian); 
    # ldiv!(ydelta,F,(ydot.-(yproposed.-abundance)./time.step))
    # display(ydelta)
    # @printf "%e\n" timestep
    # yproposed[:] .+= ydelta
    # display(jacobian[1,1])
    # display(ydelta)
    # ydelta .= jacobian \ ydot
    ydelta .= jacobian \ ((ydot .- (yproposed .- abundance) ./ time.step))
    yproposed .+= ydelta
    if check_mass_fraction_unity(yproposed, net_idx.mass_vector)
        # display(true)
        # ydelta .= yproposed .- abundance
    else 
        while !check_mass_fraction_unity(yproposed, net_idx.mass_vector) # TODO: Convert to a do-while style loop to avoid unnecessary computations
            @printf "not converged; 1 - mass fraction: %e\n" abs(1-dot(yproposed, net_idx.mass_vector))
            fill_jacobian!(jacobian, yproposed, reaction_data, trajectory, net_idx, time)
            update_ydot!(ydot, yproposed, reaction_data, net_idx, trajectory, time)

            ydelta .= jacobian \ ((ydot .- (yproposed .- abundance) ./ time.step))
            yproposed .+= ydelta
            # error("not converged")
            # break
        end
    end
    abundance .= yproposed
    time.current += time.step

    # Cap the time
    if time.current >= time.stop
        time.current = time.stop
    end
    # elseif converged == false
    #     while converged == false
    #         yproposed = jacobian \ (ydot - )


    # yproposed::Vector{Float64} = copy(abundance)
    # jacobian_inv::Matrix{Float64} = inv(jacobian)
    # yproposed = yproposed + jacobian_inv * ydot
    # yproposed::Vector{Float64} = jacobian \ ydot + yproposed
end

function update_timestep_size!(abundance::Vector{Float64}, ydelta::Vector{Float64}, time::Time)
    dely_delt::Vector{Float64} = abs.(ydelta[abundance.>1e-9]./abundance[abundance.>1e-9])
    if maximum(dely_delt)==0
        time.step *= 2
    else
        time.step = min(2*time.step, 0.25*time.step/maximum(dely_delt))
    end
    # dely_delt[isnan.(dely_delt)] .= Inf
    # display(ydelta[abundance.>1e-25])
    # display(min(2, 2*timestep, 0.1*timestep/maximum(dely_delt)))
    # @printf "min(%e, %e, %e) = %e\n" 0.1 2*timestep 0.1*timestep/maximum(dely_delt) min(0.1, 2*timestep, 0.1*timestep/maximum(dely_delt))
    # return min(1e-4,1e-4,1e-4)
end

# function newton_raphson_iteration!(abundance::Vector{Float64},yproposed::Vector{Float64},jacobian::SparseMatrixCSC{Float64, Int64},ps::MKLPardisoSolver, ydot::Vector{Float64}, ydelta::Vector{Float64} ,timestep::Float64, current_time::Float64, mass_vector::Vector{Float64})
#     yproposed .= abundance
#     # println("ok so far")
#     solve!(ps, ydelta, jacobian, ydot)
#     yproposed .+= ydelta
#     # yproposed += jacobian \ ydot

#     if check_mass_fraction_unity(yproposed, mass_vector) == true
#         # ydelta .= yproposed .- abundance
#         abundance .= yproposed
#         # current_time += timestep
#     else 

#         error("not converged")
#     end

#     # elseif converged == false
#     #     while converged == false
#     #         yproposed = jacobian \ (ydot - )


#     # yproposed::Vector{Float64} = copy(abundance)
#     # jacobian_inv::Matrix{Float64} = inv(jacobian)
#     # yproposed = yproposed + jacobian_inv * ydot
#     # yproposed::Vector{Float64} = jacobian \ ydot + yproposed
#     return current_time += timestep
# end

function SolveNetwork!(abundance::Vector{Float64}, jacobian::SparseMatrixCSC{Float64, Int64}, reaction_data::ReactionData, ydot::Vector{Float64}, time::Time, net_idx::NetworkIndex, trajectory::Trajectory, dump_ytime::Bool=false)
    ydelta = Vector{Float64}(undef,size(abundance)[1])
    yproposed = Vector{Float64}(undef,size(abundance)[1])
    F = lu(jacobian)
    print_time_step::Float64 = 10.0
    # ps = MKLPardisoSolver()
    iteration = 1
    while time.current < time.stop
        newton_raphson_iteration!(abundance, yproposed, jacobian, F, ydot, ydelta, time, reaction_data, net_idx, trajectory)
        # display(current_time)
        # display(ydot)
        update_timestep_size!(abundance, ydelta, time)
        fill_jacobian!(jacobian, abundance, reaction_data, trajectory, net_idx, time)
        # println(current_time)
        update_ydot!(ydot, abundance, reaction_data, net_idx, trajectory, time)
        # display(ydot)
        @printf "[%e, %e],\n" time.current abundance[zn_to_index(0, 1, net_idx)]

        if dump_ytime
            result = Result(abundance, net_idx)
            dump_iteration(result, trajectory, time, iteration, "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/output/wind-beta+ncap+alpha-YTime.txt")
        end

        # if current_time > print_time_step
            # println(ydot)
            # @printf "[%f, %e],\n" current_time abundance[zn_to_index_dict[[0,1]]]
            # print_time_step += 10.0
        # end
        # @printf "%e\n" current_time
        iteration += 1
    end
end

end