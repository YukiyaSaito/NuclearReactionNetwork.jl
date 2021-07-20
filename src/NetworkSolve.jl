module NetworkSolve
using ..Network
using ..ReactionTypes
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
export initialize_timestep
export read_initial_abundance
export fill_jacobian!
export initialize_ydot
export fill_probdecay_ydot!
export newton_raphson_iteration!
export initialize_and_fill_sparse_jacobian
export check_mass_fraction_unity
export SolveNetwork!

# mutable struct timestep
#     current_value::Float64
# end

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

function initialize_timestep()
    #to do: summarize this in mutable struct
    timestep::Float64 = 1e-2
    current_time::Float64 = 0
    return timestep, current_time
end

function initialize_ydot(networksize::Int64)
    ydot::Vector{Float64} = zeros(Float64,networksize)
    return ydot
end

function initialize_ydot!(ydot::Vector{Float64})
    lmul!(0,ydot)    
end


function fill_probdecay_ydot!(ydot::Vector{Float64}, abundance::Vector{Float64}, reaction_data::ReactionData, zn_to_index_dict::Dict{Vector{Int64},Int64})
    for (key, value) in reaction_data.probdecay
        if haskey(zn_to_index_dict,value.reactant[1])
            ydot[zn_to_index_dict[value.reactant[1]]] += -1.0 * value.rate * abundance[zn_to_index_dict[value.reactant[1]]]
            # @printf "%d: %.10f += -1.0 * %f * %f\n" zn_to_index_dict[value.reactant[1]] ydot[zn_to_index_dict[value.reactant[1]]] value.rate abundance[zn_to_index_dict[value.reactant[1]]]
            # println("yes")
            # display(ydot[zn_to_index_dict[value.reactant[1]]])
        else 
            # throw(DomainError(value.reactant[1],"The species outside the reaction network."))
            continue
        end
        for product_num in 1:size(value.product)[1]
            if haskey(zn_to_index_dict, value.product[product_num])
                # println(value.reactant[1])
                # println(value.product[product_num])
                ydot[zn_to_index_dict[value.product[product_num]]] += value.average_number[product_num] * value.rate * abundance[zn_to_index_dict[value.reactant[1]]]
            elseif value.rate==0.0 || value.average_number[product_num]==0.0
                # throw(DomainError(value.product[product_num],"The species outside the reaction network."))
                continue
            else
                throw(DomainError(value.product[product_num],"The species outside the reaction network."))
            end
        end
        # println(value.product[1])
        # println(zn_to_index_dict[key[1]])
        # println(key)
        # println(value.rate)
        # println(value)
    end
end

function fill_initial_abundance!(abundance_index::Matrix{Int64},abundance_vector::Vector{Float64},abundance::Vector{Float64},zn_to_index_dict::Dict{Vector{Int64},Int64})
    neutron_num::Int64 = 0
    proton_num::Int64 = 0
    A_massnum::Float64 = 0.0
    for i in 1:size(abundance_index)[1]
        A_massnum = abundance_index[i,2]
        proton_num = abundance_index[i,1]
        neutron_num = abundance_index[i,2] - abundance_index[i,1]
        # println(abundance_index[i,:])
        if haskey(zn_to_index_dict, [proton_num,neutron_num])
            # println(sum(abundance_index[i,:]))
            abundance[zn_to_index_dict[[proton_num,neutron_num]]] = abundance_vector[i]/A_massnum
        else
            neutron_subtract::Int64 = 1
            # println(haskey(zn_to_index_dict, [abundance_index[i,1],abundance_index[i,2]-neutron_subtract]))
            while !haskey(zn_to_index_dict, [abundance_index[i,1],abundance_index[i,2]-neutron_subtract])
                # println(haskey(zn_to_index_dict, [abundance_index[i,1],abundance_index[i,2]-neutron_subtract]))
                neutron_subtract += 1
                # println(neutron_subtract)
            end				
            # println(typeof([abundance_index[i,1],abundance_index[i,2]-neutron_subtract]))
            
            abundance[zn_to_index_dict[[abundance_index[i,1],abundance_index[i,2]-neutron_subtract]]] += abundance_vector[i]/A_massnum
            @printf "Abundance for [%d, %d] added to [%d, %d]" abundance_index[i,1] abundance_index[i,2] abundance_index[i,1] abundance_index[i,2]-neutron_subtract
        end
        # @printf "[%d, %d]: %e\n" proton_num neutron_num abundance[zn_to_index_dict[[proton_num,neutron_num]]]
    end        
    return abundance
end

function read_initial_abundance(path::String, abundance::Vector{Float64},zn_to_index_dict::Dict{Vector{Int64},Int64})
    raw_abundance::Matrix{Float64} = readdlm(path)
    abundance_index::Matrix{Int64} = round.(Int,raw_abundance[:,1:2])
    abundance_vector::Vector{Float64} = raw_abundance[:,3]/sum(raw_abundance[:,3])
    # println(abundance_index)
    # println(size(abundance_index)[1])
    # for i in 1:size(abundance_index)[1]
    # 	println(typeof(abundance_index[i,:]))
    # end
    fill_initial_abundance!(abundance_index,abundance_vector,abundance,zn_to_index_dict)
end

function fill_jacobian!(jacobian::Matrix{Float64}, abundance::Vector{Float64},reaction_data::ReactionData, zn_to_index_dict::Dict{Vector{Int64},Int64},timestep::Float64)
    # jacobian = Matrix{Float64}(0*I,size(jacobian)) #Jacobian coordinate: (reactant, product)
    # if typeof(jacobian)==SparseMatrixCSC{Float64, Int64}
    #     mul!(jacobian,jacobian,0)
    # end    
    for (key, value) in reaction_data.probdecay
        if haskey(zn_to_index_dict,value.reactant[1])
            jacobian[zn_to_index_dict[value.reactant[1]], zn_to_index_dict[value.reactant[1]]] += -1.0 * value.rate
            # println("yes")
        else 
            # throw(DomainError(value.reactant[1],"The species outside the reaction network."))
            continue
        end
        for product_num in 1:size(value.product)[1]
            if haskey(zn_to_index_dict, value.product[product_num])
                # println(value.reactant[1])
                # println(value.product[product_num])
                jacobian[zn_to_index_dict[value.product[product_num]],zn_to_index_dict[value.reactant[1]]] += value.average_number[product_num] * value.rate
            elseif value.rate==0.0 || value.average_number[product_num]==0.0
                # throw(DomainError(value.product[product_num],"The species outside the reaction network."))
                continue
            else
                throw(DomainError(value.product[product_num],"The species outside the reaction network."))
            end
        end
        # println(value.product[1])
        # println(zn_to_index_dict[key[1]])
        # println(key)
        # println(value.rate)
        # println(value)
    end
    if typeof(jacobian)==Matrix{Float64}
        jacobian = sparse(jacobian)
    end
    # display(jacobian)
    mul!(jacobian,jacobian,-1)
    # display(jacobian)
    jacobian[diagind(jacobian)] .+= 1/timestep
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

function initialize_and_fill_sparse_jacobian(networksize::Int64, abundance::Vector{Float64},reaction_data::ReactionData, zn_to_index_dict::Dict{Vector{Int64},Int64},timestep::Float64)
    jacobian::Matrix{Float64} = initialize_jacobian(networksize)
    fill_jacobian!(jacobian, abundance, reaction_data, zn_to_index_dict, timestep)
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
    if abs(1-dot(yproposed,mass_vector))<1e-10
        return true
    else
        return false
    end
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


function newton_raphson_iteration!(abundance::Vector{Float64},yproposed::Vector{Float64},jacobian::SparseMatrixCSC{Float64, Int64},F::UmfpackLU, ydot::Vector{Float64}, ydelta::Vector{Float64} ,timestep::Float64, current_time::Float64, mass_vector::Vector{Float64})
    yproposed .= abundance
    # println("ok so far")
    # lu_dot!(F,jacobian); 
    # ldiv!(ydelta,F,ydot)
    # display(ydelta)
    # yproposed[:] .+= ydelta
    ydelta .= jacobian \ ydot
    yproposed .+= ydelta
    if check_mass_fraction_unity(yproposed, mass_vector) == true
        # ydelta .= yproposed .- abundance
        abundance .= yproposed
        current_time += timestep
    else 
        error("not converged")
    end

    # elseif converged == false
    #     while converged == false
    #         yproposed = jacobian \ (ydot - )


    # yproposed::Vector{Float64} = copy(abundance)
    # jacobian_inv::Matrix{Float64} = inv(jacobian)
    # yproposed = yproposed + jacobian_inv * ydot
    # yproposed::Vector{Float64} = jacobian \ ydot + yproposed
    return current_time
end

function update_timestep_size(abundance::Vector{Float64}, ydelta::Vector{Float64}, timestep::Float64)
    dely_delt::Vector{Float64} = ydelta[abundance.>1e-10]./abundance[abundance.>1e-10]
    # dely_delt[isnan.(dely_delt)] .= Inf
    # display(dely_delt)
    # display(min(2, 2*timestep, 0.1*timestep/maximum(dely_delt)))
    @printf "min(%e, %e, %e) = %e\n" 2.0 2*timestep 0.1*timestep/maximum(dely_delt) min(2.0, 2*timestep, 0.1*timestep/maximum(dely_delt))
    return min(2.0, 2*timestep, 0.1*timestep/maximum(dely_delt))
end

function newton_raphson_iteration!(abundance::Vector{Float64},yproposed::Vector{Float64},jacobian::SparseMatrixCSC{Float64, Int64},ps::MKLPardisoSolver, ydot::Vector{Float64}, ydelta::Vector{Float64} ,timestep::Float64, current_time::Float64, mass_vector::Vector{Float64})
    yproposed .= abundance
    # println("ok so far")
    solve!(ps, ydelta, jacobian, ydot)
    yproposed .+= ydelta
    # yproposed += jacobian \ ydot

    if check_mass_fraction_unity(yproposed, mass_vector) == true
        # ydelta .= yproposed .- abundance
        abundance .= yproposed
        # current_time += timestep
    else 

        error("not converged")
    end

    # elseif converged == false
    #     while converged == false
    #         yproposed = jacobian \ (ydot - )


    # yproposed::Vector{Float64} = copy(abundance)
    # jacobian_inv::Matrix{Float64} = inv(jacobian)
    # yproposed = yproposed + jacobian_inv * ydot
    # yproposed::Vector{Float64} = jacobian \ ydot + yproposed
    return current_time += timestep
end

function SolveNetwork!(abundance::Vector{Float64},jacobian::SparseMatrixCSC{Float64, Int64},reaction_data::ReactionData,ydot::Vector{Float64} ,timestep::Float64, current_time::Float64, mass_vector::Vector{Float64},time_limit::Float64,zn_to_index_dict::Dict{Vector{Int64},Int64})
    ydelta = Vector{Float64}(undef,size(abundance)[1])
    yproposed = Vector{Float64}(undef,size(abundance)[1])
    F = lu(jacobian)
    print_time_step::Float64 = 0.1
    # ps = MKLPardisoSolver()
    while current_time < time_limit
        current_time = newton_raphson_iteration!(abundance,yproposed,jacobian, F ,ydot, ydelta,timestep, current_time, mass_vector)
        # display(ydot)
        initialize_ydot!(ydot)
        fill_probdecay_ydot!(ydot,abundance,reaction_data,zn_to_index_dict)
        # display(ydot)
        # timestep = update_timestep_size(abundance, ydelta, timestep)
        if current_time > print_time_step
            # println(ydot)
            @printf "[%f, %e],\n" current_time abundance[zn_to_index_dict[[54,93]]]
            print_time_step += 0.1
        end
        # @printf "%e\n" current_time
    end
end

end