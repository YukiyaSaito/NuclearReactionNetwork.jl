module NetworkSolve
using ..Astro
using ..Network
using ..ReactionTypes
using ..InOut
using ...NetworkDatas
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

export check_mass_fraction_unity
export SolveNetwork!

function check_mass_fraction_unity(nd::NetworkData, tolerance::Float64=1e-8)
    return abs(1 - dot(nd.yproposed, nd.net_idx.mass_vector)) < tolerance
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

function newton_raphson_iteration!(nd::NetworkData, Δy::Vector{Float64})
    nd.yproposed .= nd.abundance
    # lu_dot!(F,jacobian); 
    # ldiv!(Δy,F,(ydot.-(yproposed.-abundance) ./ nd.time.step))
    Δy .= nd.jacobian \ ((nd.ydot .- (nd.yproposed .- nd.abundance) ./ nd.time.step))
    nd.yproposed .+= Δy

    num_failed::Int64 = 0
    num_tries::Int64 = 0
    while !check_mass_fraction_unity(nd)
    # while num_failed < 1
        @printf "\tnot converged; 1 - mass fraction: %e\n" abs(1 - dot(nd.yproposed, nd.net_idx.mass_vector))

        # Record how many times we've failed for later
        num_failed += 1

        # Halve the time step after 1 failed attempt
        num_tries += 1
        if num_tries > 1
            @printf "\t\tHalving time step\n"
            nd.time.step /= 2
            num_tries = 0
        end

        # Update all the data
        fill_jacobian!(nd, use_yproposed=true)
        update_ydot!(nd, use_yproposed=true)

        Δy .= nd.jacobian \ ((nd.ydot .- (nd.yproposed .- nd.abundance) ./ nd.time.step))
        nd.yproposed .+= Δy
    end
    nd.abundance .= nd.yproposed # FIXME: Performance boost: This could be just regular assignment (=, not .=) because of how NetworkDatas is setup
    nd.time.current += nd.time.step

    # Cap the time
    if nd.time.current >= nd.time.stop
        # nd.time.step = nd.time.stop - (nd.time.current - nd.time.step) TODO: Add this?
        nd.time.current = nd.time.stop
    end
    return num_failed
end

function update_timestep_size!(nd::NetworkData, Δy::Vector{Float64})
    ∂y_∂t::Vector{Float64} = abs.(Δy[nd.abundance .> 1e-9] ./ nd.abundance[nd.abundance .> 1e-9])
    if iszero(maximum(∂y_∂t))
        nd.time.step *= 2
    else
        nd.time.step = min(2*nd.time.step, 0.25*nd.time.step/maximum(∂y_∂t))
    end
end

function clip_abundance!(nd::NetworkData, min_val::Float64=0.0)
    nd.abundance[nd.abundance .< min_val] .= 0.0
end

function SolveNetwork!(nd::NetworkData)
    Δy = Vector{Float64}(undef, length(nd.abundance))
    # F = lu(nd.jacobian)
    # ps = MKLPardisoSolver()

    iteration = 0
    failed_iterations = 0
    @printf "Time: %e,\tTime step: %e,\tIteration #: %d,\tFailed Iterations: %d,\tAvg. Iterations/Timestep: %f\n" nd.time.current nd.time.step iteration failed_iterations (failed_iterations + iteration)/iteration
    dump_iteration(nd, iteration)
    while nd.time.current < nd.time.stop
        iteration += 1

        clip_abundance!(nd)
        failed_iterations += newton_raphson_iteration!(nd, Δy)
        update_timestep_size!(nd, Δy)
        fill_jacobian!(nd)
        update_ydot!(nd)
        
        @printf "Time: %e,\tTime step: %e,\tIteration #: %d,\tFailed Iterations: %d,\tAvg. Iterations/Timestep: %f\n" nd.time.current nd.time.step iteration failed_iterations (failed_iterations + iteration)/iteration
        dump_iteration(nd, iteration)
    end
end

end