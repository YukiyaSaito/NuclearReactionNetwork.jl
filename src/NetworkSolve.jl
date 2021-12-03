"""
    NetworkSolve

This module is where the actual solving code for the network is held.
"""
module NetworkSolve
using ..Astro
using ..Network
using ..ReactionTypes
using ..InOut
using ..NetworkDatas
using ..LinearSolvers
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

export SolveNetwork!
export clip_abundance!

function check_mass_fraction_unity(nd::NetworkData, tolerance::Float64=1e-8)::Bool
    nd.yproposed[nd.yproposed .< 0.0] .= 0.0
    # y = copy(nd.yproposed)
    # y[y .< 0.0] .= 0.0 # Clip negative abundances
    return abs(1.0 - dot(nd.yproposed, nd.net_idx.mass_vector)) < tolerance
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

"""
    newton_raphson_iteration!(nd::NetworkData, Δy::Vector{Float64})::Int

Performs one round of Newton-Raphson.

If the mass fraction of the network fails to be within a certain tolerance, we keep
trying by first reducing the time step and doing more Newton-Raphson iterations.
"""
function newton_raphson_iteration!(nd::NetworkData, Δy::Vector{Float64})::Int
    nd.yproposed .= nd.abundance
    # lu_dot!(F,jacobian); 
    # ldiv!(Δy,F,(ydot.-(yproposed.-abundance) ./ nd.time.step))

    A::SparseMatrixCSC{Float64, Int} = nd.jacobian
    b::Vector{Float64} = ((nd.ydot .- (nd.yproposed .- nd.abundance) ./ nd.time.step))
    Δy::Vector{Float64} .= solve_linear_system(A, b, nd.solver)
    nd.yproposed .+= Δy

    num_failed::Int = 0
    num_tries::Int = 0
    while !check_mass_fraction_unity(nd)
        # @printf "\tnot converged; 1 - mass fraction: %e\n" abs(1 - dot(nd.yproposed, nd.net_idx.mass_vector))

        # Record how many times we've failed for later
        num_failed += 1

        # Halve the time step after 1 failed attempt
        num_tries += 1
        if num_tries > 1
            # @printf "\t\tHalving time step\n"
            nd.time.step /= 2.0
            num_tries = 0
        end

        # Update all the data
        fill_jacobian!(nd, use_yproposed=true)
        update_ydot!(nd, use_yproposed=true)

        A = nd.jacobian
        b = ((nd.ydot .- (nd.yproposed .- nd.abundance) ./ nd.time.step))
        # b = nd.ydot
        Δy .= solve_linear_system(A, b, nd.solver)
        nd.yproposed .+= Δy
    end
    nd.abundance .= nd.yproposed # FIXME: Performance boost: This could be just regular assignment (=, not .=) because of how NetworkDatas is setup
    step_time!(nd.time)
    return num_failed
end

function update_timestep_size!(nd::NetworkData, Δy::Vector{Float64})::Float64
    ∂y_∂t::Vector{Float64} = abs.(Δy[nd.abundance .> 1e-9] ./ nd.abundance[nd.abundance .> 1e-9])
    if iszero(maximum(∂y_∂t))
        nd.time.step *= 2.0
    else
        nd.time.step = min(2.0*nd.time.step, 0.25*nd.time.step/maximum(∂y_∂t))
    end
end

"""
    clip_abundance!(nd::NetworkData, min_val::Float64=0.0)::Vector{Float64}

Sets all abundances less than or equal to `min_val` to 0.
"""
function clip_abundance!(nd::NetworkData, min_val::Float64=0.0)::Vector{Float64}
    nd.abundance[nd.abundance .< min_val] .= 0.0
end

"""
    SolveNetwork!(nd::NetworkData)::Nothing

Solves the network from the data in `nd`.

We currently use Newton-Raphson to integrate ``\\vec{\\dot{Y}}``. This is the main
function of the network solver.
"""
function SolveNetwork!(nd::NetworkData)::Nothing
    Δy::Vector{Float64} = Vector{Float64}(undef, length(nd.abundance))
    # F = lu(nd.jacobian)
    # ps = MKLPardisoSolver()

    iteration::Int = 0
    failed_iterations::Int = 0
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
        # @printf "Time: %e,\tTime step: %e,\tIteration #: %d,\tFailed Iterations: %d,\t1 - Mass Fraction: %e\n" nd.time.current nd.time.step iteration failed_iterations abs(1 - dot(nd.yproposed, nd.net_idx.mass_vector))
        dump_iteration(nd, iteration)
        save_checkpoint(nd)
    end
end

end