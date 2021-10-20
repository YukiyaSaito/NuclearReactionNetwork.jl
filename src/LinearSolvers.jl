"""
Provides an interface to various linear solvers.

Currently, only two solvers are supported: MKLPardiso and UMFPACK.
"""
module LinearSolvers

using SparseArrays
using Pardiso

export LinearSolver
export LS_UMFPACK
export LS_MKLPardisoSolver
export solve_linear_system

"""A generic type used to refer to any linear solver."""
abstract type LinearSolver end

"""
    LS_UMFPACK <: LinearSolver

A type to interface the UMFPACK linear solver.
"""
struct LS_UMFPACK <: LinearSolver
end

"""
    LS_MKLPardisoSolver <: LinearSolver

A type to interface the MKLParidso linear solver.

# Fields:
- `pardiso::MKLPardisoSolver`: The PARDISO solver.
"""
struct LS_MKLPardisoSolver <: LinearSolver
    """The PARDISO solver."""
    pardiso::MKLPardisoSolver
end

"""
    solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LinearSolver)

Solves the system ``\\mathcal{A}\\vec{x} = \\vec{b}`` using [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html).

This method exists so that if no solver is provided, the default is UMFPACK.
"""
function solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LinearSolver)::Vector{Float64}
    return A\b
end

"""
    solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LS_UMFPACK)

Solves the system ``\\mathcal{A}\\vec{x} = \\vec{b}`` using [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html).
"""
function solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LS_UMFPACK)::Vector{Float64}
    return A\b
end

"""
    solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LS_MKLPardisoSolver)

Solves the system ``\\mathcal{A}\\vec{x} = \\vec{b}`` using [PARDISO](https://www.pardiso-project.org/).
"""
function solve_linear_system(A::SparseMatrixCSC{Float64, Int}, b::Vector{Float64}, solver::LS_MKLPardisoSolver)::Vector{Float64}
    return solve(solver.pardiso, A, b)
end

end