module LinearSolvers

using SparseArrays
using Pardiso

export LinearSolver
export LS_UMFPACK
export LS_MKLPardisoSolver
export solve_linear_system

abstract type LinearSolver end

struct LS_UMFPACK <: LinearSolver
end

struct LS_MKLPardisoSolver <: LinearSolver
    pardiso::MKLPardisoSolver
end

function solve_linear_system(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, solver::LinearSolver)::Vector{Float64}
    println("Default")
    return A\b
end

function solve_linear_system(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, solver::LS_UMFPACK)::Vector{Float64}
    println("UMFPACK")
    return A\b
end

function solve_linear_system(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, solver::LS_MKLPardisoSolver)::Vector{Float64}
    println("Pardiso")
    return solve(solver.pardiso, A, b)
end


end