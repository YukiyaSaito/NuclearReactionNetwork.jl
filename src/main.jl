using ..ReactionTypes
using ..Network
using ..NetworkSolve
using ..Astro
using ..InOut
using ..NetworkDatas
using LinearAlgebra
using SparseArrays
using Printf
using InteractiveUtils
using Profile
using Interpolations

export main
export read_test

function main()
    reactant = [[2,8]]
    product = [[3,7],[3,6],[0,1]]
    rate = 1.0
    average_number = [0.8, 0.2, 0.2]
    somedecay = ProbDecay(reactant, product, rate, average_number)
    println(somedecay)
end

function read_test(path::String)
    println("Setting BLAS thread count...")
    BLAS.set_num_threads(4)

    println("Loading data...")
    nd = initialize_network_data(path)
    try
        println("Solving network...")
        SolveNetwork!(nd)
    finally
        dump_result(nd)
        println("Wrote results to disk")
    end 
end