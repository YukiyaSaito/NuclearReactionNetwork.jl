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

    n = Threads.nthreads()
    println("Num threads: $(n)")

    println("Loading data...")
    nd::NetworkData = initialize_network_data(path)
    rd = deepcopy(nd.rd)
    wds::Vector{RWData} = [deepcopy(nd.wd) for _ in 1:n]
    nds::Vector{NetworkData} = [NetworkData(rd, wd) for wd in wds]
    Threads.@threads for nd in nds
    # for nd in nds
        try
            println("Solving network...")
            SolveNetwork!(nd)
        finally
            clip_abundance!(nd,1e-15)
            dump_result(nd)
            println("Wrote results to disk")
        end
    end
end