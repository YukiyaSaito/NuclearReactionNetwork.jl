module ReactionTypes
export ProbDecay
export ReactionData
export initialize_reactions
export read_probdecay!


abstract type AbstractReaction end

# mutable struct NBodyReaction <: AbstractReaction
#     reactant::Vector{Vector{Int32}}
#     product::Vector{Vector{Int32}}
# end


mutable struct ProbDecay <: AbstractReaction
    reactant::Vector{Vector{Int64}} #[[Z_0, N_0], [Z_1, N_1], ... [Z_n, N_n]]
    product::Vector{Vector{Int64}}
    rate::Float64
    average_number::Vector{Float64}
end


# ProbDecay(reactant::Vector{Vector{Int64}}, 
#             product::Vector{Vector{Int64}}, 
#             rate::Float64, 
#             average_number::Vector{Float64}) = ProbDecay(reactant, Set(product), rate, average_number)


# struct ReverseReaction <: AbstractReaction

# end


mutable struct ReactionData
    probdecay::Dict{Vector{Vector{Int64}},ProbDecay}
end

function initialize_reactions()
    probdecay_dict = Dict{Vector{Vector{Int64}},ProbDecay}()
    return ReactionData(probdecay_dict)
end



function read_probdecay!(path::String, reaction_data::ReactionData) ## Add Probdecay as a content of dict!
#reading line by line could be slow. Come back after implementing other parts. Use readlines() then split.
    open(path) do file
        infile = readlines(file)
        num_entry::Int64 = parse(Int, infile[1])
        # println(num_entry)
        # index = 1
        for i in 1:num_entry
            # 
            reactant_z::Vector{Int64}=[parse(Int, ss) for ss in split(infile[6*(i-1)+2])]
            reactant_n::Vector{Int64}=[parse(Int, ss) for ss in split(infile[6*(i-1)+3])]
            product_z::Vector{Int64}=[parse(Int, ss) for ss in split(infile[6*(i-1)+4])]
            product_n::Vector{Int64}=[parse(Int, ss) for ss in split(infile[6*(i-1)+5])]
            rate::Float64=parse(Float64, infile[6*(i-1)+6])
            average_number::Vector{Float64}=[parse(Float64, ss) for ss in split(infile[6*(i-1)+7])]
            reactant = Vector{Vector{Int64}}(undef,length(reactant_z))
            product = Vector{Vector{Int64}}(undef,length(product_z))
            for j in 1:length(reactant_z)
                reactant[j] = [reactant_z[j],reactant_n[j]]
            end
            for j in 1:length(product_z)
                product[j] = [product_z[j],product_n[j]]
            end
            reaction_data.probdecay[reactant] = ProbDecay(reactant,product,rate,average_number)
            # println(typeof(reactant_z[1]))
        #     # println(reactant_n)
        #     # println(product_z)
        #     # println(product_n)
        #     # println(rate)
        #     # println(average_number)
        end
        # println(reaction_data)
    end
    return reaction_data
end
    


end