module ReactionTypes

using Interpolations
using ..Astro
using ..Network

export ProbDecay
export ReactionData
export initialize_reactions
export read_probdecay!
export read_ncap!
export read_alphadecay!
export read_photodissociation!

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

mutable struct NeutronCapture <: AbstractReaction #Temperature Dependent Reactions
    reactant::Vector{Vector{Int64}} #[[Z_0, N_0], [Z_1, N_1], ... [Z_n, N_n]]
    product::Vector{Vector{Int64}}
    rate::Interpolations.Extrapolation
    pfunc::Interpolations.Extrapolation
    current_rate::Float64 
    q::Union{Missing, Float64}
end

# TODO: Make every decay struct have the same structures (i.e. tuples, or vectors, or something els)
mutable struct AlphaDecay <: AbstractReaction
    reactant::Tuple{Int64, Int64} # (Z, N)
    product::Vector{Tuple{Int64, Int64}} # [(Z_0, N_0), (Z_1, N_1)]
    rate::Float64
end

# ProbDecay(reactant::Vector{Vector{Int64}}, 
#             product::Vector{Vector{Int64}}, 
#             rate::Float64, 
#             average_number::Vector{Float64}) = ProbDecay(reactant, Set(product), rate, average_number)


# struct ReverseReaction <: AbstractReaction

# end


# TODO: Replace all dicts to simple vectors
mutable struct ReactionData
    probdecay::Dict{Vector{Vector{Int64}},ProbDecay}
    neutroncapture::Dict{Vector{Vector{Int64}},NeutronCapture}
    alphadecay::Vector{AlphaDecay}
end

function initialize_reactions()
    probdecay_dict = Dict{Vector{Vector{Int64}},ProbDecay}()
    ncap_dict = Dict{Vector{Vector{Int64}},NeutronCapture}()
    alphadecay = Vector{AlphaDecay}()
    return ReactionData(probdecay_dict, ncap_dict, alphadecay)
end



function read_probdecay!(path::String, reaction_data::ReactionData, net_idx::NetworkIndex) 
#Read Fortran formatted PRISM input file for probabilistic decay.
#reading line by line could be slow. Come back after implementing other parts.
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

            # Check if any of the products or reactants are outside of the network
            out_of_network = false
            for (z, n) in [reactant; product]
                if !zn_in_network(z, n, net_idx)
                    out_of_network = true
                    break
                end
            end
            # If any of the species are outside of the network, skip this reaction
            if out_of_network
                continue
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
    
function read_ncap!(path::String, reaction_data::ReactionData, net_idx::NetworkIndex) 
    #Read Fortran formatted input file for temperature dependent reaction rate. 
    #reading line by line could be slow. Come back after implementing other parts.
    open(path) do file
        infile = readlines(file)
        num_entry::Int64 = parse(Int, infile[1])
        pfunc_flag::Int64 = parse(Int, infile[2])
        temperature::Vector{Float64}=[parse(Float64, ss) for ss in split(infile[3])]
        reactant_z_temp = zeros(Int64,2)
        reactant_n_temp = zeros(Int64,2)
        product_z_temp = zeros(Int64,1)
        product_n_temp = zeros(Int64,1)
        rate_temp = zeros(Float64,length(temperature))
        pfunc_temp = zeros(Float64,length(temperature))
        reactant_temp = Vector{Vector{Int64}}(undef,length(reactant_z_temp))
        product_temp = Vector{Vector{Int64}}(undef,length(product_z_temp))
        current_rate::Float64 = 0.0
        if pfunc_flag==1
            for i in 1:num_entry
                reactant_z_temp.=[parse(Int, ss) for ss in split(infile[6*(i-1)+4])]
                reactant_n_temp.=[parse(Int, ss) for ss in split(infile[6*(i-1)+5])]
                product_z_temp.=[parse(Int, ss) for ss in split(infile[6*(i-1)+6])]
                product_n_temp.=[parse(Int, ss) for ss in split(infile[6*(i-1)+7])]
                rate_temp.=[parse(Float64, ss) for ss in split(infile[6*(i-1)+8])]
                pfunc_temp.=[parse(Float64, ss) for ss in split(infile[6*(i-1)+9])]

                for j in 1:length(reactant_z_temp)
                    reactant_temp[j] = [reactant_z_temp[j],reactant_n_temp[j]]
                end
                for j in 1:length(product_z_temp)
                    product_temp[j] = [product_z_temp[j],product_n_temp[j]]
                end

                # Check if any of the products or reactants are outside of the network
                out_of_network = false
                for (z, n) in [reactant_temp; product_temp]
                    if !zn_in_network(z, n, net_idx)
                        out_of_network = true
                        break
                    end
                end
                # If any of the species are outside of the network, skip this reaction
                if out_of_network
                    continue
                end

                rate_lerp = LinearInterpolation(copy(temperature), copy(rate_temp), extrapolation_bc=Flat())
                pfunc_lerp = LinearInterpolation(copy(temperature), copy(pfunc_temp), extrapolation_bc=Flat())
                reaction_data.neutroncapture[copy(reactant_temp)] = NeutronCapture(copy(reactant_temp), copy(product_temp), rate_lerp, pfunc_lerp, copy(current_rate), missing)
    # 			#     # println(reactant_n)
    #             #     # println(product_z)
    #             #     # println(product_n)
    #             #     # println(rate)
    #             #     # println(average_number)
                # println(Set(reactant_temp))
            end
        end
    end
    return reaction_data
end

function read_alphadecay!(path::String, reaction_data::ReactionData, net_idx::NetworkIndex)
    #Read Fortran formatted input file for temperature dependent reaction rate. 
    #reading line by line could be slow. Come back after implementing other parts.
    open(path) do file
        lines = readlines(file)
        num_entries::Int64 = parse(Int, lines[1])
        for i in 1:num_entries
            z_r = parse(Int64, lines[5*(i-1) + 2])
            n_r = parse(Int64, lines[5*(i-1) + 3])
            z_p = [parse(Int64, s) for s in split(lines[5*(i-1) + 4])]
            n_p = [parse(Int64, s) for s in split(lines[5*(i-1) + 5])]
            rate = parse(Float64, lines[5*(i-1) + 6])

            reactant = (z_r, n_r)
            product = [(z_p[1], n_p[1]), (z_p[2], n_p[2])]

            # Check if any of the products or reactants are outside of the network
            out_of_network = false
            for (z, n) in [reactant; product]
                if !zn_in_network(z, n, net_idx)
                    out_of_network = true
                    break
                end
            end
            # If any of the species are outside of the network, skip this reaction
            if out_of_network
                continue
            end

            push!(reaction_data.alphadecay, AlphaDecay(reactant, product, rate))
        end
    end
    return reaction_data
end

function read_photodissociation!(reverse_reaction_file::String, reaction_data::ReactionData, net_idx::NetworkIndex)
    # Read the reverse reaction nfile
    open(reverse_reaction_file) do file
        lines = readlines(file)
        num_entries::Int64 = parse(Int, lines[1])
        for i in 1:num_entries
            z_r = parse(Int64, lines[5*(i-1) + 2])
            n_r = parse(Int64, lines[5*(i-1) + 3])
            z_ps = [parse(Int64, s) for s in split(lines[5*(i-1) + 4])]
            n_ps = [parse(Int64, s) for s in split(lines[5*(i-1) + 5])]
            q = parse(Float64, lines[5*(i-1) + 6])

            reactant = (z_r, n_r)
            products = [[z_p, n_p] for (z_p, n_p) in zip(z_ps, n_ps)]

            # Check if any of the products or reactants are outside of the network
            out_of_network = false
            for (z, n) in [reactant; products]
                if !zn_in_network(z, n, net_idx)
                    out_of_network = true
                    break
                end
            end
            # If any of the species are outside of the network, skip this reaction
            if out_of_network
                continue
            end

            key = products
            if !haskey(reaction_data.neutroncapture, key)
                key = [products[2], products[1]]
            end
            if !haskey(reaction_data.neutroncapture, key)
                # println("Reverse reaction ($(reactant[1]), $(reactant[2])) doesn't have forward reaction")
                continue
            end
            reaction_data.neutroncapture[key].q = q
        end
    end

    return reaction_data
end

end