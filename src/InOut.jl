module InOut
using DelimitedFiles
using ..ReactionTypes
using ..Network

export Result
export dump_result

# TODO: Should this be a matrix rather than a series of vectors of the same size?
struct Result
    proton_nums::Vector{Int64}
    neutron_nums::Vector{Int64}
    abundance::Vector{Float64} # TODO: Add more things that we can output, like PRISM does
    function Result(proton_nums::Vector{Int64}, neutron_nums::Vector{Int64}, abundance::Vector{Float64})
        return new(proton_nums, neutron_nums, abundance)
    end
end

function Result(abundace::Vector{Float64}, net_idx::NetworkIndex)
    boundary = net_idx.networkboundary.matrix

    zs = zeros(Int64, length(abundace))
    ns = zeros(Int64, length(abundace))

    index = 1
    for (z, n_low, n_high) in eachrow(boundary)
        dn = n_high - n_low
        zs[index:index+dn] = repeat([z], dn+1)
        ns[index:index+dn] = range(n_low, n_high, length=dn+1)
        index += dn + 1
    end
    return Result(zs, ns, abundace)
end

function dump_result(result::Result, path::String)
    open(path, "w") do out_file
        writedlm(out_file, [result.proton_nums result.neutron_nums.+result.proton_nums result.abundance])
    end
end

# begin
# 	open("moller2003.dat") do file
# 	index = 1
# 		for ln in eachline(file)
# 			println([parse(Int, ss) for ss in split(ln)])
			
# 			index += 1
# 		end
# 	end
# end


# using PyCall
# ng = fort.FortranFile("input/nuclear/Reaclib/moller2003.bin","r")
# begin
# 	ng.read_ints()
# 	typeof(ng.read_ints())
# end

end