module InOut
using DelimitedFiles
using ..ReactionTypes
using ..Network
using ..Astro
using ..Network
using ..NetworkDatas

export Result
export dump_result
export dump_iteration

# TODO: Should this be a matrix rather than a series of vectors of the same size (we would lose type inhomogeneity)?
struct Result
    proton_nums::Vector{Int64}
    neutron_nums::Vector{Int64}
    abundance::Vector{Float64} # TODO: Add more things that we can output, like PRISM does
    function Result(proton_nums::Vector{Int64}, neutron_nums::Vector{Int64}, abundance::Vector{Float64})
        return new(proton_nums, neutron_nums, abundance)
    end
end

function Result(nd::NetworkData)
    boundary = nd.net_idx.networkboundary.matrix

    zs = zeros(Int64, length(nd.abundance))
    ns = zeros(Int64, length(nd.abundance))

    index = 1
    for (z, n_low, n_high) in eachrow(boundary)
        dn = n_high - n_low
        zs[index:index+dn] = repeat([z], dn+1)
        ns[index:index+dn] = range(n_low, n_high, length=dn+1)
        index += dn + 1
    end
    return Result(zs, ns, nd.abundance)
end

function dump_result(result::Result, path::String)
    open(path, "w") do out_file
        # TODO: This should really write [Int64 Int64 Float64], but right now it gets promoted to [Float64 Float64 Float64]
        writedlm(out_file, [result.proton_nums result.neutron_nums.+result.proton_nums result.abundance])
    end
end

function dump_iteration(result::Result, trajectory::Trajectory, time::Time, iteration::Int64, path::String)
    curr_traj = get_current_trajectory(trajectory, time.current)
    mode = iteration == 1 ? "w" : "a"
    open(path, mode) do out_file
        write(out_file, "$(iteration)\t$(time.current)\t$(curr_traj.temperature)\t$(curr_traj.density)\n")
        last_idx = findall(x->!iszero(x), result.abundance)[end]
        Zs = result.proton_nums[1:last_idx]
        As = result.proton_nums[1:last_idx] .+ result.neutron_nums[1:last_idx]
        ys = result.abundance[1:last_idx]
        writedlm(out_file, [Zs As ys])
        write(out_file, "\n")
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