module Astro
using DelimitedFiles
export read_trajectory


struct Trajectory
    time::Vector{Float64}
    temperature::Vector{Float64}
    density::Vector{Float64}
end

mutable struct CurrentTrajectory
    temperature::Float64
    density::Float64
end

# function initialize_current_trajectory()

function read_trajectory(path::String)
    trajectory_matrix::Matrix{Float64} = readdlm(path,skipstart=1)
    return Trajectory(trajectory_matrix[:,1], trajectory_matrix[:,2], trajectory_matrix[:,3])
end


end