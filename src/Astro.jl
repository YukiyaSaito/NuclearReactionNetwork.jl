module Astro
using DelimitedFiles
using Interpolations

export Trajectory
export CurrentTrajectory
export read_trajectory
export get_current_trajectory

struct Trajectory
    temperatures::Interpolations.Extrapolation
    densities::Interpolations.Extrapolation
    function Trajectory(temperatures::Interpolations.Extrapolation, densities::Interpolations.Extrapolation)
        return new(temperatures, densities)
    end
end

mutable struct CurrentTrajectory
    temperature::Float64
    density::Float64
    function CurrentTrajectory(temperature::Float64, density::Float64)
        return new(temperature, density)
    end
end

# TODO: Should this be a constructor for a Trajectory object?
function read_trajectory(path::String)
    trajectory_matrix::Matrix{Float64} = readdlm(path,skipstart=1)
    times, temperatures, densities = eachcol(trajectory_matrix)
    temperatures_lerp = LinearInterpolation(times, temperatures, extrapolation_bc=Flat())
    densities_lerp = LinearInterpolation(times, densities, extrapolation_bc=Flat())
    return Trajectory(temperatures_lerp, densities_lerp)
end

# TODO: Should this be a constructor for a CurrentTrajectory object?
function get_current_trajectory(trajectory::Trajectory, time::Float64)
    return CurrentTrajectory(trajectory.temperatures(time), trajectory.densities(time))
end

end