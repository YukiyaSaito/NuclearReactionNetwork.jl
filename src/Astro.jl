module Astro
using Interpolations

export Trajectory
export CurrentTrajectory
export get_current_trajectory

struct Trajectory
    temperatures::Interpolations.Extrapolation
    densities::Interpolations.Extrapolation
end

struct CurrentTrajectory
    temperature::Float64
    density::Float64
end

# TODO: Should this be a constructor for a CurrentTrajectory object?
function get_current_trajectory(trajectory::Trajectory, time::Float64)
    return CurrentTrajectory(trajectory.temperatures(time), trajectory.densities(time))
end

end