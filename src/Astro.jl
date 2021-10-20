"""
    Astro

Provides functionality to get information about the current temperature and density.
"""
module Astro

using ..LinearInterpolations

export CurrentTrajectory
export get_current_trajectory

"""
    CurrentTrajectory

A convenient structure that holds the current temperature and density.

# Fields
- `temperature::Float64`: The current temperature
- `density::Float64`: The current density
"""
struct CurrentTrajectory
    """The current temperature"""
    temperature::Float64
    """The current density"""
    density::Float64
end

"""
    get_current_trajectory(trajectory::TrajectoryLerp, time::Float64)

Return the temperature and density at time `time` from the data in `trajectory`.
"""
function get_current_trajectory(trajectory::TrajectoryLerp, time::Float64)
    temp, density = get_temp_density(trajectory, time)
    return CurrentTrajectory(temp, density)
end

end