module Astro

using ..LinearInterpolations

export CurrentTrajectory
export get_current_trajectory

struct CurrentTrajectory
    temperature::Float64
    density::Float64
end

function get_current_trajectory(trajectory::TrajectoryLerp, time::Float64)
    temp, density = get_temp_density(trajectory, time)
    return CurrentTrajectory(temp, density)
end

end