module LinearInterpolations

export TrajectoryLerp
export NcapLerp
export get_rate
export get_pfunc
export get_temp_density

struct NcapLerp
    temperatures::Vector{Float64}
    rates::Vector{Float64}
    pfuncs::Vector{Float64}
end

struct TrajectoryLerp
    times::Vector{Float64}
    temperatures::Vector{Float64}
    densities::Vector{Float64}
end

function get_lerp_idxs(vec::Vector{Float64}, val::Float64)::Tuple{Int, Int}
    right_idx::Int = searchsortedfirst(vec, val)
    if right_idx == length(vec) + 1
        return (right_idx - 1, right_idx - 1)
    end
    if right_idx == 1
        return (1, 1)
    end
    return (right_idx - 1, right_idx)
end

function lerp(val1::Float64, val2::Float64, t::Float64)::Float64
    return val1 + t*(val2 - val1)
end

function get_lerp_param(left::Float64, right::Float64, middle::Float64)::Float64
    if left == right
        return 0.0
    end
    return (middle - left) / (right - left)
end

function get_rate(ncap_lerp::NcapLerp, temp::Float64)::Float64
    idx_left::Int, idx_right::Int = get_lerp_idxs(ncap_lerp.temperatures, temp)

    temp_left::Float64 = ncap_lerp.temperatures[idx_left]
    temp_right::Float64 = ncap_lerp.temperatures[idx_right]
    t::Float64 = get_lerp_param(temp_left, temp_right, temp)

    rate::Float64 = lerp(ncap_lerp.rates[idx_left], ncap_lerp.rates[idx_right], t)
    return rate
end

function get_pfunc(ncap_lerp::NcapLerp, temp::Float64)::Float64
    idx_left::Int, idx_right::Int = get_lerp_idxs(ncap_lerp.temperatures, temp)

    temp_left::Float64 = ncap_lerp.temperatures[idx_left]
    temp_right::Float64 = ncap_lerp.temperatures[idx_right]
    t::Float64 = get_lerp_param(temp_left, temp_right, temp)

    pfunc::Float64 = lerp(ncap_lerp.pfuncs[idx_left], ncap_lerp.pfuncs[idx_right], t)
    return pfunc
end

function get_temp_density(traj_lerp::TrajectoryLerp, time::Float64)::Tuple{Float64, Float64}
    idx_left::Int, idx_right::Int = get_lerp_idxs(traj_lerp.times, time)

    time_left::Float64 = traj_lerp.times[idx_left]
    time_right::Float64 = traj_lerp.times[idx_right]
    t::Float64 = get_lerp_param(time_left, time_right, time)

    temp::Float64 = lerp(traj_lerp.temperatures[idx_left], traj_lerp.temperatures[idx_right], t)
    density::Float64 = lerp(traj_lerp.densities[idx_left], traj_lerp.densities[idx_right], t)
    return (temp, density)
end

end