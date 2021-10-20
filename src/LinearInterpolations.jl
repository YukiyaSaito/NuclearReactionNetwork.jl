"""
    LinearInterpolations

Provides functionality to perform linear interpolations.

This module is used to interpolate between the data points provided in both the
astrophysical data and the temperature-dependent reaction rates. We do this interpolation
very efficiently by using binary search to find the "knots" (i.e., the points to
interpolate between) and then use a simple linar interpolation method.
"""
module LinearInterpolations

export TrajectoryLerp
export NcapLerp
export get_rate
export get_pfunc
export get_temp_density

"""
    NcapLerp

Data to interpolate between for the rates and the partition functions of neutron capture
reactions.

# Fields
- `temperatures::Vector{Float64}`: The temperatures at which the rates and partition functions are given.
- `rates::Vector{Float64}`: The data points for the rates of reactions for neutron captures to interpolate between.
- `pfuncs::Vector{Float64}`: The data points for the partition function for neutron captures to interpolate between.

See also: [`TrajectoryLerp`](@ref)
"""
struct NcapLerp
    """The temperatures at which the rates and partition functions are given."""
    temperatures::Vector{Float64}
    """The data points for the rates of reactions for neutron captures to interpolate between."""
    rates::Vector{Float64}
    """The data points for the partition function for neutron captures to interpolate between."""
    pfuncs::Vector{Float64}
end

"""
    TrajectoryLerp

Data to interpolate between for the temperatures and the densities of astrophysical
trajectories.

# Fields
- `times::Vector{Float64}`: The times at which the temperatures and partition functions are given.
- `temperatures::Vector{Float64}`: The data points for the temperatures of the trajectory to interpolate between.
- `densities::Vector{Float64}`: The data points for the densities of the trajectory to interpolate between.

See also: [`NcapLerp`](@ref)
"""
struct TrajectoryLerp
    """The times at which the temperatures and partition functions are given."""
    times::Vector{Float64}
    """The data points for the temperatures of the trajectory to interpolate between."""
    temperatures::Vector{Float64}
    """The data points for the densities of the trajectory to interpolate between."""
    densities::Vector{Float64}
end

function bisect_right(vec::Vector{Float64}, val::Float64)::Int
    lo::Int = 1
    hi::Int = length(vec)
    while lo < hi
        mid::Int = div(lo + hi, 2)
        if val < vec[mid]
            hi = mid
        else
            lo = mid + 1
        end
    end
    return lo
end

function get_lerp_idxs(vec::Vector{Float64}, val::Float64)::Tuple{Int, Int}
    right_idx::Int = bisect_right(vec, val)
    if right_idx == length(vec) + 1
        return (right_idx - 1, right_idx - 1)
    elseif right_idx == 1
        return (1, 1)
    end
    return (right_idx - 1, right_idx)
end

function lerp(val1::Float64, val2::Float64, t::Float64)::Float64
    # TODO: Use fma(t, val2 - val1, val1)?
    return t*(val2 - val1) + val1
end

function get_lerp_param(left::Float64, right::Float64, middle::Float64)::Float64
    if left == right
        return 0.0
    end
    return (middle - left) / (right - left)
end

"""
    get_rate(ncap_lerp::NcapLerp, temp::Float64)::Float64

Return the (linearly) interpolated rate of reaction at the temperature `temp`.

See also: [`get_pfunc`](@ref), [`get_temp_density`](@ref)
"""
function get_rate(ncap_lerp::NcapLerp, temp::Float64)::Float64
    idx_left::Int, idx_right::Int = get_lerp_idxs(ncap_lerp.temperatures, temp)

    temp_left::Float64 = ncap_lerp.temperatures[idx_left]
    temp_right::Float64 = ncap_lerp.temperatures[idx_right]
    t::Float64 = get_lerp_param(temp_left, temp_right, temp)

    rate::Float64 = lerp(ncap_lerp.rates[idx_left], ncap_lerp.rates[idx_right], t)
    return rate
end

"""
    get_pfunc(ncap_lerp::NcapLerp, temp::Float64)

Return the (linearly) interpolated partition function at the temperature `temp`.

See also: [`get_rate`](@ref), [`get_temp_density`](@ref)
"""
function get_pfunc(ncap_lerp::NcapLerp, temp::Float64)::Float64
    idx_left::Int, idx_right::Int = get_lerp_idxs(ncap_lerp.temperatures, temp)

    temp_left::Float64 = ncap_lerp.temperatures[idx_left]
    temp_right::Float64 = ncap_lerp.temperatures[idx_right]
    t::Float64 = get_lerp_param(temp_left, temp_right, temp)

    pfunc::Float64 = lerp(ncap_lerp.pfuncs[idx_left], ncap_lerp.pfuncs[idx_right], t)
    return pfunc
end

"""
    get_temp_density(traj_lerp::TrajectoryLerp, time::Float64)

Return the (linearly) interpolated temperature and density at the time `time`.

See also: [`get_rate`](@ref), [`get_pfunc`](@ref)
"""
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