println("Importing modules...")
using CairoMakie
using Makie.GeometryBasics
using CSV
using DelimitedFiles
using DataFrames
using ProgressBars

measured_path = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/Misc/species_data.csv";
data_path = "/Users/pvirally/Dropbox/Waterloo/Co-op/TRIUMF/prism1.5.0/input/output/dynamical-simple_fission-YTime.txt"
out_path = "/Users/pvirally/Desktop/wind-all.mp4"

struct Epoch
    timestep::Int
    time::Float64
    temp::Float64
    density::Float64
    zs::Vector{Int}
    ns::Vector{Int}
    abundances::Vector{Float64}
end

function read_epochs(path)
    epochs = Vector{Epoch}()
    open(path) do file
        lines = readlines(file)[2:end]

        timestep = nothing
        time = nothing
        temp = nothing
        density = nothing
        zs = Vector{Int}()
        ns = Vector{Int}()
        abundances = Vector{Float64}()

        pb_lines = ProgressBar(lines)
        set_description(pb_lines, "Reading time step data")
        for line in pb_lines
            line = strip(line)
            if line == ""
                push!(epochs, Epoch(timestep, time, temp, density, zs, ns, abundances))
                zs = Vector{Int}()
                ns = Vector{Int}()
                abundances = Vector{Float64}()
                continue
            end
            if startswith(line, "timestep")
                timestep, time, temp, density = split(line)[2:2:end]
                timestep = parse(Int, timestep)
                time = parse(Float64, time)
                temp = parse(Float64, temp)
                density = parse(Float64, density)
                continue
            end

            z, n, abundance = split(line)
            z = parse(Int, z)
            n = parse(Int, n)
            abundance = parse(Float64, abundance)

            push!(zs, z)
            push!(ns, n)
            push!(abundances, abundance)
        end
    end
    return epochs
end

function remap(x, x0, x1, y0, y1)
    return y0 .+ y1 .* (x .- x0)/(x1 .- x0)
end

function zn_to_rect(z, n)
    return Rect(n-0.4, z-0.4, 0.8, 0.8)
end

function plot_measured(df)
    measured_rects = Vector{Rect2D}()
    stable_rects = Vector{Rect2D}()
    for (z, n, half_life) in eachrow(df)
        if ismissing(z) || ismissing(n)
            continue
        end
        if !ismissing(half_life) && half_life == "STABLE"
            push!(stable_rects, zn_to_rect(z, n))
        else
            push!(measured_rects, zn_to_rect(z, n))
        end
    end
end

println("Reading stable isotope file...")
measured_df = DataFrame(CSV.File(measured_path));

stable_rects = Vector{Rect2D}()
measured_rects = Vector{Rect2D}()
for (z, n, half_life) in eachrow(measured_df)
    if ismissing(z) || ismissing(n)
        continue
    end
    if !ismissing(half_life) && half_life == "STABLE"
        push!(stable_rects, zn_to_rect(z, n))
    else
        push!(measured_rects, zn_to_rect(z, n))
    end
end

epochs = read_epochs(data_path)
times = [epoch.time for epoch in epochs]
densities = [epoch.density for epoch in epochs]
temps = [epoch.temp for epoch in epochs]

println("Setting up animation...")
measured_fig = Figure()
ax = Axis(measured_fig[1, 1], xlabel="Neutron Number (N)", ylabel="Proton Number (Z)")

temp_ax = Axis(measured_fig[1, 1], width=Relative(0.3), height=Relative(0.3), halign=0.85, valign=0.15, yticklabelcolor=:blue, ytickcolor=:blue, backgroundcolor=:lightgray, xscale=log10, yscale=log10, xlabel="Time [s]", ylabel="Tmperature [GK]", xlabelsize=12, ylabelsize=12)
density_ax = Axis(measured_fig[1, 1], width=Relative(0.3), height=Relative(0.3), halign=0.85, valign=0.15, yticklabelcolor=:red, ytickcolor=:red, backgroundcolor=:transparent, yaxisposition=:right, xscale=log10, yscale=log10, ylabel="Density [g/cm^3]", xlabelsize=12, ylabelsize=12)
hidexdecorations!(density_ax)
hidespines!(density_ax)
linkxaxes!(temp_ax, density_ax)

Colorbar(measured_fig[1, 2], limits=(0.0, 1.0), colormap=:viridis, label="Abundance")

poly!(ax, stable_rects, color=:black)
poly!(ax, measured_rects, color=:grey)

temp_points = Node(Point2[])
density_points = Node(Point2[])
lines!(temp_ax, temp_points, color=:blue)
lines!(density_ax, density_points, color=:red)
limits!(temp_ax, minimum(times), maximum(times), minimum(temps), maximum(temps))
limits!(density_ax, minimum(times), maximum(times), minimum(densities), maximum(densities))

to_remove = []
pb_epochs = ProgressBar(epochs)
set_description(pb_epochs, "Creating animation")
record(measured_fig, out_path, pb_epochs; framerate=30) do epoch
    global ax, traj_ax, to_remove, temp_points, density_points
    for plot in to_remove
        delete!(ax, plot)
    end
    to_remove = []
    
    temp_points[] = push!(temp_points[], Point2(epoch.time, epoch.temp))
    density_points[] = push!(density_points[], Point2(epoch.time, epoch.density))
    
    push!(to_remove, text!(ax, "Time: $(epoch.time) [s]\nTemperature $(epoch.temp) [GK]\nDensity: $(epoch.density) [g/cm^3]", position=(0, 100)))
    data_rects = map(zn_to_rect, epoch.zs, epoch.ns)
    colors = remap.(log.(epoch.abundances), log(1e-15), log(1.0), 0.0, 1.0)
    push!(to_remove, poly!(ax, data_rects, color=colors, colormap=:viridis))
end
