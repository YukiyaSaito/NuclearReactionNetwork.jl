import NuclearReactionNetwork as NRN
using ArgParse
using DelimitedFiles
using JLD2
using ProgressBars

function parse_cli_args()
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "--input", "-i"
            help = "Input file path (.dat file)"
            arg_type = String
            required = true
        "--type", "-t"
            help = "The type of file. Supported type: ncap, probdecay, alphadecay, photodissociation, trajectory, initial-composition, extent"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Output file path (.jld file)"
            arg_type = String
            required = true
    end

    args = parse_args(settings)
    return args
end

function read_ncap(path::String) 
    # Read Fortran formatted input file for temperature dependent reaction rate. 
    ncaps = Vector{NRN.NeutronCapture}()
    open(path) do file
        lines = readlines(file)
        num_entry::Int = parse(Int, lines[1])
        pfunc_flag::Int = parse(Int, lines[2])
        temperature::Vector{Float64} = [parse(Float64, s) for s in split(lines[3])]
        if pfunc_flag == 1
            for i in ProgressBar(1:num_entry)
                z_rs = [parse(Int, s) for s in split(lines[6*(i-1)+4])]
                n_rs = [parse(Int, s) for s in split(lines[6*(i-1)+5])]
                z_ps = [parse(Int, s) for s in split(lines[6*(i-1)+6])]
                n_ps = [parse(Int, s) for s in split(lines[6*(i-1)+7])]
                rates = [parse(Float64, s) for s in split(lines[6*(i-1)+8])]
                pfuncs = [parse(Float64, s) for s in split(lines[6*(i-1)+9])]

                reactants = collect(zip(z_rs, n_rs))
                products = collect(zip(z_ps, n_ps))

                rates_pfuncs_lerp = NRN.LinearInterpolations.NcapLerp(temperature, rates, pfuncs)

                push!(ncaps, NRN.NeutronCapture(reactants, products, rates_pfuncs_lerp, 0.0, missing))
            end
        end
    end
    return ncaps
end

function read_probdecay(path::String) 
    # Read Fortran formatted PRISM input file for probabilistic decay.
    probdecay = Vector{NRN.ProbDecay}()
    open(path) do file
        lines = readlines(file)
        num_entry::Int = parse(Int, lines[1])
        for i in ProgressBar(1:num_entry)
            z_rs::Vector{Int} = [parse(Int, s) for s in split(lines[6*(i-1)+2])]
            n_rs::Vector{Int} = [parse(Int, s) for s in split(lines[6*(i-1)+3])]
            z_ps::Vector{Int} = [parse(Int, s) for s in split(lines[6*(i-1)+4])]
            n_ps::Vector{Int} = [parse(Int, s) for s in split(lines[6*(i-1)+5])]
            rate::Float64 = parse(Float64, lines[6*(i-1)+6])
            average_number::Vector{Float64} = [parse(Float64, s) for s in split(lines[6*(i-1)+7])]

            reactants = collect(zip(z_rs, n_rs))
            products = collect(zip(z_ps, n_ps))

            push!(probdecay, NRN.ProbDecay(reactants, products, rate, average_number))
        end
    end
    return probdecay
end

function read_alphadecay(path::String)
    # Read Fortran formatted input file for temperature dependent reaction rate. 
    alphadecay = Vector{NRN.AlphaDecay}()
    open(path) do file
        lines = readlines(file)
        num_entries::Int = parse(Int, lines[1])
        for i in ProgressBar(1:num_entries)
            z_r = parse(Int, lines[5*(i-1) + 2])
            n_r = parse(Int, lines[5*(i-1) + 3])
            z_p = [parse(Int, s) for s in split(lines[5*(i-1) + 4])]
            n_p = [parse(Int, s) for s in split(lines[5*(i-1) + 5])]
            rate = parse(Float64, lines[5*(i-1) + 6])

            reactant = (z_r, n_r)
            product = [(z_p[1], n_p[1]), (z_p[2], n_p[2])]

            push!(alphadecay, NRN.AlphaDecay(reactant, product, rate))
        end
    end
    return alphadecay
end

function read_photodissociation(path::String)
    # Read the reverse reaction file
    photodissociation_dict = Dict{Tuple{Int, Int}, NRN.Photodissociation}()
    open(path) do file
        lines = readlines(file)
        num_entries::Int = parse(Int, lines[1])
        for i in ProgressBar(1:num_entries)
            z_r = parse(Int, lines[5*(i-1) + 2])
            n_r = parse(Int, lines[5*(i-1) + 3])
            q = parse(Float64, lines[5*(i-1) + 6])

            reactant = (z_r, n_r)

            photodissociation_dict[reactant] = NRN.Photodissociation(q)
        end
    end
    return photodissociation_dict
end

function read_trajectory(path::String)
    trajectory_matrix::Matrix{Float64} = readdlm(path, skipstart=1)
    times, temperatures, densities = eachcol(trajectory_matrix)
    return NRN.LinearInterpolations.TrajectoryLerp(times, temperatures, densities)
end

function read_init_comp(path::String)
    data::Matrix{Float64} = readdlm(path)
    return [(Int(round(z)), Int(round(a)), y) for (z, a, y) in eachrow(data)]
end

function read_extent(path::String)
    raw_boundary::Matrix{Int} = readdlm(path, ' ', Int)
    boundary = NRN.NetworkBoundary{Matrix{Int}}(raw_boundary)
    return NRN.NetworkIndex(boundary)
end

function main()
    args = parse_cli_args()
    input_fpath = args["input"]
    output_fpath = args["output"]
    type = args["type"]

    obj = missing
    if type == "ncap"
        obj = read_ncap(input_fpath)
    elseif type == "probdecay"
        obj = read_probdecay(input_fpath)
    elseif type == "alphadecay"
        obj = read_alphadecay(input_fpath)
    elseif type == "photodissociation"
        obj = read_photodissociation(input_fpath)
    elseif type == "trajectory"
        obj = read_trajectory(input_fpath)
    elseif type == "initial-composition"
        obj = read_init_comp(input_fpath)
    elseif type == "extent"
        obj = read_extent(input_fpath)
    end
    if !ismissing(obj)
        save_object(output_fpath, obj)
    else
        error("Uknown type: $(type)")
    end
end

main()