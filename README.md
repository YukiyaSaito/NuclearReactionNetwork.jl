# NuclearReactionNetwork.jl

[![Build Status](https://travis-ci.com/YukiyaSaito/NuclearReactionNetwork.jl.svg?branch=master)](https://travis-ci.com/YukiyaSaito/NuclearReactionNetwork.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/YukiyaSaito/NuclearReactionNetwork.jl?svg=true)](https://ci.appveyor.com/project/YukiyaSaito/NuclearReactionNetwork-jl)
[![Coverage](https://codecov.io/gh/YukiyaSaito/NuclearReactionNetwork.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/YukiyaSaito/NuclearReactionNetwork.jl)
[![Coverage](https://coveralls.io/repos/github/YukiyaSaito/NuclearReactionNetwork.jl/badge.svg?branch=master)](https://coveralls.io/github/YukiyaSaito/NuclearReactionNetwork.jl?branch=master)

## Running a Calculation

`NuclearReactionNetwork.jl` has two modes of operation. You can have it run a single calculation, or multiple calculations. If you want to run multiple calculations, you have to specify how each calculation varies from run to run. Currently, the only way of doing this is by specifying "rate modulations" (multiplicative factors that alter the reaction rates of specific types of reactions). You specify these rate modulations in the control file. Here is an example:

```json
"rate_modulations": [
    { "rxn_type": "probdecay", "id": "beta", "path": "/path/to/mod_factors", "description": "Multiplicative factors for the rates of beta decays" },
    { "rxn_type": "alphadecay", "id": "alpha", "path" : "/path/to/mod_factors", "description": "Multiplicative factors for the rates of alpha decays" }
]
```

The actual `mod_factors` file from the above contains one floating point value per row. Below is an example:
```
1e-2
1e-1
1e0
1e1
1e2
```

`NuclearReactioinNetwork.jl` has the ability to run each calculation on separate threads. See the [Multi-Threaded Calculation](#multi-threaded-calculation) section.

### Single-Threaded Calculation

To actually run the calculation, simply open julia, activate the environment, and run the `test/runtests.jl` file. Here is an example session:
```
bash-5.1$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.4 (2021-11-19)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.6) pkg> activate .
  Activating environment at `~/Dropbox/Waterloo/Co-op/TRIUMF/NuclearReactionNetwork.jl/Project.toml`

julia> include("test/runtests.jl")
[ Info: Precompiling NuclearReactionNetwork [a7ba2b4c-5c8a-4117-add0-163f5ea5e167]
```

Notice that julia was run with no command line arguments. To activate the environment, from the julia REPL, type `]activate .`. The `]` character puts you in the package mode and allows you to add and remove packages. Running `activate .` inside the package mode will ensure the file structure of the program is properly loaded in. Then, you should press backspace to get back to the regular julia REPL and run `include("test/runtests.jl")` to actually run the code.


### Multi-Threaded Calculation

To run the calculation in multi-threaded mode, you can follow the same steps as in the [Single-Threaded Calculation](#single-threaded-calculation) section, but by passing the `-t 4` (or whatever number of threads you want to use) flag to julia. Julia runs in single-threaded mode by default, but by running `julia -t <n>` you can specify the number of threads julia will use. Note that the number of threads need not match the number of modulation factors used by the program.

```
bash-5.1$ julia -t 4
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.4 (2021-11-19)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.6) pkg> activate .
  Activating environment at `~/Dropbox/Waterloo/Co-op/TRIUMF/NuclearReactionNetwork.jl/Project.toml`

julia> include("test/runtests.jl")
[ Info: Precompiling NuclearReactionNetwork [a7ba2b4c-5c8a-4117-add0-163f5ea5e167]
```

Notice that julia was run with the flag `-t 4`. To activate the environment, from the julia REPL, type `]activate .`. The `]` character puts you in the package mode and allows you to add and remove packages. Running `activate .` inside the package mode will ensure the file structure of the program is properly loaded in. Then, you should press backspace to get back to the regular julia REPL and run `include("test/runtests.jl")` to actually run the code.

## Creating `.jld` files

There are two acceptable file types as inputs to `NuclearReactionNetwork.jl`. The first is simple ASCII data. This format is human readable, but slower to read in for the program. The file extension of this ASCII data must be `.dat`, as `NuclearReactionNetwork.jl` differentiates between file types exclusively by the extension of the file name.

The second format is an HDF5-like format that is created with the library [JLD2](https://github.com/JuliaIO/JLD2.jl). This format is not human readable, but much faster for the program to read in. The file extension of this data format is normally `.jld`.

To create a `.jld` file from a `.dat` file, you can use the script `dat2jld.jl` located in the `scripts` directory. To use this script you must specify the input `.dat` file, the location of the output `.jld` file, and the type of file the input file is.

Below is an example of how to use the program:
```sh
$ julia dat2jld.jl scripts/dat2jld.jl -i /path/to/input/my_decay.dat -o /path/to/output/my_decay.jld -t decay
```

The `-i` or `--input` flag is used to specify the input file path. The `-o` or `--output` flag is used to specify the output file path. The `-t` or `--type` flag is used to specify the type of file the input file is. The supported types are `decay`, `rxn`, `probrxn`, `ncap`, `probdecay`, `alphadecay`, `photodissociation`, `trajectory`, `initial-composition`, `extent`.

For more help, you can always run `scripts/dat2jld.jl --help`.
```sh
$ julia scripts/dat2jld.jl --help
usage: dat2jld.jl -i INPUT -t TYPE -o OUTPUT [-h]

optional arguments:
  -i, --input INPUT    Input file path (.dat file)
  -t, --type TYPE      The type of file. Supported type: decay, rxn,
                       probrxn, ncap, probdecay, alphadecay,
                       photodissociation, trajectory,
                       initial-composition, extent
  -o, --output OUTPUT  Output file path (.jld file)
  -h, --help           show this help message and exit
```